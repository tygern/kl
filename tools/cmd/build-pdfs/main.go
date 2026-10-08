// Command build-pdfs exports the project's LaTeX manuscripts under results/
// to PDF.
//
// It ports research/build_pdfs.py to Go (standard library only; the only
// external program is the installed latexmk, run with shell escape disabled
// and no packages installed). Only results/*.tex is built; the archived
// research notes under research/archive/notes/ and the frozen manuscript copy
// inside the supplement are not project deliverables and are deliberately
// skipped.
//
// Invocation, from the repository root:
//
//	go run -C tools ./cmd/build-pdfs
//
// `go run -C tools` runs the program with tools/ as the working directory, so
// the repository root is found by walking upward from the working directory
// until a directory containing results/ and tools/go.mod is reached (the -root
// flag overrides this). For each results/<stem>.tex it runs latexmk from the
// repository root with repository-relative arguments, builds in
// tmp/pdfs/build/<stem>, copies the PDF to output/pdf/<stem>.pdf, writes the
// combined latexmk output to results/pdf-build.log and writes
// results/pdf-build-manifest.json (same schema as the Python original: source,
// pdf, source_sha256, pdf_sha256, bytes, command, cwd). Every recorded path is
// repository-relative. The pdf_sha256 values differ between builds because
// pdfTeX embeds creation timestamps and an identifier.
//
// Imports: standard library only; no other package of this repository.
package main

import (
	"bytes"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
)

// record is one manifest entry; field order is the key order of the original.
type record struct {
	Source       string   `json:"source"`
	PDF          string   `json:"pdf"`
	SourceSHA256 string   `json:"source_sha256"`
	PDFSHA256    string   `json:"pdf_sha256"`
	Bytes        int      `json:"bytes"`
	Command      []string `json:"command"`
	CWD          string   `json:"cwd"`
}

func main() {
	rootFlag := flag.String("root", "", "repository root (default: found by walking up from the working directory)")
	flag.Parse()
	if flag.NArg() != 0 {
		fatalf("unexpected arguments: %v", flag.Args())
	}
	root, err := findRoot(*rootFlag)
	if err != nil {
		fatalf("%v", err)
	}
	pdfs, err := run(root)
	if err != nil {
		fatalf("%v", err)
	}
	for _, p := range pdfs {
		fmt.Println(p)
	}
}

func fatalf(format string, args ...any) {
	fmt.Fprintf(os.Stderr, format+"\n", args...)
	os.Exit(1)
}

// findRoot returns the repository root: the explicit one if given, otherwise
// the nearest ancestor of the working directory that holds results/ and
// tools/go.mod.
func findRoot(explicit string) (string, error) {
	start := explicit
	if start == "" {
		wd, err := os.Getwd()
		if err != nil {
			return "", err
		}
		start = wd
	}
	start, err := filepath.Abs(start)
	if err != nil {
		return "", err
	}
	if resolved, err := filepath.EvalSymlinks(start); err == nil {
		start = resolved
	}
	if explicit != "" {
		if !isRoot(start) {
			return "", fmt.Errorf("%s is not a repository root (needs results/ and tools/go.mod)", explicit)
		}
		return start, nil
	}
	for dir := start; ; {
		if isRoot(dir) {
			return dir, nil
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			return "", errors.New("repository root not found: run from the repository root or below it")
		}
		dir = parent
	}
}

func isRoot(dir string) bool {
	if st, err := os.Stat(filepath.Join(dir, "results")); err != nil || !st.IsDir() {
		return false
	}
	st, err := os.Stat(filepath.Join(dir, "tools", "go.mod"))
	return err == nil && st.Mode().IsRegular()
}

// rel returns the repository-relative slash-separated path of p.
func rel(root, p string) (string, error) {
	r, err := filepath.Rel(root, p)
	if err != nil {
		return "", err
	}
	if r == ".." || strings.HasPrefix(r, ".."+string(filepath.Separator)) {
		return "", fmt.Errorf("%s is outside the repository", p)
	}
	return filepath.ToSlash(r), nil
}

// normalizeNewlines mimics Python's universal-newline decoding of the child's
// output: "\r\n" and a lone "\r" both become "\n".
func normalizeNewlines(b []byte) []byte {
	b = bytes.ReplaceAll(b, []byte("\r\n"), []byte("\n"))
	return bytes.ReplaceAll(b, []byte("\r"), []byte("\n"))
}

func sha256Hex(b []byte) string {
	sum := sha256.Sum256(b)
	return hex.EncodeToString(sum[:])
}

// copyFile copies src to dst preserving the permission bits and the
// modification time (shutil.copy2 semantics for regular files).
func copyFile(src, dst string) error {
	in, err := os.Open(src)
	if err != nil {
		return err
	}
	defer in.Close()
	st, err := in.Stat()
	if err != nil {
		return err
	}
	out, err := os.OpenFile(dst, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, st.Mode().Perm())
	if err != nil {
		return err
	}
	if _, err := io.Copy(out, in); err != nil {
		out.Close()
		return err
	}
	if err := out.Close(); err != nil {
		return err
	}
	if err := os.Chmod(dst, st.Mode().Perm()); err != nil {
		return err
	}
	return os.Chtimes(dst, st.ModTime(), st.ModTime())
}

// marshalManifest renders the manifest exactly as json.dumps(records,
// indent=2) + "\n" does for this schema.
func marshalManifest(records []record) ([]byte, error) {
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	enc.SetIndent("", "  ")
	if err := enc.Encode(records); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil // Encode appends the trailing newline
}

// run builds every results/*.tex and returns the repository-relative PDF paths.
func run(root string) ([]string, error) {
	sources, err := filepath.Glob(filepath.Join(root, "results", "*.tex"))
	if err != nil {
		return nil, err
	}
	sort.Strings(sources)
	if len(sources) == 0 {
		return nil, errors.New("No LaTeX sources found under results/.")
	}
	stems := map[string]bool{}
	for _, s := range sources {
		stems[strings.TrimSuffix(filepath.Base(s), ".tex")] = true
	}
	if len(stems) != len(sources) {
		return nil, errors.New("Duplicate TeX basenames need distinct output names.")
	}
	compiler, err := exec.LookPath("latexmk")
	if err != nil {
		return nil, errors.New("The installed latexmk executable is unavailable.")
	}
	out := filepath.Join(root, "output", "pdf")
	if err := os.MkdirAll(out, 0o755); err != nil {
		return nil, err
	}

	var records []record
	var logs [][]byte
	for _, source := range sources {
		stem := strings.TrimSuffix(filepath.Base(source), ".tex")
		build := filepath.Join(root, "tmp", "pdfs", "build", stem)
		if err := os.MkdirAll(build, 0o755); err != nil {
			return nil, err
		}
		relBuild, err := rel(root, build)
		if err != nil {
			return nil, err
		}
		relSource, err := rel(root, source)
		if err != nil {
			return nil, err
		}
		// Run from the repository root with relative paths, so the recorded
		// command is exactly the command that ran and contains no machine paths.
		command := []string{"latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error",
			"-file-line-error", "-no-shell-escape", "-outdir=" + relBuild, relSource}
		cmd := exec.Command(compiler, command[1:]...)
		cmd.Dir = root
		output, runErr := cmd.CombinedOutput()
		var exitErr *exec.ExitError
		if runErr != nil && !errors.As(runErr, &exitErr) {
			return nil, fmt.Errorf("cannot run latexmk: %v", runErr)
		}
		logs = append(logs, normalizeNewlines(output))
		logData := bytes.Join(logs, []byte("\n"))
		if err := os.WriteFile(filepath.Join(root, "results", "pdf-build.log"), logData, 0o644); err != nil {
			return nil, err
		}
		if runErr != nil {
			return nil, fmt.Errorf("Build failed: %s; see results/pdf-build.log", relSource)
		}
		generated := filepath.Join(build, stem+".pdf")
		data, err := os.ReadFile(generated)
		if err != nil {
			return nil, err
		}
		if !bytes.HasPrefix(data, []byte("%PDF-")) {
			return nil, fmt.Errorf("%s does not start with %%PDF-", generated)
		}
		target := filepath.Join(out, stem+".pdf")
		if err := copyFile(generated, target); err != nil {
			return nil, err
		}
		srcBytes, err := os.ReadFile(source)
		if err != nil {
			return nil, err
		}
		relTarget, err := rel(root, target)
		if err != nil {
			return nil, err
		}
		records = append(records, record{
			Source:       relSource,
			PDF:          relTarget,
			SourceSHA256: sha256Hex(srcBytes),
			PDFSHA256:    sha256Hex(data),
			Bytes:        len(data),
			Command:      command,
			CWD:          "repository root",
		})
	}
	manifest, err := marshalManifest(records)
	if err != nil {
		return nil, err
	}
	if err := os.WriteFile(filepath.Join(root, "results", "pdf-build-manifest.json"), manifest, 0o644); err != nil {
		return nil, err
	}
	pdfs := make([]string, len(records))
	for i, r := range records {
		pdfs[i] = r.PDF
	}
	return pdfs, nil
}
