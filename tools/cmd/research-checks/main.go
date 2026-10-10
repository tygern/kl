// Command research-checks regenerates every research-note output that the
// tools module certifies and compares each regenerated file with the
// committed one. It is the counterpart, for the exploratory certificates
// outside the proof supplement, of the supplement's proof runner.
//
// Run from the repository root (or anywhere below it):
//
//	go run -C tools ./cmd/research-checks
//
// The repository root is the nearest ancestor of the working directory that
// contains results/ and tools/go.mod (-root overrides). Nothing in the
// repository is modified: every program writes into a scratch copy of the
// affected directories under tmp/research-checks/ (tmp/ is gitignored and the
// directory is recreated on every run), and the regenerated files are compared
// with the committed ones. JSON files are compared semantically (numbers by
// value, objects regardless of key order, lists in order; no key is ignored,
// none of these certificates records timing or paths); text files are compared
// byte for byte. The programs are built with `go build` into the scratch
// bin/ directory; their stdout and stderr are kept under logs/.
//
// Steps, in order (one row per compared file):
//
//	sources-manifest    sources/manifest.json, hashing the locally archived
//	                    PDFs; the whole sources/ directory is gitignored, so the
//	                    comparison is with the local copy and the step is skipped
//	                    with a message when sources/ is absent
//	i6-shape            results/i6-shape-certificate.json and
//	                    results/i6-independent-subword.json
//	exceptional-cosets  research/exceptional-cosets.txt (the program's stdout)
//	star-family         research/broad_cells/star_family_checks.json and
//	                    research/broad_cells/affine_d4_kl_certificate.json
//	affine-d-checks     research/broad_affine/unfolding_obstruction.json and
//	                    research/broad_affine/gern_coatoms.json
//	catalogue-eligible  research/en_families/e9_full_support_eligible.json;
//	                    its input, the complete E9 FC catalogue
//	                    research/en_independent/e9-fc.json, is a gitignored
//	                    output of the supplement's fc-catalogue, so it is
//	                    generated first (into the scratch tree) with
//	                    `fc-catalogue -rank 9`. Rank 9 is the only rank with a
//	                    committed eligible file: the committed E8 conjugate
//	                    search has 20,000 trials rather than 50,000, rank 10 has
//	                    no full-support terminal, and ranks 11 and 12 have no FC
//	                    catalogue.
//
// A table of the results is printed and the exit status is non-zero when any
// file differs, a program fails, or a comparison cannot be made.
//
// Imports: standard library only; no other package of this repository.
package main

import (
	"errors"
	"flag"
	"fmt"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sort"
	"strings"
	"time"
)

// row is one line of the printed table.
type row struct {
	step, output, result, detail string
	seconds                      float64
}

type checker struct {
	root, scratch, bins, logs, goTool string
	rows                              []row
	failed                            bool
}

func fatalf(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "research-checks: "+format+"\n", args...)
	os.Exit(1)
}

func isRoot(dir string) bool {
	st, err := os.Stat(filepath.Join(dir, "results"))
	if err != nil || !st.IsDir() {
		return false
	}
	st, err = os.Stat(filepath.Join(dir, "tools", "go.mod"))
	return err == nil && st.Mode().IsRegular()
}

// findRoot returns the explicit root, or walks upward from the working
// directory to the first directory holding results/ and tools/go.mod.
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
			return "", errors.New("repository root not found (no ancestor has results/ and tools/go.mod); use -root")
		}
		dir = parent
	}
}

func findGo() string {
	if p, err := exec.LookPath("go"); err == nil {
		return p
	}
	if goroot := os.Getenv("GOROOT"); goroot != "" {
		p := filepath.Join(goroot, "bin", "go")
		if _, err := os.Stat(p); err == nil {
			return p
		}
	}
	fatalf("the go command was not found on PATH")
	return ""
}

func copyFile(src, dst string) error {
	data, err := os.ReadFile(src)
	if err != nil {
		return err
	}
	if err := os.MkdirAll(filepath.Dir(dst), 0o755); err != nil {
		return err
	}
	return os.WriteFile(dst, data, 0o644)
}

func seconds(d time.Duration) float64 { return math.Round(d.Seconds()*100) / 100 }

// run executes a program with the given working directory, keeping its
// stdout and stderr under logs/<name>.*.log, and returns stdout and the error.
func (c *checker) run(name string, dir string, env []string, command ...string) ([]byte, error) {
	cmd := exec.Command(command[0], command[1:]...)
	cmd.Dir = dir
	if env != nil {
		cmd.Env = append(os.Environ(), env...)
	}
	errFile, err := os.Create(filepath.Join(c.logs, name+".stderr.log"))
	if err != nil {
		return nil, err
	}
	defer errFile.Close()
	cmd.Stderr = errFile
	out, err := cmd.Output()
	if werr := os.WriteFile(filepath.Join(c.logs, name+".stdout.log"), out, 0o644); werr != nil && err == nil {
		err = werr
	}
	if err != nil {
		return out, fmt.Errorf("%s failed (%v); see %s", name, err, filepath.ToSlash(filepath.Join("tmp", "research-checks", "logs", name+".stderr.log")))
	}
	return out, nil
}

// build compiles one command of a module into the scratch bin/ directory.
func (c *checker) build(module, name string) string {
	exe := filepath.Join(c.bins, name)
	if runtime.GOOS == "windows" {
		exe += ".exe"
	}
	// With GOWORK=off and empty GOFLAGS, builds ignore the enclosing workspace
	// and the caller's flags.
	if _, err := c.run("build-"+name, filepath.Join(c.root, module), []string{"GOWORK=off", "GOFLAGS="},
		c.goTool, "build", "-o", exe, "./cmd/"+name); err != nil {
		fatalf("%v", err)
	}
	return exe
}

func (c *checker) record(step, output, result, detail string, started time.Time) {
	c.rows = append(c.rows, row{step, output, result, detail, seconds(time.Since(started))})
	if result != "equal" && result != "skipped" {
		c.failed = true
	}
}

// compareJSON compares the regenerated file (under the scratch tree) with the
// committed one (under the repository root), both at the same relative path
// unless committed names another file.
func (c *checker) compareJSON(step, rel, committed string, started time.Time) {
	if committed == "" {
		committed = filepath.Join(c.root, filepath.FromSlash(rel))
	}
	regenerated, err := readJSON(filepath.Join(c.scratch, filepath.FromSlash(rel)))
	if err != nil {
		c.record(step, rel, "error", err.Error(), started)
		return
	}
	shipped, err := readJSON(committed)
	if err != nil {
		c.record(step, rel, "error", err.Error(), started)
		return
	}
	if d := diffJSON(shipped, regenerated, "$"); d != "" {
		c.record(step, rel, "differs", d, started)
		return
	}
	c.record(step, rel, "equal", "", started)
}

func (c *checker) compareBytes(step, rel string, regenerated []byte, started time.Time) {
	committed, err := os.ReadFile(filepath.Join(c.root, filepath.FromSlash(rel)))
	if err != nil {
		c.record(step, rel, "error", err.Error(), started)
		return
	}
	if string(committed) != string(regenerated) {
		c.record(step, rel, "differs", fmt.Sprintf("%d committed bytes vs %d regenerated", len(committed), len(regenerated)), started)
		return
	}
	c.record(step, rel, "equal", "", started)
}

func main() {
	rootFlag := flag.String("root", "", "repository root (default: nearest ancestor of the working directory with results/ and tools/go.mod)")
	flag.Parse()
	if flag.NArg() != 0 {
		fmt.Fprintln(os.Stderr, "usage: research-checks [-root DIR]")
		os.Exit(2)
	}
	root, err := findRoot(*rootFlag)
	if err != nil {
		fatalf("%v", err)
	}
	c := &checker{root: root, goTool: findGo()}
	c.scratch = filepath.Join(root, "tmp", "research-checks")
	c.bins = filepath.Join(c.scratch, "bin")
	c.logs = filepath.Join(c.scratch, "logs")
	if err := os.RemoveAll(c.scratch); err != nil {
		fatalf("%v", err)
	}
	// The scratch directory must contain results/ and tools/go.mod for the
	// programs' -root check.
	for _, dir := range []string{c.bins, c.logs, "results", "research/broad_cells", "research/broad_affine", "research/en_families", "research/en_independent", "sources", "tools"} {
		if !filepath.IsAbs(dir) {
			dir = filepath.Join(c.scratch, filepath.FromSlash(dir))
		}
		if err := os.MkdirAll(dir, 0o755); err != nil {
			fatalf("%v", err)
		}
	}
	if err := copyFile(filepath.Join(root, "tools", "go.mod"), filepath.Join(c.scratch, "tools", "go.mod")); err != nil {
		fatalf("%v", err)
	}
	scratchRoot := []string{"-root", c.scratch}

	// sources-manifest: hashes the local PDFs; sources/ is a gitignored local
	// archive, so the comparison is with the local manifest.
	started := time.Now()
	exe := c.build("tools", "sources-manifest")
	localSources := filepath.Join(root, "sources")
	if _, err := os.Stat(filepath.Join(localSources, "manifest.json")); err != nil {
		c.record("sources-manifest", "sources/manifest.json", "skipped", "sources/ is absent (a gitignored local archive of the cited PDFs); nothing to compare", started)
	} else {
		entries, err := os.ReadDir(localSources)
		if err != nil {
			fatalf("%v", err)
		}
		for _, e := range entries {
			if e.Type().IsRegular() && strings.HasSuffix(e.Name(), ".pdf") {
				if err := copyFile(filepath.Join(localSources, e.Name()), filepath.Join(c.scratch, "sources", e.Name())); err != nil {
					fatalf("%v", err)
				}
			}
		}
		if _, err := c.run("sources-manifest", c.scratch, nil, append([]string{exe}, scratchRoot...)...); err != nil {
			c.record("sources-manifest", "sources/manifest.json", "error", err.Error(), started)
		} else {
			c.compareJSON("sources-manifest", "sources/manifest.json", "", started)
		}
	}

	// i6-shape: two certificates of the shape of Gern's D6 interval.
	started = time.Now()
	exe = c.build("tools", "i6-shape")
	if _, err := c.run("i6-shape", c.scratch, nil, append([]string{exe}, scratchRoot...)...); err != nil {
		c.record("i6-shape", "results/i6-shape-certificate.json", "error", err.Error(), started)
	} else {
		c.compareJSON("i6-shape", "results/i6-shape-certificate.json", "", started)
		c.compareJSON("i6-shape", "results/i6-independent-subword.json", "", started)
	}

	// exceptional-cosets: text on stdout, compared byte for byte.
	started = time.Now()
	exe = c.build("tools", "exceptional-cosets")
	out, err := c.run("exceptional-cosets", c.scratch, nil, exe)
	if err != nil {
		c.record("exceptional-cosets", "research/exceptional-cosets.txt", "error", err.Error(), started)
	} else {
		c.compareBytes("exceptional-cosets", "research/exceptional-cosets.txt", out, started)
	}

	// star-family: the star-family KL polynomials and the affine D4 dependency
	// certificate.
	started = time.Now()
	exe = c.build("tools", "star-family")
	if _, err := c.run("star-family", c.scratch, nil, append([]string{exe}, scratchRoot...)...); err != nil {
		c.record("star-family", "research/broad_cells/star_family_checks.json", "error", err.Error(), started)
	} else {
		c.compareJSON("star-family", "research/broad_cells/star_family_checks.json", "", started)
		c.compareJSON("star-family", "research/broad_cells/affine_d4_kl_certificate.json", "", started)
	}

	// affine-d-checks: the affine D4 unfolding obstruction and Gern coatoms.
	started = time.Now()
	exe = c.build("tools", "affine-d-checks")
	if _, err := c.run("affine-d-checks", c.scratch, nil, append([]string{exe}, scratchRoot...)...); err != nil {
		c.record("affine-d-checks", "research/broad_affine/unfolding_obstruction.json", "error", err.Error(), started)
	} else {
		c.compareJSON("affine-d-checks", "research/broad_affine/unfolding_obstruction.json", "", started)
		c.compareJSON("affine-d-checks", "research/broad_affine/gern_coatoms.json", "", started)
	}

	// catalogue-eligible for rank 9, after generating the E9 FC catalogue it
	// reads (a gitignored output of the supplement's fc-catalogue).
	started = time.Now()
	fc := c.build("supplement", "fc-catalogue")
	exe = c.build("tools", "catalogue-eligible")
	catalogue := filepath.Join(c.scratch, "research", "en_independent", "e9-fc.json")
	eligible := "research/en_families/e9_full_support_eligible.json"
	out, err = c.run("E9-FC", c.scratch, nil, fc, "-rank", "9")
	if err == nil {
		err = os.WriteFile(catalogue, out, 0o644)
	}
	if err != nil {
		c.record("catalogue-eligible", eligible, "error", err.Error(), started)
	} else if _, err := c.run("catalogue-eligible", c.scratch, nil, exe, "-rank", "9", "-root", c.scratch,
		"-conjugates", filepath.Join(root, "research", "en_families", "conjugates_E9_50000_934.json"),
		"-catalogue", catalogue, "-out", filepath.Join(c.scratch, filepath.FromSlash(eligible))); err != nil {
		c.record("catalogue-eligible", eligible, "error", err.Error(), started)
	} else {
		c.compareJSON("catalogue-eligible", eligible, "", started)
	}

	fmt.Print(table(c.rows))
	counts := map[string]int{}
	for _, r := range c.rows {
		counts[r.result]++
	}
	keys := make([]string, 0, len(counts))
	for k := range counts {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	parts := make([]string, 0, len(keys))
	for _, k := range keys {
		parts = append(parts, fmt.Sprintf("%d %s", counts[k], k))
	}
	fmt.Printf("research-checks: %s; scratch tree tmp/research-checks/\n", strings.Join(parts, ", "))
	if c.failed {
		os.Exit(1)
	}
}

// table renders the rows with aligned columns; a detail (the first JSON
// difference or an error) follows its row indented.
func table(rows []row) string {
	w := [3]int{len("step"), len("output"), len("result")}
	for _, r := range rows {
		for i, s := range []string{r.step, r.output, r.result} {
			if len(s) > w[i] {
				w[i] = len(s)
			}
		}
	}
	var b strings.Builder
	fmt.Fprintf(&b, "%-*s  %-*s  %-*s  %s\n", w[0], "step", w[1], "output", w[2], "result", "seconds")
	for _, r := range rows {
		fmt.Fprintf(&b, "%-*s  %-*s  %-*s  %.2f\n", w[0], r.step, w[1], r.output, w[2], r.result, r.seconds)
		if r.detail != "" {
			fmt.Fprintf(&b, "    %s\n", r.detail)
		}
	}
	return b.String()
}
