// Command proofs rebuilds and verifies the local proof supplement.
//
// Default mode recomputes every classification, catalogue and certificate
// that the manuscript's proofs depend on, verifies the shipped certificates
// against the regenerated ones, and compares the stable mathematical
// outcomes with expected-summary.json. -full additionally reruns the
// historical full-group E6/E7 searches, the second (geometric) model of the
// D6/D8 Kazhdan-Lusztig computation and the larger uniform-family checks,
// comparing every regenerated snapshot with the shipped one.
// -finite runs only the finite theorem's terminal classifications, printed
// table checks and D6 coefficient checks, comparing their stable outcomes
// with the corresponding part of expected-summary.json.
//
// It ports supplement/run_proofs.py of release v0.2.0: the same step names,
// the same snapshot comparisons (ignoring the same volatile timing keys), the
// same expected-summary construction from the same certificate fields, the
// same manuscript-hash check against INPUTS.json and the same status.json
// fields. Instead of compiling C++ and running Python it builds each Go
// command of this module with `go build` into the run's bin/ directory and
// executes it inside the isolated working copy of payload/.
//
// Run from the archive root (the directory holding go.mod, MANIFEST.json and
// payload/): `go run ./cmd/proofs [-finite | -full]`. Standard library only; it imports
// none of the engines or internal packages of this module, so it cannot
// share arithmetic with any program it verifies.
package main

import (
	"crypto/sha256"
	"encoding/hex"
	"errors"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sort"
	"strings"
	"time"
)

// volatileKeys are the wall-clock timing keys written by the terminal
// engines; they are ignored wherever snapshots are compared.
var volatileKeys = map[string]bool{"seconds": true, "elapsed_seconds": true}

func stripVolatile(v any) any {
	switch x := v.(type) {
	case map[string]any:
		out := make(map[string]any, len(x))
		for k, val := range x {
			if !volatileKeys[k] {
				out[k] = stripVolatile(val)
			}
		}
		return out
	case []any:
		out := make([]any, len(x))
		for i, val := range x {
			out[i] = stripVolatile(val)
		}
		return out
	}
	return v
}

// canonical reads a JSON file without its volatile keys.
func canonical(path string) any {
	v, err := readJSON(path)
	if err != nil {
		panic(failure{err.Error()})
	}
	return stripVolatile(v)
}

type step struct {
	Step       string   `json:"step"`
	Command    []string `json:"command"`
	Returncode int      `json:"returncode"`
	Seconds    float64  `json:"seconds"`
}

type status struct {
	Status                                string   `json:"status"`
	Mode                                  string   `json:"mode"`
	Seconds                               float64  `json:"seconds"`
	GoVersion                             string   `json:"go_version"`
	Steps                                 int      `json:"steps"`
	Work                                  string   `json:"work"`
	ExpectedOutputsMatch                  bool     `json:"expected_outputs_match"`
	PrintedTableMatchesGeneratedTerminals bool     `json:"printed_table_matches_generated_terminals"`
	ManuscriptSHA256MatchesBuildRecord    bool     `json:"manuscript_sha256_matches_build_record"`
	SnapshotsCompared                     []string `json:"snapshots_compared"`
}

type runner struct {
	root, run, work, logs, bins string
	goTool                      string
	full                        bool
	finite                      bool
	started                     time.Time
	steps                       []step
	snapshots                   []string
}

func seconds(d time.Duration) float64 { return math.Round(d.Seconds()*1000) / 1000 }

func digest(path string) string {
	f, err := os.Open(path)
	if err != nil {
		panic(failure{err.Error()})
	}
	defer f.Close()
	h := sha256.New()
	if _, err := io.Copy(h, f); err != nil {
		panic(failure{err.Error()})
	}
	return hex.EncodeToString(h.Sum(nil))
}

func must(err error) {
	if err != nil {
		panic(failure{err.Error()})
	}
}

// rel renders a path relative to the archive root with forward slashes; the
// run's records never contain machine-specific absolute paths.
func (r *runner) rel(path string) string {
	p, err := filepath.Rel(r.root, path)
	must(err)
	return filepath.ToSlash(p)
}

func copyFile(src, dst string) error {
	in, err := os.Open(src)
	if err != nil {
		return err
	}
	defer in.Close()
	info, err := in.Stat()
	if err != nil {
		return err
	}
	if err := os.MkdirAll(filepath.Dir(dst), 0o755); err != nil {
		return err
	}
	out, err := os.OpenFile(dst, os.O_CREATE|os.O_TRUNC|os.O_WRONLY, info.Mode().Perm())
	if err != nil {
		return err
	}
	if _, err := io.Copy(out, in); err != nil {
		out.Close()
		return err
	}
	return out.Close()
}

func copyTree(src, dst string) error {
	return filepath.WalkDir(src, func(path string, d os.DirEntry, err error) error {
		if err != nil {
			return err
		}
		relPath, err := filepath.Rel(src, path)
		if err != nil {
			return err
		}
		target := filepath.Join(dst, relPath)
		if d.IsDir() {
			return os.MkdirAll(target, 0o755)
		}
		return copyFile(path, target)
	})
}

// execute runs a command (in the work directory, or in dir when given),
// capturing stdout and stderr to logs/<name>.*.log and recording the step.
// When output is set, stdout is copied to that path inside the work tree.
func (r *runner) execute(name string, command []string, output string, dir string) {
	if dir == "" {
		dir = r.work
	}
	before := time.Now()
	outFile, err := os.Create(filepath.Join(r.logs, name+".stdout.log"))
	must(err)
	errFile, err := os.Create(filepath.Join(r.logs, name+".stderr.log"))
	must(err)
	cmd := exec.Command(command[0], command[1:]...)
	cmd.Dir = dir
	if command[0] == r.goTool {
		// The build steps must not see an enclosing Go workspace or the
		// caller's GOFLAGS: a go.work file in a directory above the extracted
		// archive would otherwise exclude this module from every go build.
		cmd.Env = append(os.Environ(), "GOWORK=off", "GOFLAGS=")
	}
	cmd.Stdout = outFile
	cmd.Stderr = errFile
	runErr := cmd.Run()
	must(outFile.Close())
	must(errFile.Close())
	code := 0
	if runErr != nil {
		var exitErr *exec.ExitError
		if errors.As(runErr, &exitErr) {
			code = exitErr.ExitCode()
		} else {
			panic(failure{name + ": " + runErr.Error()})
		}
	}
	recorded := make([]string, len(command))
	for i, c := range command {
		if c == r.goTool {
			c = "go"
		}
		recorded[i] = strings.ReplaceAll(filepath.ToSlash(c), filepath.ToSlash(r.run), "<run>")
	}
	r.steps = append(r.steps, step{Step: name, Command: recorded, Returncode: code, Seconds: seconds(time.Since(before))})
	data, err := encodeJSON(r.steps)
	must(err)
	must(os.WriteFile(filepath.Join(r.run, "steps.json"), data, 0o644))
	if code != 0 {
		panic(failure{fmt.Sprintf("%s failed; inspect %s", name, r.rel(r.logs))})
	}
	if output != "" {
		must(copyFile(filepath.Join(r.logs, name+".stdout.log"), filepath.Join(r.work, filepath.FromSlash(output))))
	}
	fmt.Printf("%s: passed\n", name)
}

// build compiles one command of this module into bin/ and returns its path.
func (r *runner) build(name string) string {
	exe := filepath.Join(r.bins, name)
	if runtime.GOOS == "windows" {
		exe += ".exe"
	}
	r.execute("build-"+name, []string{r.goTool, "build", "-o", exe, "./cmd/" + name}, "", r.root)
	return exe
}

func (r *runner) read(path string) any {
	v, err := readJSON(filepath.Join(r.work, filepath.FromSlash(path)))
	must(err)
	return v
}

// compareSnapshot requires the regenerated file to equal the shipped one up
// to volatile timing keys; select restricts the comparison to a subtree.
func (r *runner) compareSnapshot(path string, selectFn func(any) any) {
	if selectFn == nil {
		selectFn = func(v any) any { return v }
	}
	regenerated := selectFn(canonical(filepath.Join(r.work, filepath.FromSlash(path))))
	shipped := selectFn(canonical(filepath.Join(r.root, "payload", filepath.FromSlash(path))))
	if equal, diff := sameJSON(regenerated, shipped); !equal {
		panic(failure{fmt.Sprintf("Regenerated %s differs from the shipped snapshot; inspect %s (first difference %s)", path, r.rel(filepath.Join(r.work, path)), diff)})
	}
	r.snapshots = append(r.snapshots, path)
	fmt.Printf("snapshot %s: regenerated output equals the shipped file\n", path)
}

func findGo() string {
	if p, err := exec.LookPath("go"); err == nil {
		return p
	}
	if root := os.Getenv("GOROOT"); root != "" {
		p := filepath.Join(root, "bin", "go")
		if _, err := os.Stat(p); err == nil {
			return p
		}
	}
	fmt.Fprintln(os.Stderr, "The go command is required to build the programs of this supplement; it was not found on PATH.")
	os.Exit(1)
	return ""
}

func main() {
	full := flag.Bool("full", false, "also rerun the historical E6/E7 full-group searches, the geometric D6/D8 model and the larger uniform-family checks")
	finite := flag.Bool("finite", false, "verify only the finite terminal classifications, printed table and D6 coefficient")
	flag.Parse()
	if flag.NArg() != 0 {
		fmt.Fprintln(os.Stderr, "proofs: positional arguments are not used")
		os.Exit(2)
	}
	mode, err := proofMode(*full, *finite)
	if err != nil {
		fmt.Fprintln(os.Stderr, "proofs: "+err.Error())
		os.Exit(2)
	}
	root, err := os.Getwd()
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	for _, required := range []string{"go.mod", "MANIFEST.json", "INPUTS.json", "expected-summary.json", "payload"} {
		if _, err := os.Stat(filepath.Join(root, required)); err != nil {
			fmt.Fprintf(os.Stderr, "Run this command from the archive root (the directory containing go.mod, MANIFEST.json and payload/); %s is missing here.\n", required)
			os.Exit(1)
		}
	}
	r := &runner{root: root, full: *full, finite: *finite, goTool: findGo()}

	// Integrity of every distributed input, then the manuscript hash
	// recorded by the builder.
	manifestDoc, err := readJSON(filepath.Join(root, "MANIFEST.json"))
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	var manuscriptSHA256AtBuild string
	func() {
		defer func() {
			if rec := recover(); rec != nil {
				msg := fmt.Sprint(rec)
				if f, ok := rec.(failure); ok {
					msg = f.msg
				}
				fmt.Fprintln(os.Stderr, msg)
				os.Exit(1)
			}
		}()
		hashes := mapping(field(manifestDoc, "sha256"))
		for _, name := range sortedKeys(hashes) {
			file := filepath.Join(root, filepath.FromSlash(name))
			want, ok := hashes[name].(string)
			if !ok {
				panic(failure{"MANIFEST.json: non-string digest for " + name})
			}
			info, err := os.Stat(file)
			if err != nil || !info.Mode().IsRegular() || digest(file) != want {
				panic(failure{"Input checksum mismatch: " + name})
			}
		}
		inputs, err := readJSON(filepath.Join(root, "INPUTS.json"))
		must(err)
		for _, f := range items(field(inputs, "files")) {
			if field(f, "source_path") == "results/exceptional-leading.tex" {
				manuscriptSHA256AtBuild = field(f, "source_sha256").(string)
				break
			}
		}
		if manuscriptSHA256AtBuild == "" {
			panic(failure{"INPUTS.json does not record the manuscript hash."})
		}
	}()

	runs := filepath.Join(root, "runs")
	must(os.MkdirAll(runs, 0o755))
	r.run, err = os.MkdirTemp(runs, mode+"-")
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	r.work = filepath.Join(r.run, "work")
	r.logs = filepath.Join(r.run, "logs")
	r.bins = filepath.Join(r.run, "bin")
	r.started = time.Now()

	defer func() {
		if rec := recover(); rec != nil {
			msg := fmt.Sprint(rec)
			if f, ok := rec.(failure); ok {
				msg = f.msg
			}
			fmt.Fprintln(os.Stderr, msg)
			fmt.Fprintln(os.Stderr, "Verification failed. Retained working tree and logs: "+r.rel(r.run))
			os.Exit(1)
		}
	}()
	must(copyTree(filepath.Join(root, "payload"), r.work))
	must(os.Mkdir(r.logs, 0o755))
	must(os.Mkdir(r.bins, 0o755))

	r.verify(manuscriptSHA256AtBuild)

	st := status{Status: "passed", Mode: mode, Seconds: seconds(time.Since(r.started)), GoVersion: runtime.Version(),
		Steps: len(r.steps), Work: r.rel(r.work), ExpectedOutputsMatch: true, PrintedTableMatchesGeneratedTerminals: true,
		ManuscriptSHA256MatchesBuildRecord: true, SnapshotsCompared: r.snapshots}
	data, err := encodeJSON(st)
	must(err)
	must(os.WriteFile(filepath.Join(r.run, "status.json"), data, 0o644))
	fmt.Println("All proof checks and expected outputs passed. Logs and certificates: " + r.rel(r.run))
}

// verify runs every step of the supplement in the order of the original
// runner and raises a failure on the first problem.
func (r *runner) verify(manuscriptSHA256AtBuild string) {
	if r.finite {
		r.verifyFinite(manuscriptSHA256AtBuild)
		return
	}
	flat := r.build("terminals-flat")
	recursive := r.build("terminals-recursive")
	d7 := r.build("e8-d7")
	fc := r.build("fc-catalogue")
	e6mu := r.build("e6-mu-table")
	e9kl := r.build("e9-quotient-kl")
	d8kl := r.build("d8-gern-kl")
	fcmax := r.build("fc-maxima")
	// Programs that the original ran as Python scripts are compiled too.
	verifyExceptional := r.build("verify-exceptional")
	verifyOutputs := r.build("verify-outputs")
	verifyE8 := r.build("verify-e8")
	finiteDescents := r.build("finite-descents")
	d6Certificate := r.build("d6-certificate")
	uniformVerify := r.build("uniform-verify")
	affineProof := r.build("affine-proof")
	affineFCCovers := r.build("affine-fc-covers")
	e10AllK := r.build("e10-all-k")
	terminalStructure := r.build("terminal-structure")
	uniformCovers := r.build("uniform-covers")
	affineD4 := r.build("affine-d4")
	affineReflectionFamily := r.build("affine-reflection-family")
	cartanCandidates := r.build("cartan-candidates")
	uniformConstruction := r.build("uniform-construction")
	terminalDataCheck := r.build("terminal-data-check")
	d6Ambient := r.build("d6-ambient")

	for _, n := range []int{6, 7, 8, 9} {
		r.execute(fmt.Sprintf("E%d-FC", n), []string{fc, "-rank", fmt.Sprint(n)}, fmt.Sprintf("research/en_independent/e%d-fc.json", n), "")
	}
	for _, n := range []int{6, 7, 8} {
		flatFile := fmt.Sprintf("research/en_e8/e%d-validation.json", n)
		if n == 8 {
			flatFile = "research/en_e8/e8-terminals.json"
		}
		r.execute(fmt.Sprintf("E%d-flat", n), []string{flat, "-rank", fmt.Sprint(n)}, flatFile, "")
		r.execute(fmt.Sprintf("E%d-recursive", n), []string{recursive, "-rank", fmt.Sprint(n)}, fmt.Sprintf("research/en_e8/e%d-recursive.json", n), "")
	}
	r.execute("E8-D7", []string{d7}, "research/en_independent/e8-d7-terminals.json", "")
	// The E8 classification is repeated along four parabolic chains (the
	// paper's chain through E7, and chains through D7, A7 and a second E7
	// ordering); the program asserts that all four give the same 64 terminal
	// elements and records the per-chain coset and terminal counts.
	r.execute("E8-recursive-chains", []string{recursive, "-rank", "8", "-chains-certificate", "research/en_e8/e8-recursive-chains.json"}, "", "")
	r.compareSnapshot("research/en_e8/e8-recursive-chains.json", nil)
	if r.full {
		rootids := r.build("enumerate-bad")
		matrix := r.build("matrix-search")
		for _, n := range []int{6, 7} {
			path := fmt.Sprintf("research/broad_exceptional/e%d_bad.json", n)
			r.execute(fmt.Sprintf("E%d-historical-rootids", n), []string{rootids, "-rank", fmt.Sprint(n)}, path, "")
			r.compareSnapshot(path, nil)
		}
		// The E6 search writes its certificate itself; the E7 search prints it.
		r.execute("E6-historical-matrix", []string{matrix, "-rank", "6"}, "", "")
		r.compareSnapshot("results/e6-independent-certificate.json", nil)
		r.execute("E7-historical-matrix", []string{matrix, "-rank", "7"}, "results/e7-independent-certificate.json", "")
		r.compareSnapshot("results/e7-independent-certificate.json", nil)
	}
	r.execute("finite-source-snapshot-audit", []string{verifyExceptional}, "", "")
	r.execute("finite-method-comparison", []string{verifyOutputs}, "", "")
	r.execute("E8-independent-verification", []string{verifyE8}, "", "")
	r.execute("finite-support-matchings", []string{finiteDescents}, "", "")
	r.compareSnapshot("research/ai-review-notes/finite-descents.json", nil)
	// Third check of the Table 1 data from the E6/E7 integer-matrix
	// certificates: Poincare polynomials, commuting terminal counts, the E7
	// word transcriptions and the D6 identification.
	r.execute("terminal-data-check", []string{terminalDataCheck}, "research/exceptional_referee/checks.json", "")
	r.compareSnapshot("research/exceptional_referee/checks.json", nil)
	r.execute("D6-recurrence-certificate", []string{d6Certificate}, "", "")
	// The D6 polynomial of the length-15 pair recomputed from scratch by the
	// KL recursion inside ambient E7, ambient E8 and the D6 parabolic.
	r.execute("D6-ambient-E7-E8", []string{d6Ambient, "-out", "results/d6-ambient-certificate.json"}, "", "")
	r.compareSnapshot("results/d6-ambient-certificate.json", nil)
	r.execute("uniform-symbolic-and-matrix-checks", []string{uniformVerify}, "", "")
	// The uniform verifier records the hash of the manuscript it audited; it
	// must be the manuscript snapshot recorded by the builder in INPUTS.json.
	audited := field(r.read("research/en_uniform/referee-certificate.json"), "audited_source_sha256")
	if audited != manuscriptSHA256AtBuild {
		panic(failure{"Uniform verifier audited a manuscript other than the one recorded at build time."})
	}
	// Finite counterchecks of the rank-uniform construction for r = 3..30.
	r.execute("uniform-construction", []string{uniformConstruction, "-max-r", "30"}, "", "")
	r.compareSnapshot("research/en_uniform/construction.json", nil)
	r.execute("affine-proved-family", []string{affineProof}, "", "")
	r.compareSnapshot("research/en_affine_referee/proved-affine-certificate.json", nil)
	r.execute("affine-FC-closure-and-covers", []string{affineFCCovers}, "", "")
	// The length-33 reflection family r_{beta+k delta} of affine E8: terminal
	// with full support for k = 0..10, length 33 + 58k.
	r.execute("affine-reflection-family", []string{affineReflectionFamily}, "", "")
	r.compareSnapshot("research/en_families/affine_reflection_family.json", nil)
	// The E10 seed (the real terminal roots with pairings bounded by 2 on
	// maximum independent sets) is regenerated and compared before e10-all-k
	// reads its beta0 row from the regenerated file.
	r.execute("cartan-candidates", []string{cartanCandidates, "-rank", "10", "-max-entry", "2", "-max-only"}, "", "")
	r.compareSnapshot("research/en_families/cartan_E10_m2_max1.json", nil)
	r.execute("E10-all-k", []string{e10AllK}, "", "")
	// Certificates added after the review of 7 October 2026 (Section 5 of the
	// manuscript and the remarks of Sections 3 and 4). Each program asserts its
	// expected values itself; the runner additionally requires equality with
	// the shipped certificate.
	r.execute("E6-complete-mu-table", []string{e6mu, "-type", "E", "-rank", "6", "-out", "results/e6-mu-table.json"}, "", "")
	r.compareSnapshot("results/e6-mu-table.json", nil)
	r.execute("E9-odd-gap-polynomials", []string{e9kl, "-mode", "certify", "-out", "results/e9-odd-gap-certificate.json"}, "", "")
	r.compareSnapshot("results/e9-odd-gap-certificate.json", nil)
	models := "sp"
	if r.full {
		models = "both"
	}
	r.execute("D6-D8-Gern-polynomials", []string{d8kl, "-out", "results/d8-gern-certificate.json", "-threads", "0", "-models", models}, "", "")
	r.compareSnapshot("results/d8-gern-certificate.json", func(d any) any { return []any{field(d, "ranks"), field(d, "status")} })
	r.execute("FC-maxima-E6-E13", []string{fcmax, "-out", "results/fc-maxima-certificate.json", "-ranks", "6,7,8,9,10,11,12,13"}, "", "")
	r.compareSnapshot("results/fc-maxima-certificate.json", nil)
	r.execute("terminal-structure", []string{terminalStructure}, "", "")
	r.compareSnapshot("results/terminal-structure-certificate.json", nil)
	uniformCommand := []string{uniformCovers}
	if r.full {
		uniformCommand = append(uniformCommand, "-full")
	}
	r.execute("uniform-family-covers", uniformCommand, "", "")
	if r.full {
		// The full run checks more pairs; the shipped (default) rows must reappear unchanged.
		regenerated := uniformRows(canonical(filepath.Join(r.work, "results/uniform-family-certificate.json")))
		shipped := uniformRows(canonical(filepath.Join(r.root, "payload/results/uniform-family-certificate.json")))
		ok := regenerated.status == "passed" && shipped.status == "passed"
		for _, table := range []int{0, 1} {
			for key, value := range shipped.rows[table] {
				got, present := regenerated.rows[table][key]
				if !present {
					ok = false
					continue
				}
				if equal, _ := sameJSON(got, value); !equal {
					ok = false
				}
			}
		}
		if !ok {
			panic(failure{"Full uniform-family run disagrees with the shipped certificate."})
		}
		r.snapshots = append(r.snapshots, "results/uniform-family-certificate.json (shipped rows)")
	} else {
		r.compareSnapshot("results/uniform-family-certificate.json", nil)
	}
	r.execute("affine-D4-mu-2", []string{affineD4}, "", "")
	r.compareSnapshot("results/affine-d4-independent.json", nil)

	r.compareSummary(r.summary(), nil)
	r.verifyPrintedTable(verifyOutputs)
}

// compareSummary saves only the outcomes recomputed by this mode. A nil
// selection compares the full summary; finite mode selects its proof inputs.
func (r *runner) compareSummary(summary any, selection []string) {
	expected, err := readJSON(filepath.Join(r.root, "expected-summary.json"))
	must(err)
	if selection != nil {
		expected = pick(expected, selection...)
	}
	data, err := encodeJSON(summary)
	must(err)
	must(os.WriteFile(filepath.Join(r.run, "summary.json"), data, 0o644))
	if equal, diff := sameJSON(summary, expected); !equal {
		panic(failure{"Stable proof outputs differ from expected-summary.json; inspect " + r.rel(filepath.Join(r.run, "summary.json")) + " (first difference " + diff + ")"})
	}
}

func (r *runner) verifyPrintedTable(verifyOutputs string) {
	// Also match the table words checked by the independent matching verifier
	// to the actual generated terminal matrices, not merely their lengths.
	matching := r.read("research/ai-review-notes/finite-descents.json")
	for _, n := range []int{6, 7, 8} {
		name := fmt.Sprintf("E%d-terminal-matrices", n)
		r.execute(name, []string{verifyOutputs, "-matrices", fmt.Sprint(n)}, "", "")
		generated, err := readJSON(filepath.Join(r.logs, name+".stdout.log"))
		must(err)
		if !tableMatchesGenerated(matching, generated, n) {
			panic(failure{fmt.Sprintf("The printed E%d table words do not generate the regenerated terminal matrices.", n)})
		}
	}
}

// uniformTables indexes the pairs and Bruhat cover rows of the uniform-family
// certificate by (r, k), with the certificate's status.
type uniformTables struct {
	rows   [2]map[string]any
	status string
}

func uniformRows(d any) uniformTables {
	var t uniformTables
	for i, key := range []string{"pairs", "bruhat_covers"} {
		t.rows[i] = map[string]any{}
		for _, row := range items(field(d, key)) {
			t.rows[i][fmt.Sprintf("%d,%d", integer(field(row, "r")), integer(field(row, "k")))] = row
		}
	}
	t.status, _ = field(d, "status").(string)
	return t
}

// tableMatchesGenerated recomputes the integer matrix (tuple of columns in
// the simple-root basis) of every printed E_n table word and requires the set
// to equal the generated terminal matrices printed by verify-outputs.
func tableMatchesGenerated(matching, generated any, n int) bool {
	edges := make([][2]int, 0, n-1)
	for j := 0; j < n-2; j++ {
		edges = append(edges, [2]int{j, j + 1})
	}
	edges = append(edges, [2]int{2, n - 1})
	matrix := func(word string) string {
		cols := make([][]int64, n)
		for j := range cols {
			cols[j] = make([]int64, n)
			cols[j][j] = 1
		}
		for _, ch := range word {
			s := int(ch - '0')
			if s < 0 || s >= n {
				panic(failure{fmt.Sprintf("table word %q has a letter outside 0..%d", word, n-1)})
			}
			old := cols
			out := make([][]int64, n)
			copy(out, old)
			neg := make([]int64, n)
			for i, v := range old[s] {
				neg[i] = -v
			}
			out[s] = neg
			for _, e := range edges {
				if s == e[0] || s == e[1] {
					t := e[0]
					if s == e[0] {
						t = e[1]
					}
					sum := make([]int64, n)
					for i := range sum {
						sum[i] = old[t][i] + old[s][i]
					}
					out[t] = sum
				}
			}
			cols = out
		}
		return fmt.Sprint(cols)
	}
	fromTable := map[string]bool{}
	for _, row := range items(matching) {
		if field(row, "type") == fmt.Sprintf("E%d", n) {
			fromTable[matrix(field(row, "word").(string))] = true
		}
	}
	fromEngine := map[string]bool{}
	for _, m := range items(generated) {
		cols := make([][]int64, 0, n)
		for _, col := range items(m) {
			entries := make([]int64, 0, n)
			for _, v := range items(col) {
				entries = append(entries, integer(v))
			}
			cols = append(cols, entries)
		}
		fromEngine[fmt.Sprint(cols)] = true
	}
	if len(fromTable) != len(fromEngine) {
		return false
	}
	for k := range fromTable {
		if !fromEngine[k] {
			return false
		}
	}
	return true
}

// sortedInt64 returns a sorted copy.
func sortedInt64(values []int64) []int64 {
	out := append([]int64(nil), values...)
	sort.Slice(out, func(i, j int) bool { return out[i] < out[j] })
	return out
}
