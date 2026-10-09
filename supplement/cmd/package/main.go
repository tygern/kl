// Command package collects the explicit, minimal proof supplement archive
// for release v0.5.0; no network, no installations.
//
// It ports supplement/build_package.py: the same payload inventory of data
// files (shipped certificates, snapshots, the manuscript snapshot and its
// rendered PDF) with their purposes, the same seed reduction, LICENSE,
// the archive README (supplement/README.md), INPUTS.json with version, scope
// and provenance, expected-summary.json computed from the shipped
// certificates with the same assertions, MANIFEST.json over every archive
// file but itself, and a deterministic ZIP (sorted paths, timestamps
// 2000-01-01, mode 0644). Instead of Python/C++ programs the archive carries
// this Go module: go.mod and every source file under cmd/ and internal/
// except this packager, which only runs in the repository.
//
// Run from the repository root as `go run -C supplement ./cmd/package` (the module
// lives in supplement/, so the go command changes into it; the packager then
// treats the parent of its module directory as the repository root). It writes
// supplement/dist/exceptional-leading-proof/ (staging),
// supplement/dist/exceptional-leading-proof.zip and
// supplement/dist/build-summary.json. Standard library only; it imports no
// other package of this module.
package main

import (
	"archive/zip"
	"compress/flate"
	"crypto/sha256"
	"encoding/hex"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"time"
)

const version = "v0.5.0"

// renamed maps the current location of the internal review notes to their
// pre-rename location; the archive always uses the new path.
var renamed = map[string]string{"research/ai-review-notes/": "research/journal-review/"}

// files is the payload inventory: archive path (under payload/) and purpose.
var files = map[string]string{
	"results/proof-study-companion.md":                          "Study companion with worked examples, exercises, and solution sketches; not a computational proof input.",
	"output/pdf/exceptional-leading.pdf":                        "Rendered manuscript snapshot supplied alongside the TeX; not executed by the proof runner.",
	"results/exceptional-leading.tex":                           "Manuscript snapshot; the uniform verifier records its SHA256, which the runner compares with the value recorded here at build time.",
	"research/broad_exceptional/e6_bad.json":                    "Full E6 root-index snapshot written by enumerate-bad -rank 6; regenerated and compared with -full.",
	"research/broad_exceptional/e7_bad.json":                    "Full E7 root-index snapshot written by enumerate-bad -rank 7; regenerated and compared with -full.",
	"results/e6-independent-certificate.json":                   "Full E6 integer-matrix snapshot written by matrix-search -rank 6; regenerated and compared with -full.",
	"results/e7-independent-certificate.json":                   "Full E7 integer-matrix snapshot written by matrix-search -rank 7; regenerated and compared with -full.",
	"research/en_uniform/referee.txt":                           "AI-generated internal review note (not peer review) explaining the all-rank argument referenced by the uniform verifier.",
	"research/en_affine_referee/proved-affine-certificate.json": "Shipped certificate of the proved affine family written by affine-proof (matrix of c_0 by rows, column slopes, inversion data, eligible lower endpoint); regenerated and compared.",
	"research/ai-review-notes/finite-descents.json":             "Expected finite table matching certificate; regenerated and compared.",
	"results/d6-recurrence-certificate.json":                    "All 24,245 saved D6 dependency records, checked without calling a KL evaluator.",
	"results/affine-d4-independent.json":                        "Shipped affine D4 certificate; regenerated and compared.",
	"results/e6-mu-table.json":                                  "Shipped E6 mu-table certificate; regenerated and compared.",
	"results/e9-odd-gap-certificate.json":                       "Shipped affine E8 odd-gap certificate; regenerated and compared.",
	"results/d8-gern-certificate.json":                          "Shipped D6/D8 Gern-element certificate (both models); regenerated and compared.",
	"results/fc-maxima-certificate.json":                        "Shipped FC maxima certificate for E6-E13; regenerated and compared.",
	"results/terminal-structure-certificate.json":               "Shipped structure certificate for the exceptional terminals; regenerated and compared.",
	"results/uniform-family-certificate.json":                   "Shipped uniform-family certificate (default pairs); regenerated and compared.",
	"research/en_families/affine_reflection_family.json":        "Shipped certificate of the affine E8 reflection family r_{beta+k delta}, beta = (1,2,3,3,2,2,1,1,2), written by affine-reflection-family; regenerated and compared.",
	"research/en_families/cartan_E10_m2_max1.json":              "E10 real terminal roots with pairings bounded by 2 on maximum independent sets, written by cartan-candidates -rank 10 -max-entry 2 -max-only; regenerated and compared, and e10-all-k reads its beta0 row.",
	"research/en_uniform/construction.json":                     "Finite counterchecks of the rank-uniform construction for r = 3..30 written by uniform-construction -max-r 30; regenerated and compared.",
	"research/exceptional_referee/checks.json":                  "Third check of the Table 1 data from the E6/E7 integer-matrix certificates, printed by terminal-data-check; regenerated and compared.",
	"research/en_e8/e8-recursive-chains.json":                   "E8 terminal classification repeated along four parabolic chains, written by terminals-recursive -chains-certificate; regenerated and compared.",
	"results/d6-ambient-certificate.json":                       "Kazhdan-Lusztig polynomial of the length-15 pair of Table 1 computed by the KL recursion in ambient E7, ambient E8 and D6, written by d6-ambient; regenerated and compared.",
}

// purposes describes each package of the Go module, keyed by its directory
// relative to go/.
var purposes = map[string]string{
	"cmd/proofs":                   "Proof runner: builds every program below, rebuilds every classification and certificate in an isolated working copy of payload/, compares the regenerated snapshots with the shipped ones and the stable outcomes with expected-summary.json (go run ./cmd/proofs [-full]).",
	"cmd/fc-catalogue":             "Generates complete FC catalogues of E6-E9 (ports fc_catalogue.cpp); no catalogue snapshots shipped.",
	"cmd/terminals-flat":           "Flat finite terminal engine (ports parabolic_terminals.cpp); shares its arithmetic with the recursive engine through internal/parabolic.",
	"cmd/terminals-recursive":      "Recursive complete finite parabolic pruning engine (ports recursive_terminals.cpp); -chains-certificate repeats the E8 classification along four parabolic chains; shares its arithmetic with the flat engine through internal/parabolic.",
	"internal/parabolic":           "Arithmetic shared by the flat and recursive finite terminal engines, so those two are not independent implementations.",
	"cmd/e8-d7":                    "Separate complete E8/D7 terminal enumeration (ports e8_d7_cosets.cpp); imports no internal package.",
	"cmd/enumerate-bad":            "Full E6/E7 root-index enumeration (ports enumerate_bad.cpp; rerun with -full); imports no internal package.",
	"cmd/matrix-search":            "Full E6 and E7 independent integer-matrix enumerations (ports verify_e6_independent.py and verify_e7_matrices.cpp; rerun with -full); imports no internal package.",
	"cmd/verify-exceptional":       "Compares the four historical E6/E7 full-group snapshots (root-index and integer-matrix enumerations) and identifies the signed D6 pair (ports verify_exceptional.py).",
	"cmd/verify-outputs":           "Exact integer-matrix cross-comparison of all finite terminal outputs and baseline snapshots (ports verify_outputs.py); -matrices N prints the generated terminal matrices that the runner compares with the printed table words.",
	"cmd/verify-e8":                "Independent integer-matrix closure, inversion-length, descent and endpoint verifier (ports verify_e8.py); hashes JSON inputs without their timing fields; imports neither enumeration engine.",
	"cmd/finite-descents":          "Independent table descent and support matching verifier (ports check_finite_descents.py from the AI-generated internal review notes).",
	"cmd/d6-certificate":           "D6 saved recurrence certificate checker (ports verify_certificate.py); uses the primitives of internal/d6 and calls no KL evaluator.",
	"internal/d6":                  "D6 signed permutation, length, Bruhat, lower ideal, R-polynomial and polynomial primitives (ports sparse_kl.py and coxeter.py without any KL evaluator); shared with the generator's design, so the checker is not an independent implementation.",
	"cmd/uniform-verify":           "Independent symbolic uniform verifier and finite matrix diagnostics (ports referee_verify.py); records the manuscript hash.",
	"internal/affine":              "Row-matrix affine arithmetic and all-k inversion/terminal verifier (ports the library part of verify_families.py); shared by affine-proof and affine-fc-covers.",
	"cmd/affine-proof":             "Proof-only entry point for the proved affine conjugate family (ports affine_proof.py); writes proved-affine-certificate.json.",
	"cmd/affine-fc-covers":         "E9 exact FC closure and complete base cover verification (ports verify_fc_catalogue.py); independent of the FC enumerator but shares affine arithmetic with affine-proof.",
	"internal/indefinite":          "Arithmetic library of the E10 verifier (ports the library part of verify_indefinite.py; the E13 main is not ported).",
	"cmd/e10-all-k":                "All-k E10 root/reflection proof verifier (ports verify_indefinite_e10.py).",
	"cmd/e6-mu-table":              "Complete E6 Kazhdan-Lusztig table with every mu value: maximum 10 attained by 8 pairs, histogram of mu >= 2, and mu in {0,1} for every FC lower endpoint (manuscript Section 5; ports e6_mu_table.cpp).",
	"cmd/e9-quotient-kl":           "Parabolic-quotient KL engine for the length-33 affine E8 reflection r_beta and its two length-34 extensions; polynomials of the odd-gap eligible bottoms (Section 5; ports e9_quotient_kl.cpp); validated against a naive recursion.",
	"cmd/d8-gern-kl":               "KL polynomials on the lower ideals of Gern's w_6 and w_8 in two independent models of D_n: the D6 value 1+6q+11q^2+6q^3+q^4+q^5 and the D8 polynomial of the archived transfer note (geometric model with -full; ports d8_gern_kl.cpp).",
	"cmd/fc-maxima":                "Independent FC enumeration of E6-E13 by height vectors: counts and maximum lengths (55, 66, 78, 92 for E10-E13; Section 4 remark; ports fc_maxima.cpp).",
	"cmd/terminal-structure":       "Gern-plus-two matrix identities, exhaustive D5/D6/D7 enumeration, layered palindromes, orthogonal reflections, weak-order chain and w_0(J) factorizations of w_4, w_6, w_7, w_8 (Section 3 remarks; ports terminal_structure.py).",
	"cmd/uniform-covers":           "Finite checks of the cover lemma ingredients for b_{r,k}: descents and non-FC status of bs, sbs = r_{s beta}, lengths 8r^2+3+116k, no FC Bruhat covers (Section 4; ports uniform_family_checks.py; more rows with -full).",
	"cmd/affine-d4":                "Independent affine D4 witness P_{x,xcx} = 1+3q+2q^2, mu = 2, by direct subwords and R-polynomial reciprocity (Section 5; ports verify_affine_d4_r.py).",
	"cmd/affine-reflection-family": "Certifies the affine E8 reflection family r_{beta+k delta} (length 33 + 58k, terminal with full support, both one-generator extensions with a full-support FC bottom) for k = 0..10 (Section 5; ports affine_reflection_family.py).",
	"cmd/cartan-candidates":        "Exhausts the terminal reflections of E_n with bounded root pairings on an independent set and certifies each as a real root with the stated descents; regenerates the E10 seed (ports cartan_candidates.py).",
	"cmd/uniform-construction":     "Finite counterchecks of the rank-uniform E_(4r+1) construction for r = 3..30: seed certificate, translation identities, norm-2 positive roots whose reflections are right terminal (Section 4; ports construction.py).",
	"cmd/terminal-data-check":      "Third check of Table 1 from the E6/E7 integer-matrix certificates: Poincare polynomials, commuting terminal counts, E7 word transcriptions and the D6 identification (ports verify_terminal_data.py).",
	"cmd/d6-ambient":               "KL recursion on the lower ideal of the length-15 terminal of Table 1 in ambient E7, ambient E8 and the D6 parabolic: P = 1+6q+11q^2+6q^3+q^4+q^5, mu = 1, ideal 3,184, interval 1,676 (Appendix A; self-contained).",
}

// excludedPackages are module directories that stay in the repository.
var excludedPackages = map[string]bool{"cmd/package": true}

type entry struct {
	ArchivePath    string `json:"archive_path"`
	SourcePath     string `json:"source_path"`
	SourceSHA256   string `json:"source_sha256"`
	Transformation string `json:"transformation"`
	Purpose        string `json:"purpose"`
}

type packager struct{ root, here, dist, stage string }

func fatal(format string, args ...any) {
	fmt.Fprintf(os.Stderr, format+"\n", args...)
	os.Exit(1)
}

func must(err error) {
	if err != nil {
		fatal("%v", err)
	}
}

func digest(path string) string {
	f, err := os.Open(path)
	must(err)
	defer f.Close()
	h := sha256.New()
	_, err = io.Copy(h, f)
	must(err)
	return hex.EncodeToString(h.Sum(nil))
}

func copyFile(src, dst string) {
	must(os.MkdirAll(filepath.Dir(dst), 0o755))
	data, err := os.ReadFile(src)
	must(err)
	must(os.WriteFile(dst, data, 0o644))
}

func save(path string, value any) {
	data, err := encodeJSON(value)
	must(err)
	must(os.MkdirAll(filepath.Dir(path), 0o755))
	must(os.WriteFile(path, data, 0o644))
}

// source resolves a repository file for an archive path, tolerating the
// pre-rename location of the review notes.
func (p *packager) source(path string) string {
	candidate := filepath.Join(p.root, filepath.FromSlash(path))
	if info, err := os.Stat(candidate); err == nil && info.Mode().IsRegular() {
		return candidate
	}
	for newPrefix, oldPrefix := range renamed {
		if strings.HasPrefix(path, newPrefix) {
			old := filepath.Join(p.root, filepath.FromSlash(oldPrefix+path[len(newPrefix):]))
			if info, err := os.Stat(old); err == nil && info.Mode().IsRegular() {
				fmt.Printf("note: %s taken from the pre-rename location %s\n", path, oldPrefix)
				return old
			}
		}
	}
	fatal("Missing source file: %s", path)
	return ""
}

func (p *packager) read(path string) any {
	v, err := readJSON(p.source(path))
	must(err)
	return v
}

// listFiles returns every regular file under dir, as sorted slash paths
// relative to dir.
func listFiles(dir string) []string {
	var out []string
	must(filepath.WalkDir(dir, func(path string, d os.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if d.IsDir() && path != dir {
			switch d.Name() {
			case "dist", "bin", "validation", "runs":
				return filepath.SkipDir
			}
		}
		if d.Type().IsRegular() {
			rel, err := filepath.Rel(dir, path)
			if err != nil {
				return err
			}
			out = append(out, filepath.ToSlash(rel))
		}
		return nil
	}))
	sort.Strings(out)
	return out
}

func main() {
	if len(os.Args) > 1 {
		fatal("package: arguments are not used")
	}
	root := repositoryRoot()
	p := &packager{root: root, here: filepath.Join(root, "supplement")}
	p.dist = filepath.Join(p.here, "dist")
	p.stage = filepath.Join(p.dist, "exceptional-leading-proof")
	defer func() {
		if rec := recover(); rec != nil {
			if f, ok := rec.(failure); ok {
				fatal("%s", f.msg)
			}
			panic(rec)
		}
	}()

	// Only replace this builder's generated staging directory. Runs live in
	// extracted archives, so no proof-review logs are erased here.
	must(os.RemoveAll(p.stage))
	must(os.MkdirAll(p.stage, 0o755))
	var inventory []entry
	archivePaths := make([]string, 0, len(files))
	for path := range files {
		archivePaths = append(archivePaths, path)
	}
	sort.Strings(archivePaths)
	for _, archivePath := range archivePaths {
		origin := p.source(archivePath)
		target := filepath.Join(p.stage, "payload", filepath.FromSlash(archivePath))
		copyFile(origin, target)
		inventory = append(inventory, entry{"payload/" + archivePath, archivePath, digest(origin), "unchanged", files[archivePath]})
	}

	// One discovery-result file is reduced to the exact seed row the proofs use.
	type selection struct {
		path    string
		selects func(any) any
		purpose string
	}
	selections := []selection{
		{"research/en_families/e9_full_support_eligible.json",
			func(data any) any {
				for _, row := range items(data) {
					if integer(field(row, "length")) == 27 {
						return []any{row}
					}
				}
				panic(failure{"no length-27 row in e9_full_support_eligible.json"})
			},
			"Retain only the length-27 base word row used by the proved affine family."},
	}
	for _, s := range selections {
		target := filepath.Join(p.stage, "payload", filepath.FromSlash(s.path))
		save(target, s.selects(p.read(s.path)))
		inventory = append(inventory, entry{"payload/" + s.path, s.path, digest(filepath.Join(p.root, filepath.FromSlash(s.path))), s.purpose,
			"Required fixed proof seed; discovery code and other results excluded."})
	}

	// The Go module is this supplement directory itself: go.mod and every
	// source file under cmd/ and internal/ except the packager, placed at the
	// archive root. dist/, bin/ and validation/ hold build outputs and records,
	// not module sources, and are skipped by the walk.
	goDir := p.here
	for _, rel := range listFiles(goDir) {
		dir := filepath.ToSlash(filepath.Dir(rel))
		var purpose string
		switch {
		case rel == "go.mod":
			purpose = "Go module definition (module github.com/tygern/kl/supplement, Go 1.22 or later, standard library only)."
		case strings.HasSuffix(rel, ".go") && (strings.HasPrefix(rel, "cmd/") || strings.HasPrefix(rel, "internal/")):
			if excludedPackages[dir] {
				continue
			}
			purpose = purposes[dir]
			if purpose == "" {
				fatal("no purpose recorded for Go package %s", dir)
			}
			if strings.HasSuffix(rel, "_test.go") {
				purpose = "Unit test of " + dir + " (go test ./...). " + purpose
			}
		default:
			continue
		}
		origin := filepath.Join(goDir, filepath.FromSlash(rel))
		copyFile(origin, filepath.Join(p.stage, filepath.FromSlash(rel)))
		inventory = append(inventory, entry{rel, "supplement/" + rel, digest(origin), "unchanged", purpose})
	}
	present := map[string]bool{}
	for _, pkg := range listPackages(goDir) {
		present[pkg] = true
		if !excludedPackages[pkg] && purposes[pkg] == "" {
			fatal("no purpose recorded for Go package %s", pkg)
		}
	}
	for pkg := range purposes {
		if !present[pkg] {
			fatal("purpose recorded for a Go package that does not exist: %s", pkg)
		}
	}

	readme := filepath.Join(p.here, "README.md")
	copyFile(readme, filepath.Join(p.stage, "README.md"))
	inventory = append(inventory, entry{"README.md", "supplement/README.md", digest(readme), "unchanged", "Supplement documentation: how to run the proof runner and what it verifies."})
	licenseFile := filepath.Join(p.root, "LICENSE")
	if info, err := os.Stat(licenseFile); err != nil || !info.Mode().IsRegular() {
		fatal("LICENSE (MIT) is missing from the repository root; it must be shipped in the archive.")
	}
	copyFile(licenseFile, filepath.Join(p.stage, "LICENSE"))
	inventory = append(inventory, entry{"LICENSE", "LICENSE", digest(licenseFile), "unchanged", "MIT License covering all code in this archive (copyright 2026 Tyson Gern)."})

	save(filepath.Join(p.stage, "INPUTS.json"), newObject().set("version", version).set("files", inventory).
		set("scope", fmt.Sprintf("Current-checkout proof snapshot using the %s program suite of https://github.com/tygern/kl; source hashes identify this build.", version)).
		set("provenance", "Programs and notes were developed with assistance from OpenAI and Anthropic models under the direction of Tyson Gern, who takes responsibility for the manuscript. AI-generated review notes are not peer review. The archive README describes the computational checks and implementation dependencies."))
	save(filepath.Join(p.stage, "expected-summary.json"), expected(p.read))

	hashes := newObject()
	staged := listFiles(p.stage)
	for _, rel := range staged {
		hashes.set(rel, digest(filepath.Join(p.stage, filepath.FromSlash(rel))))
	}
	save(filepath.Join(p.stage, "MANIFEST.json"), newObject().set("algorithm", "SHA256").
		set("scope", "Every distributed input file except this manifest itself.").set("sha256", hashes))

	archive := filepath.Join(p.dist, "exceptional-leading-proof.zip")
	writeZip(archive, p.stage, "exceptional-leading-proof/")
	info, err := os.Stat(archive)
	must(err)
	summary := newObject().set("version", version).set("archive", "supplement/dist/exceptional-leading-proof.zip").
		set("bytes", info.Size()).set("sha256", digest(archive)).set("files", len(staged)+1).
		set("staging", "supplement/dist/exceptional-leading-proof").
		set("manuscript_sha256", digest(filepath.Join(p.root, "results", "exceptional-leading.tex")))
	save(filepath.Join(p.dist, "build-summary.json"), summary)
	out, err := encodeJSON(summary)
	must(err)
	fmt.Print(string(out))
}

// repositoryRoot locates the repository: the working directory when it holds
// supplement/go.mod (a workspace invocation), otherwise the parent of the module
// directory (`go run -C supplement ./cmd/package`, or `cd supplement && go run ./cmd/package`).
func repositoryRoot() string {
	cwd, err := os.Getwd()
	must(err)
	for _, candidate := range []string{cwd, filepath.Dir(cwd)} {
		ok := true
		for _, required := range []string{"supplement/go.mod", "supplement/README.md", "results/exceptional-leading.tex"} {
			if _, err := os.Stat(filepath.Join(candidate, filepath.FromSlash(required))); err != nil {
				ok = false
			}
		}
		if ok {
			return candidate
		}
	}
	fatal("Run this command from the repository root as `go run -C supplement ./cmd/package`; supplement/go.mod, supplement/README.md and results/exceptional-leading.tex were not found.")
	return ""
}

// listPackages returns the directories under cmd/ and internal/ holding Go files.
func listPackages(goDir string) []string {
	seen := map[string]bool{}
	for _, rel := range listFiles(goDir) {
		if strings.HasSuffix(rel, ".go") && (strings.HasPrefix(rel, "cmd/") || strings.HasPrefix(rel, "internal/")) {
			seen[filepath.ToSlash(filepath.Dir(rel))] = true
		}
	}
	out := make([]string, 0, len(seen))
	for k := range seen {
		out = append(out, k)
	}
	sort.Strings(out)
	return out
}

// dosDate encodes a calendar date in the MS-DOS format used by ZIP headers.
func dosDate(year int, month time.Month, day int) uint16 {
	return uint16(day + int(month)<<5 + (year-1980)<<9)
}

// writeZip writes every file under stage into a deterministic archive:
// sorted paths under prefix, Deflate at best compression, timestamp
// 2000-01-01 00:00:00 UTC and Unix mode 0644.
func writeZip(archive, stage, prefix string) {
	f, err := os.Create(archive)
	must(err)
	w := zip.NewWriter(f)
	w.RegisterCompressor(zip.Deflate, func(out io.Writer) (io.WriteCloser, error) {
		return flate.NewWriter(out, flate.BestCompression)
	})
	for _, rel := range listFiles(stage) {
		data, err := os.ReadFile(filepath.Join(stage, filepath.FromSlash(rel)))
		must(err)
		header := &zip.FileHeader{Name: prefix + rel, Method: zip.Deflate}
		// A single fixed MS-DOS timestamp (2000-01-01 00:00:00) and no
		// extended-timestamp extra field, as in the v0.2.0 archive; setting
		// Modified instead would add a timezone-displayed extra field.
		header.ModifiedDate, header.ModifiedTime = dosDate(2000, time.January, 1), 0
		header.SetMode(0o644)
		entryWriter, err := w.CreateHeader(header)
		must(err)
		_, err = entryWriter.Write(data)
		must(err)
	}
	must(w.Close())
	must(f.Close())
}
