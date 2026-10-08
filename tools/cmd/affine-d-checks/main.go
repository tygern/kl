// Command affine-d-checks regenerates two exploratory certificates of the
// broad affine-D study.
//
//  1. research/broad_affine/unfolding_obstruction.json: the lower Bruhat
//     interval of the length-10 element with reduced word (0,1,3,4,2)^2 in
//     the affine Weyl group of type D4, computed in the faithful affine
//     signed-permutation model; its rank vector, and the four rank-3
//     two-generator-support edges that obstruct a principal finite-D
//     unfolding (plus a check of the Coxeter relations in the model).
//  2. research/broad_affine/gern_coatoms.json: for the Gern bad pairs
//     (x, w) in finite D_n, n = 4, 6, ..., 24, the lower covers of w by all
//     signed transpositions and those lying above x in Bruhat order.
//
// It ports research/broad_affine/check_unfolding.py (together with the class
// AffineD of research/broad_affine/affine_d_crowns.py) and
// research/broad_affine/gern_coatoms.py (together with the Bruhat order,
// length formula and right action of SparseCoxeter in
// research/archive/python-cpp/computations/sparse_kl.py, repository history
// at commit d997526; the Kazhdan-Lusztig
// parts of that class are never used by gern_coatoms.py and are not ported).
// Output files are written with the same layout as Python's json.dump with
// indent=2 and are byte-identical to the originals. All arithmetic is exact
// small-integer arithmetic.
//
// Invocation (from the repository root; "go run -C tools" runs the program
// inside tools/, so the repository root is found as the nearest ancestor of
// the working directory that holds results/ and tools/go.mod; -root overrides):
//
//	go run -C tools ./cmd/affine-d-checks
//
// Imports: standard library only; no other package of this repository.
package main

import (
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"path/filepath"
)

const (
	unfoldingPath = "research/broad_affine/unfolding_obstruction.json"
	coatomsPath   = "research/broad_affine/gern_coatoms.json"
)

func fail(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "assertion failed: "+format+"\n", args...)
	os.Exit(1)
}

// writeJSON writes v in the layout of Python's json.dump(v, f, indent=2):
// two-space indentation, no trailing newline.
func writeJSON(path string, v any) {
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	enc.SetIndent("", "  ")
	if err := enc.Encode(v); err != nil {
		fail("encode %s: %v", path, err)
	}
	out := bytes.TrimRight(buf.Bytes(), "\n")
	if err := os.WriteFile(path, out, 0o644); err != nil {
		fail("write %s: %v", path, err)
	}
}

// isRoot reports whether dir is the repository root.
func isRoot(dir string) bool {
	if st, err := os.Stat(filepath.Join(dir, "results")); err != nil || !st.IsDir() {
		return false
	}
	st, err := os.Stat(filepath.Join(dir, "tools", "go.mod"))
	return err == nil && st.Mode().IsRegular()
}

func findRoot(explicit string) string {
	start := explicit
	if start == "" {
		wd, err := os.Getwd()
		if err != nil {
			fail("%v", err)
		}
		start = wd
	}
	start, err := filepath.Abs(start)
	if err != nil {
		fail("%v", err)
	}
	if explicit != "" {
		if !isRoot(start) {
			fail("%s is not a repository root (needs results/ and tools/go.mod)", explicit)
		}
		return start
	}
	for dir := start; ; {
		if isRoot(dir) {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			fail("repository root not found (no ancestor has results/ and tools/go.mod); use -root")
		}
		dir = parent
	}
}

func main() {
	rootFlag := flag.String("root", "", "repository root (default: nearest ancestor of the working directory with results/ and tools/go.mod)")
	flag.Parse()
	if flag.NArg() > 0 {
		fmt.Fprintln(os.Stderr, "usage: affine-d-checks [-root dir]")
		os.Exit(2)
	}
	root := findRoot(*rootFlag)
	unfolding := runUnfolding()
	writeJSON(filepath.Join(root, unfoldingPath), unfolding)
	fmt.Printf("wrote %s: interval size %d, rank vector %v, edges %v\n",
		unfoldingPath, unfolding.IntervalSize, unfolding.RankVector, unfolding.DetectedRank3Edges)

	rows := runGernCoatoms(4, 24)
	writeJSON(filepath.Join(root, coatomsPath), rows)
	fmt.Printf("wrote %s: %d rows\n", coatomsPath, len(rows))
}
