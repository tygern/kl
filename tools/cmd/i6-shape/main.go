// Command i6-shape certifies the shape of Gern's interval I6 in the Coxeter
// group of type D6 (signed permutations with an even number of sign changes):
// the interval [x, w] with x = (-1,-2,4,3,6,5), w = (-1,-6,3,-4,5,-2) has 1676
// elements and rank 11, but 12 atoms, so it is not a principal lower interval
// (an identity-bottom interval of rank d has at most d atoms), and it is not a
// lattice (two atoms have two incomparable minimal common upper bounds).
//
// It ports two Python programs into one command:
//
//   - research/check_i6_shape.py: uses the "lifting property" Bruhat order
//     (SparseCoxeter.leq of research/archive/python-cpp/computations/sparse_kl.py,
//     repository history at commit d997526).
//     Output: results/i6-shape-certificate.json (same schema, key order and
//     formatting as the original: indent 2).
//   - research/verify_i6_subword.py: an independent check that generates every
//     lower ideal by direct subwords and uses no order comparison and no
//     lifting property. Output: results/i6-independent-subword.json (the one
//     compact JSON line the original printed; compare with
//     results/i6-independent-subword.log).
//
// The two parts live in shape.go and subword.go and share no code except the
// element type and output helpers: each has its own signed-permutation
// arithmetic (length, right multiplication), exactly as the two Python
// programs were independent of each other.
//
// Invocation, from the repository root:
//
//	go run -C tools ./cmd/i6-shape
//
// `go run -C tools` runs with tools/ as the working directory, so the
// repository root is found by walking upward until a directory containing
// results/ and tools/go.mod is reached (the -root flag overrides this).
//
// Imports: standard library only; no other package of this repository.
// Every assertion of the two originals is preserved and aborts with exit
// status 1 and a message when it fails.
package main

import (
	"bytes"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"os"
	"path/filepath"
)

// n is the rank of the Coxeter group D6.
const n = 6

// elem is a signed permutation of 1..6 in window notation, as a Python tuple.
type elem [n]int8

// lessLex compares two elements as Python compares tuples of ints.
func lessLex(a, b elem) bool {
	for i := 0; i < n; i++ {
		if a[i] != b[i] {
			return a[i] < b[i]
		}
	}
	return false
}

// toLists converts elements to [][]int for JSON output.
func toLists(es []elem) [][]int {
	out := make([][]int, len(es))
	for i, e := range es {
		out[i] = make([]int, n)
		for j, v := range e {
			out[i][j] = int(v)
		}
	}
	return out
}

func fail(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "i6-shape: FAILED: "+format+"\n", args...)
	os.Exit(1)
}

func isRoot(dir string) bool {
	if st, err := os.Stat(filepath.Join(dir, "results")); err != nil || !st.IsDir() {
		return false
	}
	st, err := os.Stat(filepath.Join(dir, "tools", "go.mod"))
	return err == nil && !st.IsDir()
}

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
			return "", errors.New("repository root not found: run from the repository root or below it")
		}
		dir = parent
	}
}

// writeJSON writes v to path. indent selects Python's json.dumps(indent=2)
// layout (plus the trailing newline of write_text(... + '\n')) or the compact
// separators=(",", ":") layout (plus the newline of print).
func writeJSON(path string, v any, indent bool) error {
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	if indent {
		enc.SetIndent("", "  ")
	}
	if err := enc.Encode(v); err != nil {
		return err
	}
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	return os.WriteFile(path, buf.Bytes(), 0o644)
}

func main() {
	rootFlag := flag.String("root", "", "repository root (default: found by walking up from the working directory)")
	shapeOut := flag.String("shape-out", "results/i6-shape-certificate.json", "shape certificate, relative to the repository root")
	subwordOut := flag.String("subword-out", "results/i6-independent-subword.json", "subword verification, relative to the repository root")
	flag.Parse()
	root, err := findRoot(*rootFlag)
	if err != nil {
		fail("%v", err)
	}
	shape := computeShape()
	if err := writeJSON(filepath.Join(root, *shapeOut), shape, true); err != nil {
		fail("%v", err)
	}
	fmt.Printf("%s: interval %d, rank %d, atoms %d, common bounds checked %d\n",
		*shapeOut, shape.IntervalSize, shape.Rank, shape.AtomCount, shape.NonLatticeWitness.AllCommonUpperBoundsChecked)
	sub := computeSubword()
	if err := writeJSON(filepath.Join(root, *subwordOut), sub, false); err != nil {
		fail("%v", err)
	}
	fmt.Printf("%s: %s\n", *subwordOut, sub.Status)
}
