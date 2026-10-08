// Command catalogue-eligible lists, for every full-support terminal element w
// found by the conjugate search, the fully commutative (FC) elements x <= w
// of the complete FC catalogue that could be Kazhdan-Lusztig-relevant (the
// "eligible" elements): x lies in the Bruhat interval below w and has at
// least the right and the left descents of w.
//
// It ports research/en_families/catalogue_eligible.py to Go (the Coxeter
// arithmetic used there is research/en_families/targeted.py, class Coxeter and
// en(n), reimplemented below). Elements of the Coxeter group of type E_n (the
// generalized, possibly infinite group with diagram 0-1-...-(n-2), and node
// n-1 attached to node 2) are stored as integer matrices by columns in the
// basis of simple roots; all arithmetic is exact (int64 with an overflow
// guard; entries are root coordinates, at most 7 in the committed data).
//
// Inputs (paths relative to the repository root; `go run -C tools` runs with
// tools/ as the working directory, so the root is found by walking upward to
// the first directory containing results/ and tools/go.mod; -root overrides):
//
//	research/en_families/conjugates_E{n}_50000_934.json   committed (-conjugates overrides)
//	research/en_independent/e{n}-fc.json                  the complete FC catalogue
//
// The catalogue files research/en_independent/e{n}-fc.json are gitignored
// outputs of the FC-catalogue enumerator; generate them first with
//
//	go run -C supplement ./cmd/fc-catalogue -rank N > research/en_independent/eN-fc.json
//
// (the enumerator supports N = 6..10; ranks 11 and 12 have no catalogue, so
// this command cannot be run for them). The committed conjugate-search file for
// rank 8 has 20000 trials, not 50000: use -conjugates for it.
//
// Output: research/en_families/e{n}_full_support_eligible.json (-out
// overrides), byte-identical to the Python original's json.dumps(indent=2)
// output; one JSON line per full-support row is also printed to standard
// output (the original printed the Python dict repr, a progress log only).
// Every assertion of the original (w terminal, x fully commutative, the
// catalogue masks equal the computed descent masks) is an explicit check that
// exits non-zero with a message.
//
// Invocation (from the repository root or from tools/):
//
//	go run -C tools ./cmd/catalogue-eligible -rank 9
//
// Imports: standard library only; no other package of this repository.
package main

import (
	"bytes"
	"encoding/binary"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"io/fs"
	"os"
	"path/filepath"
	"strings"
)

type fatal struct{ msg string }

func failf(format string, args ...any) {
	panic(fatal{fmt.Sprintf(format, args...)})
}

func check(ok bool, format string, args ...any) {
	if !ok {
		failf("assertion failed: "+format, args...)
	}
}

// elt is a group element as columns: elt[j] is the image of simple root j.
type elt [][]int64

const coordLimit = int64(1) << 40

// coxeter is the Coxeter group of E_n with the matrix model of targeted.py.
type coxeter struct {
	n        int
	adj      [][]int
	identity elt
	fcMemo   map[string]bool
}

func newEn(n int) *coxeter {
	check(n >= 4 && n <= 30, "rank %d outside 4..30", n)
	g := &coxeter{n: n, adj: make([][]int, n), fcMemo: map[string]bool{}}
	var edges [][2]int
	for i := 0; i < n-2; i++ {
		edges = append(edges, [2]int{i, i + 1})
	}
	edges = append(edges, [2]int{2, n - 1})
	for _, e := range edges {
		g.adj[e[0]] = append(g.adj[e[0]], e[1])
		g.adj[e[1]] = append(g.adj[e[1]], e[0])
	}
	g.identity = make(elt, n)
	for j := range g.identity {
		g.identity[j] = make([]int64, n)
		g.identity[j][j] = 1
	}
	return g
}

func (g *coxeter) key(w elt) string {
	var b bytes.Buffer
	for _, col := range w {
		for _, v := range col {
			_ = binary.Write(&b, binary.LittleEndian, v)
		}
	}
	return b.String()
}

func (g *coxeter) equal(a, b elt) bool {
	for j := range a {
		for i := range a[j] {
			if a[j][i] != b[j][i] {
				return false
			}
		}
	}
	return true
}

// right returns w s: column s negated, neighbours t get w[s] + w[t].
func (g *coxeter) right(w elt, s int) elt {
	out := make(elt, g.n)
	copy(out, w)
	neg := make([]int64, g.n)
	for i, v := range w[s] {
		neg[i] = -v
	}
	out[s] = neg
	for _, t := range g.adj[s] {
		col := make([]int64, g.n)
		for i := range col {
			v := w[s][i] + w[t][i]
			check(v > -coordLimit && v < coordLimit, "matrix entry overflow guard")
			col[i] = v
		}
		out[t] = col
	}
	return out
}

// descMask returns the bit mask of right descents (columns with all entries <= 0).
func (g *coxeter) descMask(w elt) int {
	m := 0
	for i, col := range w {
		all := true
		for _, v := range col {
			if v > 0 {
				all = false
				break
			}
		}
		if all {
			m |= 1 << i
		}
	}
	return m
}

func (g *coxeter) desc(w elt) []int {
	var out []int
	m := g.descMask(w)
	for i := 0; i < g.n; i++ {
		if m>>i&1 == 1 {
			out = append(out, i)
		}
	}
	return out
}

// word returns the reduced word obtained by repeatedly removing the smallest
// right descent (the original's recursive word()).
func (g *coxeter) word(w elt) []int {
	var rev []int
	for !g.equal(w, g.identity) {
		ds := g.desc(w)
		check(len(ds) > 0, "element without descent")
		rev = append(rev, ds[0])
		w = g.right(w, ds[0])
	}
	out := make([]int, len(rev))
	for i, s := range rev {
		out[len(rev)-1-i] = s
	}
	return out
}

func (g *coxeter) length(w elt) int { return len(g.word(w)) }

func (g *coxeter) elt(word []int) elt {
	w := g.identity
	for _, s := range word {
		check(s >= 0 && s < g.n, "letter %d out of range for rank %d", s, g.n)
		w = g.right(w, s)
	}
	return w
}

func (g *coxeter) inv(w elt) elt {
	word := g.word(w)
	out := g.identity
	for i := len(word) - 1; i >= 0; i-- {
		out = g.right(out, word[i])
	}
	return out
}

func (g *coxeter) fc(w elt) bool {
	k := g.key(w)
	if v, ok := g.fcMemo[k]; ok {
		return v
	}
	mask := g.descMask(w)
	ds := g.desc(w)
	res := true
	for _, s := range ds {
		for _, t := range g.adj[s] {
			if mask>>t&1 == 1 {
				res = false
			}
		}
	}
	if res {
		for _, s := range ds {
			if !g.fc(g.right(w, s)) {
				res = false
				break
			}
		}
	}
	g.fcMemo[k] = res
	return res
}

func (g *coxeter) terminal(w elt) bool {
	for _, v := range []elt{w, g.inv(w)} {
		for _, s := range g.desc(v) {
			m := g.descMask(g.right(v, s))
			for _, t := range g.adj[s] {
				if m>>t&1 == 1 {
					return false
				}
			}
		}
	}
	return true
}

// leq is the Bruhat order test of the original (subword property via the
// smallest right descent of w; lengths are tracked instead of recomputed).
func (g *coxeter) leq(x, w elt) bool {
	lx, lw := g.length(x), g.length(w)
	for !g.equal(x, w) {
		if lx >= lw {
			return false
		}
		s := g.desc(w)[0]
		if g.descMask(x)>>s&1 == 1 {
			x = g.right(x, s)
			lx--
		}
		w = g.right(w, s)
		lw--
	}
	return true
}

// trace is the sum of the diagonal entries.
func trace(w elt) int64 {
	var t int64
	for i := range w {
		t += w[i][i]
	}
	return t
}

// floorDiv2 is Python's v // 2.
func floorDiv2(v int64) int64 {
	if v >= 0 {
		return v / 2
	}
	return -((-v + 1) / 2)
}

type conjRow struct {
	Rank        int   `json:"rank"`
	Word        []int `json:"word"`
	Length      int   `json:"length"`
	Support     []int `json:"support"`
	R           []int `json:"R"`
	FullSupport bool  `json:"full_support"`
	Trial       int   `json:"trial"`
}

type eligibleItem struct {
	Word   []int `json:"word"`
	Length int   `json:"length"`
	Rank   int   `json:"rank"`
}

type outRow struct {
	conjRow
	ReflectionRank int64          `json:"reflection_rank"`
	MaskCandidates int            `json:"mask_candidates"`
	Eligible       []eligibleItem `json:"eligible"`
}

// catalogueFile is the part of the FC catalogue used here.
type catalogueFile struct {
	Elements []struct {
		Word  []int `json:"word"`
		Rmask int   `json:"Rmask"`
		Lmask int   `json:"Lmask"`
	} `json:"elements"`
}

// reason describes a file error without echoing the (absolute) path.
func reason(err error) string {
	var pe *fs.PathError
	if errors.As(err, &pe) {
		return pe.Err.Error()
	}
	return err.Error()
}

func loadConjugates(path string) []conjRow {
	data, err := os.ReadFile(path)
	if err != nil {
		failf("cannot read %s: %s (for rank 8 the committed file has 20000 trials; use -conjugates)", filepath.Base(path), reason(err))
	}
	var top struct {
		Terminals []json.RawMessage `json:"terminals"`
	}
	check(json.Unmarshal(data, &top) == nil, "%s: invalid JSON", path)
	rows := make([]conjRow, 0, len(top.Terminals))
	for _, raw := range top.Terminals {
		var r conjRow
		dec := json.NewDecoder(bytes.NewReader(raw))
		dec.DisallowUnknownFields()
		if err := dec.Decode(&r); err != nil {
			failf("%s: unexpected terminal row: %v", path, err)
		}
		rows = append(rows, r)
	}
	return rows
}

func loadCatalogue(path string, n int) *catalogueFile {
	f, err := os.Open(path)
	if err != nil {
		failf("cannot read %s: %s\nthe FC catalogue is a gitignored output of fc-catalogue; generate it with (from the repository root):\n  go run -C supplement ./cmd/fc-catalogue -rank %d > research/en_independent/e%d-fc.json\n(the enumerator supports ranks 6..10 only)", filepath.Base(path), reason(err), n, n)
	}
	defer f.Close()
	var c catalogueFile
	if err := json.NewDecoder(f).Decode(&c); err != nil {
		failf("%s: invalid catalogue: %v", path, err)
	}
	return &c
}

func run(n int, conjPath, cataloguePath, outPath string, log io.Writer) {
	rows := loadConjugates(conjPath)
	catalogue := loadCatalogue(cataloguePath, n)
	g := newEn(n)
	out := []outRow{}
	for _, row := range rows {
		if !row.FullSupport {
			continue
		}
		w := g.elt(row.Word)
		check(g.terminal(w), "conjugate-search element is not terminal (word %v)", row.Word)
		rm := g.descMask(w)
		lm := g.descMask(g.inv(w))
		matches := 0
		eligible := []eligibleItem{}
		lw := g.length(w)
		for _, item := range catalogue.Elements {
			if item.Rmask&rm != rm || item.Lmask&lm != lm {
				continue
			}
			matches++
			x := g.elt(item.Word)
			check(g.fc(x), "catalogue element %v is not fully commutative", item.Word)
			check(g.descMask(x) == item.Rmask, "right-descent mask mismatch for %v", item.Word)
			check(g.descMask(g.inv(x)) == item.Lmask, "left-descent mask mismatch for %v", item.Word)
			if g.leq(x, w) {
				lx := g.length(x)
				eligible = append(eligible, eligibleItem{Word: item.Word, Length: lx, Rank: lw - lx})
			}
		}
		o := outRow{conjRow: row, ReflectionRank: floorDiv2(int64(n) - trace(w)), MaskCandidates: matches, Eligible: eligible}
		line, err := json.Marshal(o)
		check(err == nil, "marshal")
		fmt.Fprintln(log, string(line))
		out = append(out, o)
	}
	data, err := json.MarshalIndent(out, "", "  ")
	check(err == nil, "marshal output")
	data = append(data, '\n')
	check(os.WriteFile(outPath, data, 0o644) == nil, "cannot write %s", outPath)
}

// findRoot walks upward from the working directory to the first directory that
// contains both results/ and tools/go.mod.
func findRoot() string {
	dir, err := os.Getwd()
	check(err == nil, "cannot determine working directory")
	for {
		_, e1 := os.Stat(filepath.Join(dir, "results"))
		_, e2 := os.Stat(filepath.Join(dir, "tools", "go.mod"))
		if e1 == nil && e2 == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		check(parent != dir, "repository root not found (needs results/ and tools/go.mod); use -root")
		dir = parent
	}
}

func main() {
	root := flag.String("root", "", "repository root (default: found by walking up from the working directory)")
	rank := flag.Int("rank", 0, "rank n of E_n (required; catalogues exist for 6..10)")
	conj := flag.String("conjugates", "", "conjugate-search file (default research/en_families/conjugates_E{n}_50000_934.json)")
	cat := flag.String("catalogue", "", "FC catalogue (default research/en_independent/e{n}-fc.json; gitignored, generate with fc-catalogue)")
	out := flag.String("out", "", "output file (default research/en_families/e{n}_full_support_eligible.json)")
	flag.Parse()
	if flag.NArg() != 0 || *rank == 0 {
		fmt.Fprintln(os.Stderr, "usage: catalogue-eligible -rank N [-root DIR] [-conjugates F] [-catalogue F] [-out F]")
		os.Exit(2)
	}
	n := *rank
	defer func() {
		if r := recover(); r != nil {
			if f, ok := r.(fatal); ok {
				fmt.Fprintln(os.Stderr, "catalogue-eligible: "+strings.TrimSpace(f.msg))
				os.Exit(1)
			}
			panic(r)
		}
	}()
	if *root == "" {
		*root = findRoot()
	}
	rel := func(p string) string { return filepath.Join(*root, p) }
	if *conj == "" {
		*conj = rel(fmt.Sprintf("research/en_families/conjugates_E%d_50000_934.json", n))
	}
	if *cat == "" {
		*cat = rel(fmt.Sprintf("research/en_independent/e%d-fc.json", n))
	}
	if *out == "" {
		*out = rel(fmt.Sprintf("research/en_families/e%d_full_support_eligible.json", n))
	}
	run(n, *conj, *cat, *out, os.Stdout)
}
