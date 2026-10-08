// Command affine-d4 is an independent affine-D4 witness check via subwords
// and R-polynomial reciprocity.
//
// It uses faithful affine signed permutations, no KL descent recurrence and
// no lifting comparison. All lower ideals are direct subword sets. It
// verifies that the Kazhdan-Lusztig polynomial P_{x,w} = 1+3q+2q^2 for the
// four-fork bottom word x = 0134 and the top word w = 0134 2 0134 (labels
// Forks0,1--center2--forks3,4; affine generator 4), by solving the
// reciprocity identity q^d P(q^-1) - P(q) = sum_{x<z<=w} R_{x,z}(q) P_{z,w}(q)
// for every element of the lower ideal of w.
//
// Ports: research/verify_affine_d4_r.py
// Imports: standard library only; no other engine and no internal package.
//
// Usage (inside the work directory): affine-d4
// Writes results/affine-d4-independent.json and prints the same JSON.
package main

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

// elem is an affine signed permutation: a signed permutation of 1..4 and a
// translation vector.
type elem struct {
	P [4]int
	T [4]int
}

type poly []int64

var identity = elem{P: [4]int{1, 2, 3, 4}}

var qx = []int{0, 1, 3, 4}
var qw = []int{0, 1, 3, 4, 2, 0, 1, 3, 4}

func abs(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

func simple(w elem, s int) elem {
	p, t := w.P, w.T
	switch {
	case s == 0:
		p[0], p[1] = -p[1], -p[0]
	case s < 4:
		p[s-1], p[s] = p[s], p[s-1]
	default:
		for _, a := range p[2:] {
			if a > 0 {
				t[abs(a)-1]++
			} else {
				t[abs(a)-1]--
			}
		}
		p[2], p[3] = -p[3], -p[2]
	}
	return elem{P: p, T: t}
}

func evaluate(word []int) elem {
	w := identity
	for _, s := range word {
		w = simple(w, s)
	}
	return w
}

// subwords returns every subword, in the order of itertools.product((0,1)).
func subwords(word []int) [][]int {
	n := len(word)
	out := make([][]int, 0, 1<<uint(n))
	for m := 0; m < 1<<uint(n); m++ {
		sub := make([]int, 0, n)
		for i := 0; i < n; i++ {
			if m>>uint(n-1-i)&1 == 1 {
				sub = append(sub, word[i])
			}
		}
		out = append(out, sub)
	}
	return out
}

func trim(a []int64) poly {
	for len(a) > 0 && a[len(a)-1] == 0 {
		a = a[:len(a)-1]
	}
	return poly(a)
}

// addScaled adds scale*c to result[k] with int64 overflow guards.
func addScaled(result []int64, k int, c int64, scale int64) {
	prod := scale * c
	if scale != 0 && prod/scale != c {
		panic("int64 overflow")
	}
	sum := result[k] + prod
	if (prod > 0 && sum < result[k]) || (prod < 0 && sum > result[k]) {
		panic("int64 overflow")
	}
	result[k] = sum
}

func add(a, b []int64, shift int, scale int64) poly {
	n := len(a)
	if len(b)+shift > n {
		n = len(b) + shift
	}
	result := make([]int64, n)
	copy(result, a)
	for k, c := range b {
		addScaled(result, k+shift, c, scale)
	}
	return trim(result)
}

func multiply(a, b poly) poly {
	result := poly{}
	for k, c := range b {
		result = add(result, a, k, c)
	}
	return result
}

func polyEqual(a, b poly) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func floorDiv(a, b int) int {
	q := a / b
	if (a%b != 0) && ((a < 0) != (b < 0)) {
		q--
	}
	return q
}

type pair struct{ x, w elem }

func fail(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "affine-d4: assertion failed: "+format+"\n", args...)
	os.Exit(1)
}

func wordString(w []int) string {
	parts := make([]string, len(w))
	for i, s := range w {
		parts[i] = strconv.Itoa(s)
	}
	return "(" + strings.Join(parts, ",") + ")"
}

type rankCounts struct {
	keys   []int
	counts map[int]int
}

func (r rankCounts) MarshalJSON() ([]byte, error) {
	var sb strings.Builder
	sb.WriteByte('{')
	for i, k := range r.keys {
		if i > 0 {
			sb.WriteByte(',')
		}
		fmt.Fprintf(&sb, "\"%d\":%d", k, r.counts[k])
	}
	sb.WriteByte('}')
	return []byte(sb.String()), nil
}

type report struct {
	Method           string     `json:"method"`
	SimpleLabels     string     `json:"simple_labels"`
	BottomWord       []int      `json:"bottom_word"`
	TopWord          []int      `json:"top_word"`
	Bottom           [2][]int   `json:"bottom"`
	Top              [2][]int   `json:"top"`
	Lengths          [2]int     `json:"lengths"`
	LowerIdealSize   int        `json:"lower_ideal_size"`
	IntervalRank     rankCounts `json:"interval_rank_counts"`
	Polynomial       []int64    `json:"polynomial"`
	Mu               int        `json:"mu"`
	FCBottom         string     `json:"fc_bottom"`
	ReciprocityCount int        `json:"reciprocity_identities_checked"`
	Status           string     `json:"status"`
}

func main() {
	// words: shortest subword (first found in product order) of each element.
	words := map[elem][]int{}
	var order []elem // insertion order, as Python dict
	for _, word := range subwords(qw) {
		w := evaluate(word)
		cur, ok := words[w]
		if !ok {
			order = append(order, w)
			words[w] = word
		} else if len(word) < len(cur) {
			words[w] = word
		}
	}
	lengths := map[elem]int{}
	for w, word := range words {
		lengths[w] = len(word)
	}
	ideals := map[elem]map[elem]bool{}
	for w, q := range words {
		set := map[elem]bool{}
		for _, sw := range subwords(q) {
			set[evaluate(sw)] = true
		}
		ideals[w] = set
	}
	below := func(x, w elem) bool {
		return ideals[w][x] // nil map lookup is false for unknown w
	}
	lenOr100 := func(e elem) int {
		if l, ok := lengths[e]; ok {
			return l
		}
		return 100
	}
	desc := func(w elem) []int {
		var d []int
		for s := 0; s < 5; s++ {
			if lenOr100(simple(w, s)) < lengths[w] {
				d = append(d, s)
			}
		}
		return d
	}
	cache := map[pair]poly{}
	var r func(x, w elem) poly
	r = func(x, w elem) poly {
		key := pair{x, w}
		if v, ok := cache[key]; ok {
			return v
		}
		var res poly
		switch {
		case !below(x, w):
			res = poly{}
		case x == w:
			res = poly{1}
		default:
			s := desc(w)[0]
			xs, ws := simple(x, s), simple(w, s)
			if lenOr100(xs) < lengths[x] {
				res = r(xs, ws)
			} else {
				rxw := r(x, ws)
				res = add(add(poly{}, rxw, 0, -1), add(rxw, r(xs, ws), 0, 1), 1, 1)
			}
		}
		cache[key] = res
		return res
	}
	w := evaluate(qw)
	x := evaluate(qx)
	if !(lengths[x] == 4 && lengths[w] == 9) {
		fail("lengths[x]=%d lengths[w]=%d", lengths[x], lengths[w])
	}
	p := map[elem]poly{w: {1}}
	// q^d P(q^-1) - P(q) = sum_{x<z<=w} R_{x,z}(q) P_{z,w}(q).
	ys := append([]elem(nil), order...)
	sort.SliceStable(ys, func(i, j int) bool { return lengths[ys[i]] > lengths[ys[j]] })
	for _, y := range ys {
		if y == w {
			continue
		}
		d := lengths[w] - lengths[y]
		rhs := poly{}
		for _, z := range order {
			if lengths[z] > lengths[y] && below(y, z) {
				rhs = add(rhs, multiply(r(y, z), p[z]), 0, 1)
			}
		}
		bound := floorDiv(d-1, 2)
		value := make([]int64, 0, bound+1)
		for k := 0; k <= bound; k++ {
			if k < len(rhs) {
				value = append(value, -rhs[k])
			} else {
				value = append(value, 0)
			}
		}
		val := trim(value)
		p[y] = val
		dual := make([]int64, d+1)
		for k, c := range val {
			dual[d-k] += c
		}
		if !polyEqual(add(dual, val, 0, -1), rhs) {
			fail("reciprocity %s value=%v rhs=%v", wordString(words[y]), []int64(val), []int64(rhs))
		}
		if len(val) == 0 || val[0] != 1 {
			fail("P(0) != 1 for %s", wordString(words[y]))
		}
		for _, c := range val {
			if c < 0 {
				fail("negative coefficient for %s", wordString(words[y]))
			}
		}
	}
	px := p[x]
	if !polyEqual(px, poly{1, 3, 2}) {
		fail("p[x] = %v, expected (1,3,2)", []int64(px))
	}
	counts := map[int]int{}
	for z := range ideals[w] {
		if below(x, z) {
			counts[lengths[z]-lengths[x]]++
		}
	}
	var keys []int
	for k := range counts {
		keys = append(keys, k)
	}
	sort.Ints(keys)

	rep := report{
		Method:           "Affine signed actions, direct subwords, R-polynomial reciprocity",
		SimpleLabels:     "Forks0,1--center2--forks3,4; affine generator4",
		BottomWord:       qx,
		TopWord:          qw,
		Bottom:           [2][]int{x.P[:], x.T[:]},
		Top:              [2][]int{w.P[:], w.T[:]},
		Lengths:          [2]int{lengths[x], lengths[w]},
		LowerIdealSize:   len(words),
		IntervalRank:     rankCounts{keys: keys, counts: counts},
		Polynomial:       []int64(px),
		Mu:               2,
		FCBottom:         "Four pairwise commuting generators, each once",
		ReciprocityCount: len(p) - 1,
		Status:           "passed",
	}
	out, err := json.MarshalIndent(rep, "", "  ")
	if err != nil {
		fail("marshal: %v", err)
	}
	out = append(out, '\n')
	path := filepath.Join("results", "affine-d4-independent.json")
	if err := os.MkdirAll("results", 0o755); err != nil {
		fail("mkdir: %v", err)
	}
	if err := os.WriteFile(path, out, 0o644); err != nil {
		fail("write: %v", err)
	}
	os.Stdout.Write(out)
}
