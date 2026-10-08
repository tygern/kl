// Command verify-outputs cross-checks all terminal outputs of the E6/E7/E8
// computations using exact integer matrices.
//
// It does not enumerate a Coxeter group. For every certificate it verifies
// reduced word lengths, descent masks, both terminal conditions (weak-order
// condition on w and w^-1, and the commuting-bottom condition), and the
// equivalence of the outputs obtained from the flat E7 cosets, the recursive
// cosets, and the independent D7 implementation.
//
// Ports research/en_e8/verify_outputs.py. Standard library only; it imports
// none of the other engines or internal packages (self-contained). It reads
// its inputs relative to the current working directory (the work directory).
//
// Usage:
//
//	verify-outputs             print the finite method comparison summary (JSON)
//	verify-outputs -matrices N check research/en_e8/eN-recursive.json and print
//	                           the generated terminal matrices as a JSON list
//	                           of N x N integer matrices. Convention: a matrix
//	                           is a list of N columns; entry [s] is the image
//	                           vector of the simple root alpha_s (the same
//	                           "tuple of columns" layout as the original). The
//	                           list is sorted lexicographically.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"
)

func fail(format string, args ...any) {
	panic(fmt.Sprintf("assertion failed: "+format, args...))
}

func assert(cond bool, format string, args ...any) {
	if !cond {
		fail(format, args...)
	}
}

type matrix [][]int64

func (m matrix) key() string {
	var sb strings.Builder
	for _, col := range m {
		for _, c := range col {
			sb.WriteString(strconv.FormatInt(c, 10))
			sb.WriteByte(',')
		}
		sb.WriteByte(';')
	}
	return sb.String()
}

type word []int64

type bottom struct {
	Word         word   `json:"word"`
	Length       *int64 `json:"length"`
	Rank         *int64 `json:"rank"`
	IntervalRank *int64 `json:"interval_rank"`
}

type badEntry struct {
	Word          word     `json:"word"`
	Length        *int64   `json:"length"`
	Rmask         *int64   `json:"Rmask"`
	Lmask         *int64   `json:"Lmask"`
	RightDescents []int64  `json:"right_descents"`
	LeftDescents  []int64  `json:"left_descents"`
	Bottoms       []bottom `json:"bottoms"`
	Candidates    []bottom `json:"candidates"`
}

type certificate struct {
	Bad                json.RawMessage `json:"bad"`
	NonFCWeakBad       json.RawMessage `json:"non_fc_weak_bad_terminals"`
	RightTerminalCount *int64          `json:"right_terminal_count"`
	CommutingTerminals *int64          `json:"commuting_terminals"`
	FullRightTerminals *int64          `json:"full_right_terminals"`
	ParabolicHistogram []int64         `json:"parabolic_histogram"`
	CosetHistogram     []int64         `json:"coset_histogram"`
	hasParabolicHist   bool
	hasCosetHist       bool
}

// entry is the value stored per terminal matrix w: (length, Rmask, Lmask, set
// of (matrix, length, gap)). The bottoms are kept as a sorted, deduplicated
// list of canonical strings so that entries compare by value.
type entry struct {
	length  int64
	rd, ld  int64
	bottoms []string
}

func (e entry) equal(o entry) bool {
	if e.length != o.length || e.rd != o.rd || e.ld != o.ld || len(e.bottoms) != len(o.bottoms) {
		return false
	}
	for i := range e.bottoms {
		if e.bottoms[i] != o.bottoms[i] {
			return false
		}
	}
	return true
}

type result map[string]entry // keyed by matrix.key()

func resultsEqual(a, b result) bool {
	if len(a) != len(b) {
		return false
	}
	for k, v := range a {
		w, ok := b[k]
		if !ok || !v.equal(w) {
			return false
		}
	}
	return true
}

func readFile(path string) []byte {
	data, err := os.ReadFile(path)
	if err != nil {
		panic(fmt.Sprintf("cannot read %s: %v", path, err))
	}
	return data
}

// checkCertificate mirrors check_certificate of the original. It returns the
// parsed certificate, the result keyed by matrix key, and the matrices.
func checkCertificate(path string, n int) (*certificate, result, map[string]matrix) {
	raw := readFile(path)
	var data certificate
	if err := json.Unmarshal(raw, &data); err != nil {
		panic(fmt.Sprintf("cannot parse %s: %v", path, err))
	}
	// Presence of the histogram keys (the original indexes them with [] only
	// for the D7 certificate); recorded so callers can assert presence.
	var probe map[string]json.RawMessage
	if err := json.Unmarshal(raw, &probe); err != nil {
		panic(fmt.Sprintf("cannot parse %s: %v", path, err))
	}
	_, data.hasParabolicHist = probe["parabolic_histogram"]
	_, data.hasCosetHist = probe["coset_histogram"]

	type edge struct{ a, b int }
	var edges []edge
	for i := 0; i < n-2; i++ {
		edges = append(edges, edge{i, i + 1})
	}
	edges = append(edges, edge{2, n - 1})
	adj := make([][]int, n)
	has := make([][]bool, n)
	for i := range has {
		has[i] = make([]bool, n)
	}
	for _, e := range edges {
		if !has[e.a][e.b] {
			has[e.a][e.b] = true
			adj[e.a] = append(adj[e.a], e.b)
		}
		if !has[e.b][e.a] {
			has[e.b][e.a] = true
			adj[e.b] = append(adj[e.b], e.a)
		}
	}

	positive := func(v []int64) bool {
		any := false
		allNonNeg, allNonPos := true, true
		for _, c := range v {
			if c != 0 {
				any = true
			}
			if c < 0 {
				allNonNeg = false
			}
			if c > 0 {
				allNonPos = false
			}
		}
		assert(any, "zero vector")
		assert(allNonNeg || allNonPos, "mixed-sign root %v", v)
		return allNonNeg
	}
	add := func(x, y []int64) []int64 {
		out := make([]int64, len(x))
		for i := range x {
			out[i] = x[i] + y[i]
			if out[i] > 1<<40 || out[i] < -(1<<40) {
				fail("coefficient overflow guard")
			}
		}
		return out
	}
	element := func(w word) matrix {
		a := make(matrix, n)
		for j := 0; j < n; j++ {
			a[j] = make([]int64, n)
			a[j][j] = 1
		}
		for _, s64 := range w {
			assert(s64 >= 0 && int(s64) < n, "%s letter %d out of range", path, s64)
			s := int(s64)
			assert(positive(a[s]), "%s nonreduced word %v", path, w)
			b := make(matrix, n)
			copy(b, a)
			neg := make([]int64, n)
			for i, c := range a[s] {
				neg[i] = -c
			}
			b[s] = neg
			for _, t := range adj[s] {
				b[t] = add(a[t], a[s])
			}
			a = b
		}
		return a
	}
	desc := func(a matrix) int64 {
		var m int64
		for s := 0; s < n; s++ {
			if !positive(a[s]) {
				m += 1 << uint(s)
			}
		}
		return m
	}
	weak := func(a matrix) {
		for s := 0; s < n; s++ {
			if !positive(a[s]) {
				for _, t := range adj[s] {
					assert(positive(add(a[s], a[t])), "%s weak-order condition fails at s=%d t=%d", path, s, t)
				}
			}
		}
	}
	reversed := func(w word) word {
		r := make(word, len(w))
		for i, c := range w {
			r[len(w)-1-i] = c
		}
		return r
	}
	maskOf := func(l []int64) int64 {
		var m int64
		for _, s := range l {
			m += 1 << uint(s)
		}
		return m
	}
	isEdge := func(s, t int64) bool {
		for _, e := range edges {
			if (int64(e.a) == s && int64(e.b) == t) || (int64(e.a) == t && int64(e.b) == s) {
				return true
			}
		}
		return false
	}

	// bads = data.get("bad", data.get("non_fc_weak_bad_terminals"))
	badRaw := data.Bad
	if _, ok := probe["bad"]; !ok {
		badRaw = data.NonFCWeakBad
	}
	if len(badRaw) == 0 || string(badRaw) == "null" {
		panic(fmt.Sprintf("%s: no bad terminal list", path))
	}
	var bads []badEntry
	if err := json.Unmarshal(badRaw, &bads); err != nil {
		panic(fmt.Sprintf("%s: cannot parse bad list: %v", path, err))
	}
	// Detect per-entry presence of "bottoms" (vs "candidates").
	var rawBads []map[string]json.RawMessage
	if err := json.Unmarshal(badRaw, &rawBads); err != nil {
		panic(fmt.Sprintf("%s: cannot parse bad list: %v", path, err))
	}

	res := result{}
	mats := map[string]matrix{}
	for bi, b := range bads {
		w := element(b.Word)
		wi := element(reversed(b.Word))
		assert(b.Length != nil && int64(len(b.Word)) == *b.Length, "word length mismatch")
		weak(w)
		weak(wi)
		var rd, ld int64
		if b.Rmask != nil {
			rd = *b.Rmask
		} else {
			rd = maskOf(b.RightDescents)
		}
		if b.Lmask != nil {
			ld = *b.Lmask
		} else {
			ld = maskOf(b.LeftDescents)
		}
		assert(rd == desc(w) && ld == desc(wi), "descent masks mismatch for %v", b.Word)
		xs := b.Candidates
		if _, ok := rawBads[bi]["bottoms"]; ok {
			xs = b.Bottoms
		}
		wsupport := map[int64]bool{}
		for _, s := range b.Word {
			wsupport[s] = true
		}
		set := map[string]bool{}
		for _, x := range xs {
			a := element(x.Word)
			ai := element(reversed(x.Word))
			assert(x.Length != nil && int64(len(x.Word)) == *x.Length, "bottom length mismatch")
			assert(rd&desc(a) == rd && ld&desc(ai) == ld, "bottom does not dominate descents")
			support := map[int64]bool{}
			for _, s := range x.Word {
				support[s] = true
			}
			assert(len(support) == len(x.Word), "repeated letter in bottom")
			for s := range support {
				for t := range support {
					assert(!isEdge(s, t), "non-commuting bottom support")
				}
			}
			// A product of distinct commuting generators is below w whenever
			// its support is contained in w: select one of each as a subword.
			for s := range support {
				assert(wsupport[s], "bottom support not contained in word support")
			}
			gap := *b.Length - *x.Length
			var rank *int64
			if x.Rank != nil {
				rank = x.Rank
			} else {
				rank = x.IntervalRank
			}
			assert(rank != nil && gap == *rank, "gap mismatch")
			set[fmt.Sprintf("%s|%d|%d", a.key(), *x.Length, gap)] = true
		}
		keys := make([]string, 0, len(set))
		for k := range set {
			keys = append(keys, k)
		}
		sort.Strings(keys)
		wk := w.key()
		_, dup := res[wk]
		assert(!dup, "duplicate terminal matrix")
		res[wk] = entry{length: *b.Length, rd: rd, ld: ld, bottoms: keys}
		mats[wk] = w
	}
	return &data, res, mats
}

func convolve(a, b []int64) []int64 {
	out := make([]int64, len(a)+len(b)-1)
	for i, x := range a {
		for j, y := range b {
			p := x * y
			if x != 0 && p/x != y {
				fail("overflow in convolve")
			}
			out[i+j] += p
			if out[i+j] < 0 {
				fail("overflow in convolve")
			}
		}
	}
	return out
}

func poincare(degrees []int) []int64 {
	r := []int64{1}
	for _, d := range degrees {
		ones := make([]int64, d)
		for i := range ones {
			ones[i] = 1
		}
		r = convolve(r, ones)
	}
	return r
}

func sliceEq(a, b []int64) bool {
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

func need(p *int64, name string) int64 {
	assert(p != nil, "missing key %s", name)
	return *p
}

func main() {
	matrices := flag.Int("matrices", 0, "print the generated terminal matrices of research/en_e8/eN-recursive.json (N = 6, 7 or 8) as JSON")
	flag.Parse()
	if flag.NArg() != 0 {
		fmt.Fprintln(os.Stderr, "verify-outputs: no positional arguments")
		os.Exit(2)
	}
	defer func() {
		if r := recover(); r != nil {
			fmt.Fprintf(os.Stderr, "verify-outputs: FAILED: %v\n", r)
			os.Exit(1)
		}
	}()
	if *matrices != 0 {
		n := *matrices
		if n < 4 || n > 8 {
			fmt.Fprintln(os.Stderr, "verify-outputs: -matrices N needs N in 4..8")
			os.Exit(2)
		}
		_, _, mats := checkCertificate(fmt.Sprintf("research/en_e8/e%d-recursive.json", n), n)
		list := make([]matrix, 0, len(mats))
		for _, m := range mats {
			list = append(list, m)
		}
		sort.Slice(list, func(i, j int) bool { return list[i].key() < list[j].key() })
		out, err := json.Marshal(list)
		if err != nil {
			panic(err)
		}
		fmt.Println(string(out))
		return
	}

	for _, n := range []int{6, 7} {
		_, reference, _ := checkCertificate(fmt.Sprintf("results/e%d-independent-certificate.json", n), n)
		_, flat, _ := checkCertificate(fmt.Sprintf("research/en_e8/e%d-validation.json", n), n)
		_, recursive, _ := checkCertificate(fmt.Sprintf("research/en_e8/e%d-recursive.json", n), n)
		assert(resultsEqual(reference, flat) && resultsEqual(flat, recursive), "E%d methods disagree", n)
	}

	flatData, flat, _ := checkCertificate("research/en_e8/e8-terminals.json", 8)
	recursiveData, recursive, _ := checkCertificate("research/en_e8/e8-recursive.json", 8)
	d7Data, d7, _ := checkCertificate("research/en_independent/e8-d7-terminals.json", 8)
	assert(resultsEqual(flat, recursive) && resultsEqual(recursive, d7), "E8 methods disagree")
	fr := need(flatData.RightTerminalCount, "right_terminal_count")
	rr := need(recursiveData.RightTerminalCount, "right_terminal_count")
	dr := need(d7Data.FullRightTerminals, "full_right_terminals")
	assert(fr == rr && rr == dr && dr == 2160, "right terminal counts %d %d %d", fr, rr, dr)
	fc := need(flatData.CommutingTerminals, "commuting_terminals")
	rc := need(recursiveData.CommutingTerminals, "commuting_terminals")
	dc := need(d7Data.CommutingTerminals, "commuting_terminals")
	assert(fc == rc && rc == dc && dc == 58, "commuting terminal counts %d %d %d", fc, rc, dc)
	assert(len(flat) == 6, "noncommuting terminal count %d", len(flat))
	assert(d7Data.hasParabolicHist && d7Data.hasCosetHist, "missing D7 histograms")
	assert(sliceEq(d7Data.ParabolicHistogram, poincare([]int{2, 4, 6, 7, 8, 10, 12})), "D7 parabolic histogram")
	assert(sliceEq(convolve(d7Data.ParabolicHistogram, d7Data.CosetHistogram), poincare([]int{2, 8, 12, 14, 18, 20, 24, 30})), "E8 Poincare histogram")

	gaps := make([]int64, 0, len(flat))
	for _, e := range flat {
		// next(iter(entry[3]))[2]: the gap of the (single) bottom.
		assert(len(e.bottoms) > 0, "terminal without bottoms")
		parts := strings.Split(e.bottoms[0], "|")
		g, err := strconv.ParseInt(parts[2], 10, 64)
		if err != nil {
			panic(err)
		}
		gaps = append(gaps, g)
	}
	sort.Slice(gaps, func(i, j int) bool { return gaps[i] < gaps[j] })

	summary := struct {
		Status        string  `json:"status"`
		E6E7          string  `json:"E6_E7"`
		E8            string  `json:"E8"`
		RightTerminal int     `json:"right_terminal_count"`
		Commuting     int     `json:"commuting_terminal_count"`
		Noncommuting  int     `json:"noncommuting_terminal_count"`
		EligibleGaps  []int64 `json:"eligible_gaps"`
		Histograms    string  `json:"D7_and_E8_Poincare_histograms"`
	}{
		"All assertions passed",
		"Flat and recursive methods match independent complete enumerations",
		"Flat E7, recursive parabolic, and independent D7 methods agree as integer matrices",
		2160, 58, 6, gaps,
		"verified coefficient by coefficient",
	}
	out, err := json.MarshalIndent(summary, "", "  ")
	if err != nil {
		panic(err)
	}
	fmt.Println(string(out))
}
