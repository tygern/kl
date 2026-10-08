// Command finite-descents checks the eleven finite table rows of the
// ai-review-notes: for each (type E_n, word, lower set) it verifies that the
// common left/right descent set of the word equals the stated independent set
// (computed with integer-column reflection arithmetic), that the independent
// set contains no edge of the Coxeter diagram, and that a maximum matching on
// the support edges gives alpha(support) <= |support| - |matching| with
// equality |support| - |matching| = |independent set|.
//
// Ports research/ai-review-notes/check_finite_descents.py.
// Writes research/ai-review-notes/finite-descents.json (relative to the
// current directory) and prints one summary line per row.
// Imports only the Go standard library; no other engine or internal package.
package main

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"sort"
)

type row struct {
	n     int
	word  string
	lower string
}

var rows = []row{
	{6, "1325213", "135"},
	{7, "1326213", "136"},
	{7, "13256213", "1356"},
	{7, "132543621324356", "1356"},
	{7, "1325436210321432543621324356", "1356"},
	{8, "1327213", "137"},
	{8, "13257213", "1357"},
	{8, "61327213", "1367"},
	{8, "132543721324357", "1357"},
	{8, "1325437210321432543721324357", "1357"},
	{8, "7534231270123456210321432" + "5437210321432543721324357", "1357"},
}

type result struct {
	Type               string   `json:"type"`
	Word               string   `json:"word"`
	Length             int      `json:"length"`
	Support            []int    `json:"support"`
	CommonDescents     []int    `json:"common_descents"`
	Matching           [][2]int `json:"matching"`
	IndependenceNumber int      `json:"independence_number"`
}

func fail(format string, a ...any) {
	fmt.Fprintf(os.Stderr, "finite-descents: assertion failed: "+format+"\n", a...)
	os.Exit(1)
}

func digits(s string) []int {
	out := make([]int, len(s))
	for i, c := range s {
		if c < '0' || c > '9' {
			fail("non-digit letter %q in word %q", c, s)
		}
		out[i] = int(c - '0')
	}
	return out
}

func toSet(a []int) map[int]bool {
	m := map[int]bool{}
	for _, x := range a {
		m[x] = true
	}
	return m
}

func sortedKeys(m map[int]bool) []int {
	out := []int{}
	for k := range m {
		out = append(out, k)
	}
	sort.Ints(out)
	return out
}

func multiply(n int, edges [][2]int, seq []int) [][]int64 {
	columns := make([][]int64, n)
	for j := 0; j < n; j++ {
		columns[j] = make([]int64, n)
		columns[j][j] = 1
	}
	for _, s := range seq {
		if s < 0 || s >= n {
			fail("letter %d out of range for rank %d", s, n)
		}
		old := append([]int64(nil), columns[s]...)
		for i := range old {
			columns[s][i] = -old[i]
		}
		for _, e := range edges {
			a, b := e[0], e[1]
			if s == a || s == b {
				t := b
				if s != a {
					t = a
				}
				for i := range old {
					columns[t][i] += old[i]
				}
			}
		}
	}
	return columns
}

func desc(columns [][]int64) map[int]bool {
	out := map[int]bool{}
	for j, c := range columns {
		anyNonzero, allNonneg, allNonpos, anyNeg := false, true, true, false
		for _, v := range c {
			if v != 0 {
				anyNonzero = true
			}
			if v < 0 {
				allNonneg = false
				anyNeg = true
			}
			if v > 0 {
				allNonpos = false
			}
		}
		if !anyNonzero || !(allNonneg || allNonpos) {
			fail("column %d is zero or not sign-definite", j)
		}
		if anyNeg {
			out[j] = true
		}
	}
	return out
}

func setsEqual(a, b map[int]bool) bool {
	if len(a) != len(b) {
		return false
	}
	for k := range a {
		if !b[k] {
			return false
		}
	}
	return true
}

func check(r row) result {
	n := r.n
	word := digits(r.word)
	support := toSet(word)
	edges := [][2]int{}
	for i := 0; i < n-2; i++ {
		edges = append(edges, [2]int{i, i + 1})
	}
	edges = append(edges, [2]int{2, n - 1})
	supportEdges := [][2]int{}
	for _, e := range edges {
		if support[e[0]] && support[e[1]] {
			supportEdges = append(supportEdges, e)
		}
	}
	indep := toSet(digits(r.lower))
	for _, e := range edges {
		if indep[e[0]] && indep[e[1]] {
			fail("independent set %q contains edge %v (E%d)", r.lower, e, n)
		}
	}

	rev := make([]int, len(word))
	for i, s := range word {
		rev[len(word)-1-i] = s
	}
	d1 := desc(multiply(n, edges, word))
	d2 := desc(multiply(n, edges, rev))
	if !setsEqual(d1, d2) || !setsEqual(d2, indep) {
		fail("descent sets differ for E%d word %s: left %v right %v expected %v",
			n, r.word, sortedKeys(d1), sortedKeys(d2), sortedKeys(indep))
	}

	// First maximum-size matching in increasing bitmask order (matches Python max()).
	var best [][2]int
	found := false
	for bits := 0; bits < 1<<len(supportEdges); bits++ {
		chosen := [][2]int{}
		for i, e := range supportEdges {
			if bits&(1<<i) != 0 {
				chosen = append(chosen, e)
			}
		}
		seen := map[int]bool{}
		ok := true
		for _, e := range chosen {
			for _, s := range e {
				if seen[s] {
					ok = false
				}
				seen[s] = true
			}
		}
		if ok && (!found || len(chosen) > len(best)) {
			best = chosen
			found = true
		}
	}
	if len(support)-len(best) != len(indep) {
		fail("|support|-|matching| = %d != |independent set| = %d for E%d word %s",
			len(support)-len(best), len(indep), n, r.word)
	}
	return result{
		Type:               fmt.Sprintf("E%d", n),
		Word:               r.word,
		Length:             len(word),
		Support:            sortedKeys(support),
		CommonDescents:     sortedKeys(indep),
		Matching:           best,
		IndependenceNumber: len(indep),
	}
}

func main() {
	results := make([]result, 0, len(rows))
	for _, r := range rows {
		results = append(results, check(r))
	}
	data, err := json.MarshalIndent(results, "", "  ")
	if err != nil {
		fail("marshal: %v", err)
	}
	target := filepath.Join("research", "ai-review-notes", "finite-descents.json")
	if err := os.WriteFile(target, append(data, '\n'), 0o644); err != nil {
		fmt.Fprintf(os.Stderr, "finite-descents: %v\n", err)
		os.Exit(1)
	}
	for _, r := range results {
		fmt.Println(r.Type, r.Length, ints(r.Support), ints(r.CommonDescents), pairs(r.Matching))
	}
}

func ints(a []int) string {
	s := "["
	for i, x := range a {
		if i > 0 {
			s += ", "
		}
		s += fmt.Sprint(x)
	}
	return s + "]"
}

func pairs(a [][2]int) string {
	s := "["
	for i, e := range a {
		if i > 0 {
			s += ", "
		}
		s += fmt.Sprintf("(%d, %d)", e[0], e[1])
	}
	return s + "]"
}
