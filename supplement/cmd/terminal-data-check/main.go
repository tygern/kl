// Command terminal-data-check is a referee check of the terminal data of the
// exceptional groups that is independent of both exceptional-group engines.
// It checks the full Poincare polynomials of E6 and E7 against the length
// distributions of results/e{6,7}-independent-certificate.json, the commuting
// terminal counts (independent sets of the diagram), the root-system size and
// maximum absolute root coefficient, the transcription of the four E7 bad
// terminal words and their bottoms, and the identification of the D6 pair
// with the signed-permutation model.
//
// Ports research/exceptional_referee/verify_terminal_data.py.
// Reads results/e6-independent-certificate.json and
// results/e7-independent-certificate.json relative to the current directory
// (the repository root or the payload copy) and prints the JSON summary to
// stdout; the runner captures it to research/exceptional_referee/checks.json.
// Exits non-zero on any failed assertion.
// Imports only the Go standard library; no other engine or internal package.
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

// poincare returns the coefficients of prod_d (1 + q + ... + q^(d-1)).
func poincare(degrees []int) []int64 {
	result := []int64{1}
	for _, d := range degrees {
		next := make([]int64, len(result)+d-1)
		for i, c := range result {
			for j := 0; j < d; j++ {
				next[i+j] += c
			}
		}
		result = next
	}
	return result
}

// independentSets counts the subsets of {0..n-1} containing no edge.
func independentSets(n int, edges [][2]int) int {
	count := 0
	for mask := 0; mask < 1<<n; mask++ {
		ok := true
		for _, e := range edges {
			if mask&(1<<e[0]) != 0 && mask&(1<<e[1]) != 0 {
				ok = false
				break
			}
		}
		if ok {
			count++
		}
	}
	return count
}

// rootBounds returns the number of roots of the simply laced diagram (as the
// orbit of the simple roots under the simple reflections) and the maximum
// absolute coefficient; every root must be nonnegative or nonpositive.
func rootBounds(n int, edges [][2]int) (int, int, error) {
	adjacency := make([][]int, n)
	for _, e := range edges {
		adjacency[e[0]] = append(adjacency[e[0]], e[1])
		adjacency[e[1]] = append(adjacency[e[1]], e[0])
	}
	key := func(r []int) string {
		parts := make([]string, len(r))
		for i, c := range r {
			parts[i] = strconv.Itoa(c)
		}
		return strings.Join(parts, ",")
	}
	roots := map[string][]int{}
	var todo [][]int
	for j := 0; j < n; j++ {
		r := make([]int, n)
		r[j] = 1
		roots[key(r)] = r
		todo = append(todo, r)
	}
	for idx := 0; idx < len(todo); idx++ {
		root := todo[idx]
		for s := 0; s < n; s++ {
			image := append([]int(nil), root...)
			sum := 0
			for _, t := range adjacency[s] {
				sum += root[t]
			}
			image[s] = sum - root[s]
			k := key(image)
			if _, ok := roots[k]; !ok {
				roots[k] = image
				todo = append(todo, image)
			}
		}
	}
	maxAbs := 0
	for _, r := range roots {
		nonneg, nonpos := true, true
		for _, c := range r {
			if c < 0 {
				nonneg = false
			}
			if c > 0 {
				nonpos = false
			}
			a := c
			if a < 0 {
				a = -a
			}
			if a > maxAbs {
				maxAbs = a
			}
		}
		if !nonneg && !nonpos {
			return 0, 0, fmt.Errorf("root %v has mixed signs", r)
		}
	}
	return len(roots), maxAbs, nil
}

// d6Product multiplies the generators of the word on the right, in Gern's
// one-based labelling of D6 as signed permutations.
func d6Product(word []int) ([6]int, error) {
	w := [6]int{1, 2, 3, 4, 5, 6}
	for _, s := range word {
		switch {
		case s == 1:
			w[0], w[1] = -w[1], -w[0]
		case s >= 2 && s <= 6:
			w[s-2], w[s-1] = w[s-1], w[s-2]
		default:
			return w, fmt.Errorf("D6 generator label %d out of range", s)
		}
	}
	return w, nil
}

func d6Length(w [6]int) int {
	n := 0
	for i := 0; i < 6; i++ {
		for j := i + 1; j < 6; j++ {
			if w[i] > w[j] {
				n++
			}
			if -w[i] > w[j] {
				n++
			}
		}
	}
	return n
}

type bottom struct {
	Word []int `json:"word"`
	Rank int   `json:"rank"`
}

type badTerminal struct {
	Word    []int    `json:"word"`
	Bottoms []bottom `json:"bottoms"`
}

type cert struct {
	LengthDistribution json.RawMessage `json:"length_distribution"`
	Bad                []badTerminal   `json:"bad"`
}

func readCert(root string, n int) (*cert, error) {
	path := filepath.Join(root, "results", fmt.Sprintf("e%d-independent-certificate.json", n))
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var c cert
	if err := json.Unmarshal(data, &c); err != nil {
		return nil, fmt.Errorf("%s: %v", path, err)
	}
	return &c, nil
}

// histogram accepts either a list or an object keyed by the decimal index.
func histogram(raw json.RawMessage) ([]int64, error) {
	var list []int64
	if err := json.Unmarshal(raw, &list); err == nil {
		return list, nil
	}
	var m map[string]int64
	if err := json.Unmarshal(raw, &m); err != nil {
		return nil, fmt.Errorf("length_distribution is neither list nor object: %v", err)
	}
	out := make([]int64, len(m))
	for i := range out {
		v, ok := m[strconv.Itoa(i)]
		if !ok {
			return nil, fmt.Errorf("length_distribution lacks key %d", i)
		}
		out[i] = v
	}
	return out, nil
}

func equalInts64(a, b []int64) bool {
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

func equalInts(a, b []int) bool {
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

func wordString(w []int) string {
	var sb strings.Builder
	for _, s := range w {
		sb.WriteString(strconv.Itoa(s))
	}
	return sb.String()
}

func mapWord(w []int, eToD map[int]int) ([]int, error) {
	out := make([]int, len(w))
	for i, s := range w {
		d, ok := eToD[s]
		if !ok {
			return nil, fmt.Errorf("E7 label %d has no D6 image", s)
		}
		out[i] = d
	}
	return out, nil
}

// output is the JSON summary in the key order of the original.
type output struct {
	Status                   string         `json:"status"`
	PoincarePolynomials      string         `json:"poincare_polynomials"`
	CommutingTerminalCounts  []int          `json:"commuting_terminal_counts"`
	MaximumAbsoluteRootCoeff []int          `json:"maximum_absolute_root_coefficients"`
	E7TerminalRanks          []int          `json:"E7_terminal_ranks"`
	D6MapGernToE7            map[string]int `json:"D6_map_Gern_to_E7"`
	D6Bottom                 []int          `json:"D6_bottom"`
	D6Top                    []int          `json:"D6_top"`
}

func run(root string) (*output, error) {
	certs := map[int]*cert{}
	for _, n := range []int{6, 7} {
		c, err := readCert(root, n)
		if err != nil {
			return nil, err
		}
		certs[n] = c
	}
	degrees := map[int][]int{6: {2, 5, 6, 8, 9, 12}, 7: {2, 6, 8, 10, 12, 14, 18}}
	wantIndep := map[int]int{6: 22, 7: 36}
	wantRoots := map[int][2]int{6: {72, 3}, 7: {126, 4}}
	for _, n := range []int{6, 7} {
		hist, err := histogram(certs[n].LengthDistribution)
		if err != nil {
			return nil, err
		}
		if !equalInts64(hist, poincare(degrees[n])) {
			return nil, fmt.Errorf("E%d length distribution differs from the Poincare polynomial", n)
		}
		var edges [][2]int
		for i := 0; i < n-2; i++ {
			edges = append(edges, [2]int{i, i + 1})
		}
		edges = append(edges, [2]int{2, n - 1})
		if got := independentSets(n, edges); got != wantIndep[n] {
			return nil, fmt.Errorf("E%d independent sets %d, want %d", n, got, wantIndep[n])
		}
		count, maxAbs, err := rootBounds(n, edges)
		if err != nil {
			return nil, err
		}
		if [2]int{count, maxAbs} != wantRoots[n] {
			return nil, fmt.Errorf("E%d root bounds (%d, %d), want %v", n, count, maxAbs, wantRoots[n])
		}
	}

	bad := certs[7].Bad
	expectedWords := []string{"1326213", "13256213", "132543621324356", "1325436210321432543621324356"}
	if len(bad) != len(expectedWords) {
		return nil, fmt.Errorf("E7 has %d bad terminals, want %d", len(bad), len(expectedWords))
	}
	wantLengths := []int{7, 8, 15, 28}
	wantRanks := []int{4, 4, 11, 24}
	var ranks []int
	for i, b := range bad {
		if wordString(b.Word) != expectedWords[i] {
			return nil, fmt.Errorf("E7 bad word %d is %s, want %s", i, wordString(b.Word), expectedWords[i])
		}
		if len(b.Word) != wantLengths[i] {
			return nil, fmt.Errorf("E7 bad word %d has length %d, want %d", i, len(b.Word), wantLengths[i])
		}
		if len(b.Bottoms) != 1 {
			return nil, fmt.Errorf("E7 bad word %d has %d bottoms, want 1", i, len(b.Bottoms))
		}
		ranks = append(ranks, b.Bottoms[0].Rank)
	}
	if !equalInts(ranks, wantRanks) {
		return nil, fmt.Errorf("E7 terminal ranks %v, want %v", ranks, wantRanks)
	}
	support := map[int]bool{}
	for _, s := range bad[3].Word {
		support[s] = true
	}
	if len(support) != 7 {
		return nil, fmt.Errorf("the longest E7 bad word has support of size %d, want 7", len(support))
	}
	for s := 0; s < 7; s++ {
		if !support[s] {
			return nil, fmt.Errorf("the longest E7 bad word misses generator %d", s)
		}
	}

	// In Gern, the D6 edges are 1--3, 2--3, 3--4--5--6.
	dToE := map[int]int{1: 1, 2: 6, 3: 2, 4: 3, 5: 4, 6: 5}
	eToD := map[int]int{}
	for d, e := range dToE {
		eToD[e] = d
	}
	dEdges := [][2]int{{1, 3}, {2, 3}, {3, 4}, {4, 5}, {5, 6}}
	var mapped []string
	for _, e := range dEdges {
		a, b := dToE[e[0]], dToE[e[1]]
		if a > b {
			a, b = b, a
		}
		mapped = append(mapped, fmt.Sprintf("%d-%d", a, b))
	}
	sort.Strings(mapped)
	wantMapped := []string{"1-2", "2-3", "2-6", "3-4", "4-5"}
	if strings.Join(mapped, " ") != strings.Join(wantMapped, " ") {
		return nil, fmt.Errorf("mapped D6 edges %v, want %v", mapped, wantMapped)
	}
	b2 := bad[2]
	topWord, err := mapWord(b2.Word, eToD)
	if err != nil {
		return nil, err
	}
	bottomWord, err := mapWord(b2.Bottoms[0].Word, eToD)
	if err != nil {
		return nil, err
	}
	top, err := d6Product(topWord)
	if err != nil {
		return nil, err
	}
	bot, err := d6Product(bottomWord)
	if err != nil {
		return nil, err
	}
	if top != [6]int{-1, -6, 3, -4, 5, -2} {
		return nil, fmt.Errorf("D6 top %v unexpected", top)
	}
	if bot != [6]int{-1, -2, 4, 3, 6, 5} {
		return nil, fmt.Errorf("D6 bottom %v unexpected", bot)
	}
	if d6Length(top) != 15 || d6Length(bot) != 4 {
		return nil, fmt.Errorf("D6 lengths (%d, %d), want (15, 4)", d6Length(top), d6Length(bot))
	}
	mapOut := map[string]int{}
	for d, e := range dToE {
		mapOut[strconv.Itoa(d)] = e
	}
	return &output{
		Status:                   "All assertions passed",
		PoincarePolynomials:      "E6 and E7 match coefficient by coefficient",
		CommutingTerminalCounts:  []int{22, 36},
		MaximumAbsoluteRootCoeff: []int{3, 4},
		E7TerminalRanks:          ranks,
		D6MapGernToE7:            mapOut,
		D6Bottom:                 bot[:],
		D6Top:                    top[:],
	}, nil
}

func main() {
	out, err := run(".")
	if err != nil {
		fmt.Fprintf(os.Stderr, "terminal-data-check: %v\n", err)
		os.Exit(1)
	}
	data, err := json.MarshalIndent(out, "", "  ")
	if err != nil {
		fmt.Fprintf(os.Stderr, "terminal-data-check: %v\n", err)
		os.Exit(1)
	}
	fmt.Println(string(data))
}
