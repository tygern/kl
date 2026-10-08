package main

import (
	"fmt"
	"sort"
)

// maxN bounds the rank of the affine model (arrays are fixed-size so that
// elements are comparable and usable as map keys).
const maxN = 8

// affElem is a faithful affine signed-permutation action: the one-line signed
// permutation p and the translation vector t (only the first n entries used).
type affElem struct {
	p [maxN]int
	t [maxN]int
}

func abs(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

// affineD is the ball of radius L in the affine Weyl group of type D_n^(1)
// (rank n+1, generators 0..n) in the Python class AffineD of
// affine_d_crowns.py. Simple 0 flips/swaps the first pair; i=1..n-1 swaps
// positions i-1,i; the affine generator n flips/swaps the last pair and adds
// e_(n-1)+e_n to the translation part.
type affineD struct {
	n       int
	e       affElem
	els     []affElem
	ids     map[affElem]int
	words   [][]int
	lengths []int
	down    [][]int
}

func newAffineD(n, L int) *affineD {
	if n < 3 || n > maxN {
		fail("unsupported rank %d", n)
	}
	g := &affineD{n: n, ids: map[affElem]int{}}
	for i := 0; i < n; i++ {
		g.e.p[i] = i + 1
	}
	g.els = []affElem{g.e}
	g.ids[g.e] = 0
	g.words = [][]int{{}}
	g.lengths = []int{0}
	for wi := 0; wi < len(g.els); wi++ { // the slice grows during the loop, as in the original
		if g.lengths[wi] == L {
			continue
		}
		for s := 0; s <= n; s++ {
			z := g.simple(g.els[wi], s)
			if _, ok := g.ids[z]; !ok {
				g.ids[z] = len(g.els)
				g.els = append(g.els, z)
				word := append(append([]int{}, g.words[wi]...), s)
				g.words = append(g.words, word)
				g.lengths = append(g.lengths, g.lengths[wi]+1)
			}
		}
	}
	g.down = make([][]int, len(g.words))
	for wi, Q := range g.words {
		seen := map[int]bool{}
		var row []int
		for j := range Q {
			sub := make([]int, 0, len(Q)-1)
			sub = append(sub, Q[:j]...)
			sub = append(sub, Q[j+1:]...)
			id, ok := g.ids[g.word(sub)]
			if !ok {
				fail("subword element outside the ball")
			}
			if g.lengths[id] == len(Q)-1 && !seen[id] {
				seen[id] = true
				row = append(row, id)
			}
		}
		sort.Ints(row)
		g.down[wi] = row
	}
	return g
}

// simple returns w * s_s (right multiplication by the simple reflection s).
func (g *affineD) simple(w affElem, s int) affElem {
	n := g.n
	switch {
	case s == 0:
		w.p[0], w.p[1] = -w.p[1], -w.p[0]
	case s < n:
		w.p[s-1], w.p[s] = w.p[s], w.p[s-1]
	default:
		for _, a := range []int{w.p[n-2], w.p[n-1]} {
			if a > 0 {
				w.t[abs(a)-1]++
			} else {
				w.t[abs(a)-1]--
			}
		}
		w.p[n-2], w.p[n-1] = -w.p[n-1], -w.p[n-2]
	}
	return w
}

func (g *affineD) word(Q []int) affElem {
	w := g.e
	for _, s := range Q {
		w = g.simple(w, s)
	}
	return w
}

type unfoldingResult struct {
	TopWord               []int    `json:"top_word"`
	Length                int      `json:"length"`
	IntervalSize          int      `json:"interval_size"`
	RankVector            []int    `json:"rank_vector"`
	DetectedRank3Edges    [][2]int `json:"detected_rank3_edges"`
	MaximumDetectedDegree int      `json:"maximum_detected_degree"`
}

// runUnfolding certifies the degree-four obstruction (check_unfolding.py).
func runUnfolding() unfoldingResult {
	W := newAffineD(4, 10)
	Q := []int{0, 1, 3, 4, 2, 0, 1, 3, 4, 2}
	w, ok := W.ids[W.word(Q)]
	if !ok {
		fail("top element not in the ball")
	}
	if W.lengths[w] != 10 {
		fail("length of top element is %d, expected 10", W.lengths[w])
	}
	nodes := map[int]bool{w: true}
	stack := []int{w}
	for len(stack) > 0 {
		z := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		for _, a := range W.down[z] {
			if !nodes[a] {
				nodes[a] = true
				stack = append(stack, a)
			}
		}
	}
	ranks := make([]int, 11)
	edgeSet := map[[2]int]bool{}
	for z := range nodes {
		ranks[W.lengths[z]]++
		if W.lengths[z] == 3 {
			support := map[int]bool{}
			for _, s := range W.words[z] {
				support[s] = true
			}
			if len(support) == 2 {
				var pair []int
				for s := range support {
					pair = append(pair, s)
				}
				sort.Ints(pair)
				edgeSet[[2]int{pair[0], pair[1]}] = true
			}
		}
	}
	var edges [][2]int
	for e := range edgeSet {
		edges = append(edges, e)
	}
	sort.Slice(edges, func(i, j int) bool {
		if edges[i][0] != edges[j][0] {
			return edges[i][0] < edges[j][0]
		}
		return edges[i][1] < edges[j][1]
	})
	want := [][2]int{{0, 2}, {1, 2}, {2, 3}, {2, 4}}
	if fmt.Sprint(edges) != fmt.Sprint(want) {
		fail("detected rank-3 edges %v, expected %v", edges, want)
	}
	// Coxeter relations in the faithful affine model.
	for s := 0; s < 5; s++ {
		if W.word([]int{s, s}) != W.e {
			fail("s_%d^2 is not the identity", s)
		}
		for t := 0; t < s; t++ {
			m := 2
			if s == 2 || t == 2 {
				m = 3
			}
			var rel []int
			for k := 0; k < m; k++ {
				rel = append(rel, s, t)
			}
			if W.word(rel) != W.e {
				fail("(s_%d s_%d)^%d is not the identity", s, t, m)
			}
		}
	}
	return unfoldingResult{
		TopWord:               Q,
		Length:                10,
		IntervalSize:          len(nodes),
		RankVector:            ranks,
		DetectedRank3Edges:    edges,
		MaximumDetectedDegree: 4,
	}
}
