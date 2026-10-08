package main

import (
	"fmt"
	"math/bits"
	"os"
	"sort"
)

// fatal prints a message and exits non-zero (the original called exit(1)).
func fatal(msg string) {
	fmt.Fprintln(os.Stderr, msg)
	os.Exit(1)
}

// Model is a Bruhat lower ideal of a Coxeter group of rank <= 8 together with
// its multiplication tables and descent masks.
type Model struct {
	n, N    int
	rmul    [][8]int32 // -1 if outside the ideal
	lmul    [][8]int32
	length  []uint8
	rdes    []uint8 // descent masks: bit g = Gern's s_{g+1}
	ldes    []uint8
	coatoms [][]int32 // z < y with l(z) = l(y) - 1
	idE     int
	idTop   int
	idX     int
	canon   []string // canonical reduced word (smallest left descent first)
}

func newModel(n, count int) *Model {
	m := &Model{n: n, N: count}
	m.rmul = make([][8]int32, count)
	m.lmul = make([][8]int32, count)
	m.length = make([]uint8, count)
	m.rdes = make([]uint8, count)
	m.ldes = make([]uint8, count)
	m.coatoms = make([][]int32, count)
	return m
}

func canonicalWords(m *Model) {
	order := make([]int, m.N)
	for i := range order {
		order[i] = i
	}
	sort.SliceStable(order, func(a, b int) bool { return m.length[order[a]] < m.length[order[b]] })
	m.canon = make([]string, m.N)
	for _, id := range order {
		if m.length[id] == 0 {
			continue
		}
		g := bits.TrailingZeros8(m.ldes[id])
		u := m.lmul[id][g]
		if u < 0 || m.length[u] != m.length[id]-1 {
			fatal("canonical word error")
		}
		m.canon[id] = string(rune('1'+g)) + m.canon[u]
	}
}
