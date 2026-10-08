package main

// E6 integer-matrix terminal certificate: port of research/verify_e6_independent.py.

import (
	"crypto/sha256"
	_ "embed"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"math/bits"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

//go:embed e6.go
var e6Source []byte

const e6N = 6

var e6Edges = [][2]int{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {2, 5}}

// e6Mat holds the columns of an integral matrix: column c occupies [c*6, c*6+6).
type e6Mat [e6N * e6N]int64

type lengthDist []int64

// MarshalJSON writes the distribution as an object keyed by length, in numeric order.
func (d lengthDist) MarshalJSON() ([]byte, error) {
	var sb strings.Builder
	sb.WriteByte('{')
	for k, n := range d {
		if k > 0 {
			sb.WriteByte(',')
		}
		sb.WriteString(strconv.Quote(strconv.Itoa(k)))
		sb.WriteByte(':')
		sb.WriteString(strconv.FormatInt(n, 10))
	}
	sb.WriteByte('}')
	return []byte(sb.String()), nil
}

type e6Candidate struct {
	Word         []int `json:"word"`
	Length       int   `json:"length"`
	IntervalRank int   `json:"interval_rank"`
}

type e6Terminal struct {
	Word          []int         `json:"word"`
	Length        int           `json:"length"`
	RightDescents []int         `json:"right_descents"`
	LeftDescents  []int         `json:"left_descents"`
	Candidates    []e6Candidate `json:"candidates"`
}

type e6Certificate struct {
	Type               string       `json:"type"`
	Edges              [][2]int     `json:"edges"`
	ElementCount       int          `json:"element_count"`
	FcCount            int          `json:"fc_count"`
	LengthDistribution lengthDist   `json:"length_distribution"`
	FcStarChecks       int          `json:"fc_star_checks"`
	CommutingTerminals int          `json:"commuting_terminals"`
	NonFcWeakBad       []e6Terminal `json:"non_fc_weak_bad_terminals"`
	Status             string       `json:"status"`
	CodeSHA256         string       `json:"code_sha256"`
}

func e6Adjacency() [e6N][]int {
	var adj [e6N][]int
	for i := 0; i < e6N; i++ {
		for _, e := range e6Edges {
			if e[0] == i {
				adj[i] = append(adj[i], e[1])
			} else if e[1] == i {
				adj[i] = append(adj[i], e[0])
			}
		}
	}
	return adj
}

// e6Right is the right action of the simple reflection s on the column-stored matrix a.
func e6Right(a *e6Mat, s int, adj *[e6N][]int) e6Mat {
	b := *a
	for k := 0; k < e6N; k++ {
		b[s*e6N+k] = -a[s*e6N+k]
	}
	for _, t := range adj[s] {
		for k := 0; k < e6N; k++ {
			b[t*e6N+k] = a[t*e6N+k] + a[s*e6N+k]
		}
	}
	return b
}

// e6Negative reports whether the column is nonpositive and asserts it has a single sign.
func e6Negative(a *e6Mat, c int) bool {
	allNonneg, allNonpos, anyNeg := true, true, false
	for k := 0; k < e6N; k++ {
		v := a[c*e6N+k]
		if v < 0 {
			allNonneg = false
			anyNeg = true
		}
		if v > 0 {
			allNonpos = false
		}
	}
	check(allNonneg || allNonpos, "column %d has mixed signs", c)
	return anyNeg
}

func hasBit(mask uint, s int) bool { return mask>>uint(s)&1 == 1 }

func maskList(mask uint) []int {
	out := make([]int, 0, e6N)
	for s := 0; s < e6N; s++ {
		if hasBit(mask, s) {
			out = append(out, s)
		}
	}
	return out
}

func e6Compute() e6Certificate {
	adj := e6Adjacency()
	var identity e6Mat
	for i := 0; i < e6N; i++ {
		identity[i*e6N+i] = 1
	}
	elements := []e6Mat{identity}
	index := map[e6Mat]int{identity: 0}
	parent := []int{-1}
	last := []int{-1}
	length := []int{0}
	var action [][e6N]int
	for k := 0; k < len(elements); k++ {
		a := elements[k]
		var row [e6N]int
		for s := 0; s < e6N; s++ {
			b := e6Right(&a, s, &adj)
			id, ok := index[b]
			if !ok {
				id = len(elements)
				index[b] = id
				elements = append(elements, b)
				parent = append(parent, k)
				last = append(last, s)
				length = append(length, length[k]+1)
			}
			row[s] = id
		}
		action = append(action, row)
	}
	n := len(elements)
	check(n == 51840, "E6 element count %d != 51840", n)

	word := func(k int) []int {
		w := make([]int, 0, length[k])
		for v := k; v != 0; v = parent[v] {
			w = append(w, last[v])
		}
		for i, j := 0, len(w)-1; i < j; i, j = i+1, j-1 {
			w[i], w[j] = w[j], w[i]
		}
		return w
	}

	inv := make([]int, n)
	for k := 0; k < n; k++ {
		z := 0
		for v := k; v != 0; v = parent[v] {
			z = action[z][last[v]]
		}
		inv[k] = z
	}
	desc := make([]uint, n)
	for k := 0; k < n; k++ {
		for s := 0; s < e6N; s++ {
			if e6Negative(&elements[k], s) {
				desc[k] |= 1 << uint(s)
			}
		}
	}
	leftDesc := make([]uint, n)
	left := make([][e6N]int, n)
	for k := 0; k < n; k++ {
		leftDesc[k] = desc[inv[k]]
		for s := 0; s < e6N; s++ {
			left[k][s] = inv[action[inv[k]][s]]
		}
	}
	fc := make([]bool, n)
	fcCount := 0
	for k := 0; k < n; k++ {
		// A reduced braid suffix is detected after stripping descents. Matsumoto
		// gives the equivalence with full commutativity in simply laced type.
		ok := true
		for s := 0; s < e6N; s++ {
			if hasBit(desc[k], s) {
				t := action[k][s]
				check(t < k, "descent product not earlier in order")
				if !fc[t] {
					ok = false
				}
			}
		}
		if ok {
			for _, e := range e6Edges {
				if hasBit(desc[k], e[0]) && hasBit(desc[k], e[1]) {
					ok = false
				}
			}
		}
		fc[k] = ok
		if ok {
			fcCount++
		}
	}
	check(fcCount == 662, "E6 fc count %d != 662", fcCount)

	leq := func(x, w int) bool {
		for x != w {
			if length[x] >= length[w] {
				return false
			}
			s := bits.TrailingZeros(desc[w])
			if hasBit(desc[x], s) {
				x = action[x][s]
			}
			w = action[w][s]
		}
		return true
	}

	var terminal []e6Terminal
	starChecks := 0
	commuting := 0
	for k := 0; k < n; k++ {
		weakBad := true
		for side := 0; side < 2; side++ {
			table := func(k, s int) int {
				if side == 0 {
					return action[k][s]
				}
				return left[k][s]
			}
			ds := func(k int) uint {
				if side == 0 {
					return desc[k]
				}
				return leftDesc[k]
			}
			for s := 0; s < e6N; s++ {
				if !hasBit(ds(k), s) {
					continue
				}
				for _, t := range adj[s] {
					if hasBit(ds(table(k, s)), t) {
						weakBad = false
					}
				}
			}
			if fc[k] {
				// Check the complete length-three star domain, including stars
				// that raise length. There is exactly one domain-preserving move.
				for _, e := range e6Edges {
					s, t := e[0], e[1]
					pair := uint(1)<<uint(s) | uint(1)<<uint(t)
					if bits.OnesCount(ds(k)&pair) == 1 {
						var moves []int
						for _, u := range []int{s, t} {
							v := table(k, u)
							if bits.OnesCount(ds(v)&pair) == 1 {
								moves = append(moves, v)
							}
						}
						check(len(moves) == 1, "star at element %d has %d moves", k, len(moves))
						check(fc[moves[0]], "star move leaves FC set")
						starChecks++
					}
				}
			}
		}
		if !weakBad {
			continue
		}
		if fc[k] {
			commuting++
			w := word(k)
			seen := uint(0)
			for _, s := range w {
				check(!hasBit(seen, s), "commuting terminal repeats a letter")
				seen |= 1 << uint(s)
			}
			for _, e := range e6Edges {
				check(!(hasBit(seen, e[0]) && hasBit(seen, e[1])), "commuting terminal contains an edge")
			}
		} else {
			cands := make([]e6Candidate, 0)
			for x := 0; x < n; x++ {
				if fc[x] && desc[k]&^desc[x] == 0 && leftDesc[k]&^leftDesc[x] == 0 && leq(x, k) {
					cands = append(cands, e6Candidate{Word: word(x), Length: length[x], IntervalRank: length[k] - length[x]})
				}
			}
			terminal = append(terminal, e6Terminal{
				Word: word(k), Length: length[k],
				RightDescents: maskList(desc[k]), LeftDescents: maskList(leftDesc[k]),
				Candidates: cands,
			})
		}
	}
	check(len(terminal) == 1, "E6 non-FC weak-bad terminal count %d != 1", len(terminal))
	for _, w := range terminal {
		for _, c := range w.Candidates {
			check(c.IntervalRank%2 == 0, "odd interval rank %d", c.IntervalRank)
		}
	}

	maxLen := 0
	for _, l := range length {
		if l > maxLen {
			maxLen = l
		}
	}
	dist := make(lengthDist, maxLen+1)
	for _, l := range length {
		dist[l]++
	}
	sum := sha256.Sum256(e6Source)
	if terminal == nil {
		terminal = []e6Terminal{}
	}
	return e6Certificate{
		Type: "E6", Edges: e6Edges, ElementCount: n, FcCount: fcCount,
		LengthDistribution: dist, FcStarChecks: starChecks, CommutingTerminals: commuting,
		NonFcWeakBad: terminal,
		Status:       "Independent finite terminal enumeration; all compatible bad ranks even.",
		CodeSHA256:   hex.EncodeToString(sum[:]),
	}
}

// runE6 computes the certificate, writes results/e6-independent-certificate.json and prints it.
func runE6() {
	cert := e6Compute()
	data, err := json.MarshalIndent(cert, "", "  ")
	if err != nil {
		fmt.Fprintln(os.Stderr, "matrix-search:", err)
		os.Exit(1)
	}
	data = append(data, '\n')
	out := filepath.Join("results", "e6-independent-certificate.json")
	if err := os.MkdirAll("results", 0o755); err != nil {
		fmt.Fprintln(os.Stderr, "matrix-search:", err)
		os.Exit(1)
	}
	if err := os.WriteFile(out, data, 0o644); err != nil {
		fmt.Fprintln(os.Stderr, "matrix-search:", err)
		os.Exit(1)
	}
	os.Stdout.Write(data)
}
