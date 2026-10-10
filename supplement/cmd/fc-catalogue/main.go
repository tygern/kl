// Command fc-catalogue enumerates all fully commutative (FC) elements of the
// generalized Coxeter group E_n (rank n; Dynkin diagram: a
// path 0-1-...-(n-2) with an extra node n-1 attached to node 2) by direct
// breadth-first search over exact integer root matrices, retaining only FC
// states. For each element it records a reduced word, length, and the right and
// left descent masks.
//
// Ports: research/en_independent/fc_catalogue.cpp
//
// Usage: fc-catalogue -rank N   (N = 6..10; the supplement uses 6..9)
// Output: one JSON document on stdout (the runner redirects it to eN-fc.json).
//
// Imports: standard library only. This command is self-contained
// and imports no internal package and none of the other engines.
package main

import (
	"bufio"
	"encoding/binary"
	"flag"
	"fmt"
	"os"
	"strconv"
)

type mat []int32

type engine struct {
	n       int
	adj     [][]int
	els     []mat
	parents []uint32
	lengths []uint32
	last    []uint8
	ids     map[string]uint32
}

func newEngine(n int) *engine {
	e := &engine{n: n, adj: make([][]int, n), ids: map[string]uint32{}}
	for s := 0; s < n-2; s++ {
		e.adj[s] = append(e.adj[s], s+1)
		e.adj[s+1] = append(e.adj[s+1], s)
	}
	e.adj[2] = append(e.adj[2], n-1)
	e.adj[n-1] = append(e.adj[n-1], 2)
	return e
}

func key(a mat) string {
	b := make([]byte, 4*len(a))
	for i, v := range a {
		binary.LittleEndian.PutUint32(b[4*i:], uint32(v))
	}
	return string(b)
}

// right applies the right multiplication by the simple reflection s on the
// root matrix a, with an overflow guard on every sum (int32 range, as in the
// original's assertion).
func (e *engine) right(a mat, s int) mat {
	n := e.n
	b := make(mat, len(a))
	copy(b, a)
	for k := 0; k < n; k++ {
		b[s*n+k] = -a[s*n+k]
	}
	for _, t := range e.adj[s] {
		for k := 0; k < n; k++ {
			c := int64(a[t*n+k]) + int64(a[s*n+k])
			if c < -2147483648 || c > 2147483647 {
				panic("fc-catalogue: int32 overflow in root matrix")
			}
			b[t*n+k] = int32(c)
		}
	}
	return b
}

// desc returns the descent mask; every row must be wholly positive or wholly
// negative (assertion of the original).
func (e *engine) desc(a mat) uint32 {
	n := e.n
	var d uint32
	for s := 0; s < n; s++ {
		p, m := false, false
		for k := 0; k < n; k++ {
			if a[s*n+k] > 0 {
				p = true
			}
			if a[s*n+k] < 0 {
				m = true
			}
		}
		if p == m {
			panic("fc-catalogue: assertion p != m failed (row neither positive nor negative)")
		}
		if m {
			d |= 1 << uint(s)
		}
	}
	return d
}

func (e *engine) independent(d uint32) bool {
	for s := 0; s < e.n; s++ {
		if d&(1<<uint(s)) != 0 {
			for _, t := range e.adj[s] {
				if d&(1<<uint(t)) != 0 {
					return false
				}
			}
		}
	}
	return true
}

func (e *engine) word(i uint32) []int {
	var w []int
	for i != 0 {
		w = append(w, int(e.last[i]))
		i = e.parents[i]
	}
	for a, b := 0, len(w)-1; a < b; a, b = a+1, b-1 {
		w[a], w[b] = w[b], w[a]
	}
	return w
}

func (e *engine) run(out *bufio.Writer) {
	n := e.n
	id := make(mat, n*n)
	for s := 0; s < n; s++ {
		id[s*n+s] = 1
	}
	e.els = append(e.els, id)
	e.ids[key(id)] = 0
	e.parents = append(e.parents, 0)
	e.lengths = append(e.lengths, 0)
	e.last = append(e.last, 0)
	var attempted uint64
	for i := uint32(0); int(i) < len(e.els); i++ {
		a := e.els[i]
		d := e.desc(a)
		for s := 0; s < n; s++ {
			if d&(1<<uint(s)) != 0 {
				continue
			}
			attempted++
			b := e.right(a, s)
			kb := key(b)
			if _, ok := e.ids[kb]; ok {
				continue
			}
			bd := e.desc(b)
			if !e.independent(bd) {
				continue
			}
			good := true
			for t := 0; t < n; t++ {
				if bd&(1<<uint(t)) != 0 {
					j, ok := e.ids[key(e.right(b, t))]
					if !ok || e.lengths[j] != e.lengths[i] {
						good = false
						break
					}
				}
			}
			if !good {
				continue
			}
			j := uint32(len(e.els))
			e.ids[kb] = j
			e.els = append(e.els, b)
			e.parents = append(e.parents, i)
			e.last = append(e.last, uint8(s))
			e.lengths = append(e.lengths, e.lengths[i]+1)
		}
	}
	maxLen := e.lengths[len(e.lengths)-1]
	hist := make([]uint64, maxLen+1)
	maxcoord := 0
	for i := range e.els {
		hist[e.lengths[i]]++
		if !e.independent(e.desc(e.els[i])) {
			panic("fc-catalogue: assertion failed: descent set of a retained element is not independent")
		}
		for _, c := range e.els[i] {
			v := int(c)
			if v < 0 {
				v = -v
			}
			if v > maxcoord {
				maxcoord = v
			}
		}
	}
	fmt.Fprintf(out, "{\"rank\":%d,\"fc_count\":%d,\"max_fc_length\":%d,\"max_root_coordinate\":%d,\"ascents_tested\":%d,\"length_distribution\":[",
		n, len(e.els), maxLen, maxcoord, attempted)
	for k, h := range hist {
		if k > 0 {
			out.WriteString(",")
		}
		out.WriteString(strconv.FormatUint(h, 10))
	}
	out.WriteString("],\"elements\":[\n")
	for i := range e.els {
		w := e.word(uint32(i))
		inverse := make(mat, len(id))
		copy(inverse, id)
		for k := len(w) - 1; k >= 0; k-- {
			inverse = e.right(inverse, w[k])
		}
		if i > 0 {
			out.WriteString(",\n")
		}
		out.WriteString("{\"word\":[")
		for k, x := range w {
			if k > 0 {
				out.WriteString(",")
			}
			out.WriteString(strconv.Itoa(x))
		}
		fmt.Fprintf(out, "],\"length\":%d,\"Rmask\":%d,\"Lmask\":%d}", e.lengths[i], e.desc(e.els[i]), e.desc(inverse))
	}
	out.WriteString("\n],\"complete\":true,\"criterion\":\"commuting descents and all descent predecessors FC; BFS closed under FC ascents\"}\n")
}

func main() {
	rank := flag.Int("rank", 8, "rank n of the generalized E_n (6..10)")
	flag.Parse()
	if flag.NArg() != 0 {
		fmt.Fprintln(os.Stderr, "fc-catalogue: positional arguments are not accepted; use -rank N")
		os.Exit(2)
	}
	if *rank < 6 || *rank > 10 {
		fmt.Fprintln(os.Stderr, "Supported ranks 6..10")
		os.Exit(2)
	}
	out := bufio.NewWriterSize(os.Stdout, 1<<20)
	newEngine(*rank).run(out)
	if err := out.Flush(); err != nil {
		fmt.Fprintln(os.Stderr, "fc-catalogue: write error:", err)
		os.Exit(1)
	}
}
