// Command enumerate-bad exhaustively enumerates the Weyl group of E6/E7
// (E8 is available but very expensive) by its action on the exact root
// system (roots as integer coordinate vectors in the simple-root basis),
// computes right descent sets, the fully commutative (FC) elements, and the
// "bad" elements: non-FC elements whose right and left descents are both
// weakly commuting-free (weakright). For each bad element it lists the FC
// "bottoms" below it with the same descent masks. The JSON certificate is
// written to standard output.
//
// It ports research/broad_exceptional/enumerate_bad.cpp. Diagram: nodes
// 0--1--2--3--...--(n-2), with node n-1 attached to node 2.
//
// This is the root-index engine. It is self-contained: it imports only the
// Go standard library and no internal package, and shares no code with the
// integer-matrix engine (cmd/matrix-search).
package main

import (
	"bufio"
	"flag"
	"fmt"
	"math/bits"
	"os"
	"strconv"
)

type root [8]int

var (
	n      int
	adj    [8][]int
	roots  []root
	ri     = map[root]int{}
	negs   [256]int
	pos    [256]bool
	adds   [256][256]int
	elems  []uint64
	rights [][8]uint32
	parent []uint32
	lastg  []uint8
	lens   []uint8
	ds     []uint8
	fc     []bool
)

func fail(format string, a ...any) {
	fmt.Fprintf(os.Stderr, "enumerate-bad: assertion failed: "+format+"\n", a...)
	os.Exit(1)
}

func check(cond bool, format string, a ...any) {
	if !cond {
		fail(format, a...)
	}
}

func rget(x uint64, s int) int { return int((x >> (8 * uint(s))) & 255) }

func rset(x uint64, s int, v int) uint64 {
	return (x &^ (uint64(255) << (8 * uint(s)))) | (uint64(v) << (8 * uint(s)))
}

func mult(w uint64, s int) uint64 {
	v := rset(w, s, negs[rget(w, s)])
	for _, t := range adj[s] {
		a := adds[rget(w, s)][rget(w, t)]
		check(a >= 0, "root sum is not a root")
		v = rset(v, t, a)
	}
	return v
}

func word(i uint32) []int {
	var w []int
	for i != 0 {
		w = append(w, int(lastg[i]))
		i = parent[i]
	}
	for a, b := 0, len(w)-1; a < b; a, b = a+1, b-1 {
		w[a], w[b] = w[b], w[a]
	}
	return w
}

func inv(i uint32) uint32 {
	w := word(i)
	var v uint32
	for k := len(w) - 1; k >= 0; k-- {
		v = rights[v][w[k]]
	}
	return v
}

func weakright(i uint32) bool {
	for s := 0; s < n; s++ {
		if ds[i]&(1<<uint(s)) != 0 {
			for _, t := range adj[s] {
				if ds[rights[i][s]]&(1<<uint(t)) != 0 {
					return false
				}
			}
		}
	}
	return true
}

func leq(x, w uint32) bool {
	for x != w {
		if lens[x] >= lens[w] {
			return false
		}
		s := bits.TrailingZeros8(ds[w])
		if ds[x]&(1<<uint(s)) != 0 {
			x = rights[x][s]
		}
		w = rights[w][s]
	}
	return true
}

func printword(out *bufio.Writer, i uint32) {
	w := word(i)
	out.WriteByte('[')
	for j, s := range w {
		if j > 0 {
			out.WriteByte(',')
		}
		out.WriteString(strconv.Itoa(s))
	}
	out.WriteByte(']')
}

func main() {
	rank := flag.Int("rank", 7, "rank of the exceptional group (6, 7 or 8)")
	flag.Parse()
	n = *rank
	check(n >= 6 && n <= 8, "rank must be in 6..8")
	for s := 0; s < n-2; s++ {
		adj[s] = append(adj[s], s+1)
		adj[s+1] = append(adj[s+1], s)
	}
	adj[2] = append(adj[2], n-1)
	adj[n-1] = append(adj[n-1], 2)

	for s := 0; s < n; s++ {
		var a root
		a[s] = 1
		ri[a] = len(roots)
		roots = append(roots, a)
	}
	for i := 0; i < len(roots); i++ {
		for s := 0; s < n; s++ {
			a := roots[i]
			a[s] = -a[s]
			for _, t := range adj[s] {
				a[s] += roots[i][t]
			}
			if _, ok := ri[a]; !ok {
				ri[a] = len(roots)
				roots = append(roots, a)
			}
		}
	}
	nroots := map[int]int{6: 72, 7: 126, 8: 240}[n]
	check(len(roots) == nroots, "root count %d, expected %d", len(roots), nroots)
	for i := 0; i < len(roots); i++ {
		a := roots[i]
		pos[i] = true
		for s := 0; s < n; s++ {
			pos[i] = pos[i] && a[s] >= 0
			a[s] = -a[s]
		}
		negs[i] = ri[a]
		for j := 0; j < len(roots); j++ {
			var b root
			for s := 0; s < n; s++ {
				b[s] = roots[i][s] + roots[j][s]
			}
			if k, ok := ri[b]; ok {
				adds[i][j] = k
			} else {
				adds[i][j] = -1
			}
		}
	}
	expected := map[int]int{6: 51840, 7: 2903040, 8: 696729600}[n]
	capHint := expected
	if capHint > 3000000 {
		capHint = 3000000
	}
	indexer := make(map[uint64]uint32, capHint)
	elems = make([]uint64, 0, capHint)
	rights = make([][8]uint32, 0, capHint)
	var e uint64
	for s := 0; s < n; s++ {
		e = rset(e, s, s)
	}
	elems = append(elems, e)
	indexer[e] = 0
	parent = append(parent, 0)
	lastg = append(lastg, 0)
	lens = append(lens, 0)
	for i := 0; i < len(elems); i++ {
		var row [8]uint32
		var d uint8
		for s := 0; s < n; s++ {
			v := mult(elems[i], s)
			k, ok := indexer[v]
			if !ok {
				k = uint32(len(elems))
				indexer[v] = k
				elems = append(elems, v)
				parent = append(parent, uint32(i))
				lastg = append(lastg, uint8(s))
				lens = append(lens, lens[i]+1)
			}
			row[s] = k
			if !pos[rget(elems[i], s)] {
				d |= 1 << uint(s)
			}
		}
		rights = append(rights, row)
		ds = append(ds, d)
		f := true
		for s := 0; s < n; s++ {
			if d&(1<<uint(s)) != 0 {
				f = f && fc[row[s]]
				for _, t := range adj[s] {
					f = f && d&(1<<uint(t)) == 0
				}
			}
		}
		fc = append(fc, f)
	}
	check(len(elems) == expected, "order %d, expected %d", len(elems), expected)

	var commutingWeak, weakfcNoncommuting, starChecks int64
	maxLen := int(lens[len(lens)-1])
	hist := make([]int64, maxLen+1)
	for i := uint32(0); int(i) < len(elems); i++ {
		hist[lens[i]]++
		for s := 0; s < n; s++ {
			l1 := int(lens[rights[i][s]])
			l0 := int(lens[i])
			dl := l1 - l0
			check(dl == 1 || dl == -1, "length changes by %d", dl)
			check((ds[i]&(1<<uint(s)) != 0) == (l1 < l0), "descent/length mismatch")
			check(rights[rights[i][s]][s] == i, "generator not an involution")
		}
		if fc[i] {
			if weakright(i) && weakright(inv(i)) {
				wd := word(i)
				comm := true
				present := map[int]bool{}
				for _, s := range wd {
					present[s] = true
				}
				for _, s := range wd {
					for _, t := range adj[s] {
						if present[t] {
							comm = false
						}
					}
				}
				comm = comm && len(present) == len(wd)
				if comm {
					commutingWeak++
				} else {
					weakfcNoncommuting++
				}
			}
			for s := 0; s < n; s++ {
				for _, t := range adj[s] {
					if s < t && (ds[i]&(1<<uint(s)) != 0) != (ds[i]&(1<<uint(t)) != 0) {
						a, b := rights[i][s], rights[i][t]
						da := (ds[a]&(1<<uint(s)) != 0) != (ds[a]&(1<<uint(t)) != 0)
						db := (ds[b]&(1<<uint(s)) != 0) != (ds[b]&(1<<uint(t)) != 0)
						check(da != db, "star operation not unique")
						pick := b
						if da {
							pick = a
						}
						check(fc[pick], "star image not fully commutative")
						starChecks++
					}
				}
			}
		}
	}
	check(weakfcNoncommuting == 0, "weakly-commuting FC element that is not commuting")

	var fcs, bads []uint32
	lds := make([]uint8, len(elems))
	for i := uint32(0); int(i) < len(elems); i++ {
		if fc[i] {
			fcs = append(fcs, i)
			lds[i] = ds[inv(i)]
		}
		if !fc[i] && weakright(i) && weakright(inv(i)) {
			bads = append(bads, i)
		}
	}
	out := bufio.NewWriterSize(os.Stdout, 1<<20)
	defer out.Flush()
	fmt.Fprintf(out, "{\"rank\":%d,\"roots\":%d,\"order\":%d,\"fc_count\":%d,\"max_length\":%d,\"commuting_weak\":%d,\"fc_right_star_checks\":%d,\"length_distribution\":[",
		n, len(roots), len(elems), len(fcs), maxLen, commutingWeak, starChecks)
	for k, h := range hist {
		if k > 0 {
			out.WriteByte(',')
		}
		out.WriteString(strconv.FormatInt(h, 10))
	}
	out.WriteString("],\"bad\":[\n")
	first := true
	for _, w := range bads {
		if !first {
			out.WriteString(",\n")
		}
		first = false
		ld := ds[inv(w)]
		fmt.Fprintf(out, "{\"id\":%d,\"word\":", w)
		printword(out, w)
		fmt.Fprintf(out, ",\"length\":%d,\"Rmask\":%d,\"Lmask\":%d,\"bottoms\":[", lens[w], ds[w], ld)
		ff := true
		for _, x := range fcs {
			if ds[x]&ds[w] == ds[w] && lds[x]&ld == ld && leq(x, w) {
				if !ff {
					out.WriteByte(',')
				}
				ff = false
				fmt.Fprintf(out, "{\"id\":%d,\"word\":", x)
				printword(out, x)
				fmt.Fprintf(out, ",\"length\":%d,\"rank\":%d}", lens[x], int(lens[w])-int(lens[x]))
			}
		}
		out.WriteString("]}")
	}
	out.WriteString("\n]}\n")
}
