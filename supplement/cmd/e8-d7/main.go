// Command e8-d7: independent exhaustive enumeration of the right-and-left
// weak "terminal" elements of the Coxeter group E8 through cosets of the
// parabolic subgroup D7 (not E7).
//
// E8 is encoded as chain 0--1--2--3--4--5--6 with a branch node 7 attached to
// node 2; J={1,...,7} is D7. Elements are stored as the 8 root ids of the
// images of the simple roots, so the full E8 root action is retained while
// only D7 (322560 elements) is enumerated. Terminals w = c*v are built from
// the 2160 minimal left coset representatives c and the D7 elements v that are
// weak on the right in D7, and are tested for the weak condition on both sides.
// The result is cross-checked against an independently generated root-id
// fully-commutative catalogue (10846 elements) and printed as JSON on stdout.
//
// Ports research/en_independent/e8_d7_cosets.cpp. Imports only the Go
// standard library: no other engine of this module (in particular no
// internal/parabolic) is used, so the enumeration stays independent.
package main

import (
	"bufio"
	"fmt"
	"math/bits"
	"os"
	"sort"
)

type root [8]int

var (
	adj      [8][]int
	roots    []root
	rid      = map[root]int{}
	negs     [240]int
	positive [240]bool
	adds     [240][240]int
)

func fail(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "e8-d7: assertion failed: "+format+"\n", args...)
	os.Exit(1)
}

func must(cond bool, what string) {
	if !cond {
		fail("%s", what)
	}
}

func get(w uint64, s int) int { return int((w >> (8 * uint(s))) & 255) }

func set(w uint64, s int, a int) uint64 {
	return (w &^ (uint64(255) << (8 * uint(s)))) | (uint64(a) << (8 * uint(s)))
}

func right(w uint64, s int) uint64 {
	v := set(w, s, negs[get(w, s)])
	for _, t := range adj[s] {
		a := adds[get(w, s)][get(w, t)]
		must(a >= 0, "root sum in right multiplication")
		v = set(v, t, a)
	}
	return v
}

func desc(w uint64) uint {
	var d uint
	for s := 0; s < 8; s++ {
		if !positive[get(w, s)] {
			d |= 1 << uint(s)
		}
	}
	return d
}

func weak(w uint64, allowed uint) bool {
	d := desc(w) & allowed
	for s := 0; s < 8; s++ {
		if d&(1<<uint(s)) != 0 {
			for _, t := range adj[s] {
				if allowed&(1<<uint(t)) != 0 {
					sum := adds[get(w, s)][get(w, t)]
					must(sum >= 0, "root sum in weak test")
					if !positive[sum] {
						return false
					}
				}
			}
		}
	}
	return true
}

func independent(d uint) bool {
	for s := 0; s < 8; s++ {
		if d&(1<<uint(s)) != 0 {
			for _, t := range adj[s] {
				if d&(1<<uint(t)) != 0 {
					return false
				}
			}
		}
	}
	return true
}

func image(w uint64, a root) root {
	var b root
	for s := 0; s < 8; s++ {
		for k := 0; k < 8; k++ {
			b[k] += a[s] * roots[get(w, s)][k]
		}
	}
	return b
}

func permutation(w uint64) [240]uint8 {
	var p [240]uint8
	for r := 0; r < 240; r++ {
		id, ok := rid[image(w, roots[r])]
		must(ok, "image of a root is a root")
		p[r] = uint8(id)
	}
	return p
}

func compose(p *[240]uint8, v uint64) uint64 {
	var w uint64
	for s := 0; s < 8; s++ {
		w = set(w, s, int(p[get(v, s)]))
	}
	return w
}

func identity() uint64 {
	var e uint64
	for s := 0; s < 8; s++ {
		e = set(e, s, s)
	}
	return e
}

func element(word []int) uint64 {
	w := identity()
	for _, s := range word {
		w = right(w, s)
	}
	return w
}

func inverse(word []int) uint64 {
	w := identity()
	for i := len(word) - 1; i >= 0; i-- {
		w = right(w, word[i])
	}
	return w
}

// wordOf follows parent pointers from index w back to the identity (index 0)
// and returns the generators in application order.
func wordOf(w uint32, parent []uint32, last []uint8) []int {
	var a []int
	for w != 0 {
		a = append(a, int(last[w]))
		w = parent[w]
	}
	for i, j := 0, len(a)-1; i < j; i, j = i+1, j-1 {
		a[i], a[j] = a[j], a[i]
	}
	return a
}

func printWord(out *bufio.Writer, a []int) {
	out.WriteByte('[')
	for i, s := range a {
		if i > 0 {
			out.WriteByte(',')
		}
		fmt.Fprintf(out, "%d", s)
	}
	out.WriteByte(']')
}

func leq(x uint64, lx int, w uint64, lw int) bool {
	for x != w {
		if lx >= lw {
			return false
		}
		d := desc(w)
		must(d != 0, "nonidentity element has a descent in leq")
		s := bits.TrailingZeros(d)
		if desc(x)&(1<<uint(s)) != 0 {
			x = right(x, s)
			lx--
		}
		w = right(w, s)
		lw--
	}
	return true
}

func lessWord(a, b []int) bool {
	for i := 0; i < len(a) && i < len(b); i++ {
		if a[i] != b[i] {
			return a[i] < b[i]
		}
	}
	return len(a) < len(b)
}

type terminal struct {
	w, inv uint64
	length int
	wd     []int
}

func supportOf(wd []int) uint {
	var support uint
	for _, s := range wd {
		support |= 1 << uint(s)
	}
	return support
}

func initRoots() {
	for s := 0; s < 6; s++ {
		adj[s] = append(adj[s], s+1)
		adj[s+1] = append(adj[s+1], s)
	}
	adj[2] = append(adj[2], 7)
	adj[7] = append(adj[7], 2)
	for s := 0; s < 8; s++ {
		var a root
		a[s] = 1
		rid[a] = len(roots)
		roots = append(roots, a)
	}
	for i := 0; i < len(roots); i++ {
		for s := 0; s < 8; s++ {
			a := roots[i]
			a[s] = -a[s]
			for _, t := range adj[s] {
				a[s] += roots[i][t]
			}
			if _, ok := rid[a]; !ok {
				rid[a] = len(roots)
				roots = append(roots, a)
			}
		}
	}
	must(len(roots) == 240, "E8 has 240 roots")
	for i := 0; i < 240; i++ {
		a := roots[i]
		positive[i] = true
		for k := 0; k < 8; k++ {
			positive[i] = positive[i] && a[k] >= 0
			a[k] = -a[k]
		}
		id, ok := rid[a]
		must(ok, "negative of a root is a root")
		negs[i] = id
		for j := 0; j < 240; j++ {
			var b root
			for k := 0; k < 8; k++ {
				b[k] = roots[i][k] + roots[j][k]
			}
			if id, ok := rid[b]; ok {
				adds[i][j] = id
			} else {
				adds[i][j] = -1
			}
		}
	}

}

func main() {
	initRoots()

	// Enumerate only D7, retaining full E8 root actions.
	els := []uint64{identity()}
	ids := make(map[uint64]uint32, 400000)
	ids[identity()] = 0
	parent := []uint32{0}
	last := []uint8{0}
	length := []uint8{0}
	for i := 0; i < len(els); i++ {
		for s := 1; s < 8; s++ {
			v := right(els[i], s)
			if _, ok := ids[v]; !ok {
				j := uint32(len(els))
				ids[v] = j
				els = append(els, v)
				parent = append(parent, uint32(i))
				last = append(last, uint8(s))
				length = append(length, length[i]+1)
			}
		}
	}
	must(len(els) == 322560 && length[len(length)-1] == 42, "D7 order 322560 and longest length 42")
	var rightweak []uint32
	var d7hist [43]uint64
	for i := 0; i < len(els); i++ {
		d7hist[length[i]]++
		if weak(els[i], 254) {
			rightweak = append(rightweak, uint32(i))
		}
	}

	// The fundamental weight at deleted node0 is a norm-four integral vector.
	beta := root{4, 7, 10, 8, 6, 4, 2, 5}
	for s := 0; s < 8; s++ {
		pairing := 2 * beta[s]
		for _, t := range adj[s] {
			pairing -= beta[t]
		}
		want := 0
		if s == 0 {
			want = 1
		}
		must(pairing == want, "beta is the fundamental weight at node 0")
	}
	orbit := []root{beta}
	orbitIDs := map[root]uint32{beta: 0}
	cp := []uint32{0}
	cg := []uint8{0}
	cl := []uint8{0}
	for i := 0; i < len(orbit); i++ {
		for s := 0; s < 8; s++ {
			b := orbit[i]
			b[s] = -b[s]
			for _, t := range adj[s] {
				b[s] += orbit[i][t]
			}
			if _, ok := orbitIDs[b]; !ok {
				orbitIDs[b] = uint32(len(orbit))
				orbit = append(orbit, b)
				cp = append(cp, uint32(i))
				cg = append(cg, uint8(s))
				cl = append(cl, cl[i]+1)
			}
		}
	}
	must(len(orbit) == 2160 && cl[len(cl)-1] == 78, "orbit of beta has 2160 points and depth 78")
	var cosets, cosetinv []uint64
	var cosetword [][]int
	var cosethist [79]uint64
	var cperm [][240]uint8
	for i := 0; i < len(orbit); i++ {
		// Left BFS: generators are recovered latest-first, with no reversal.
		var a []int
		for j := uint32(i); j != 0; j = cp[j] {
			a = append(a, int(cg[j]))
		}
		w := element(a)
		must(image(w, beta) == orbit[i], "coset representative maps beta to its orbit point")
		must(desc(w)&254 == 0, "coset representative has no descent in D7")
		cosetword = append(cosetword, a)
		cosets = append(cosets, w)
		cosetinv = append(cosetinv, inverse(a))
		cperm = append(cperm, permutation(w))
		cosethist[cl[i]]++
	}

	// Independent root-ID FC catalogue for checking all eligible bottoms.
	fcs := []uint64{identity()}
	fid := map[uint64]uint32{identity(): 0}
	fp := []uint32{0}
	fg := []uint8{0}
	fl := []uint8{0}
	for i := 0; i < len(fcs); i++ {
		a := fcs[i]
		ad := desc(a)
		for s := 0; s < 8; s++ {
			if ad&(1<<uint(s)) != 0 {
				continue
			}
			b := right(a, s)
			if _, ok := fid[b]; ok {
				continue
			}
			d := desc(b)
			if !independent(d) {
				continue
			}
			good := true
			for t := 0; t < 8; t++ {
				if d&(1<<uint(t)) != 0 {
					idx, ok := fid[right(b, t)]
					if !ok || fl[idx] != fl[i] {
						good = false
						break
					}
				}
			}
			if good {
				fid[b] = uint32(len(fcs))
				fcs = append(fcs, b)
				fp = append(fp, uint32(i))
				fg = append(fg, uint8(s))
				fl = append(fl, fl[i]+1)
			}
		}
	}
	must(len(fcs) == 10846, "FC catalogue has 10846 elements")
	fleft := make([]uint, len(fcs))
	for i := range fcs {
		fleft[i] = desc(inverse(wordOf(uint32(i), fp, fg)))
	}

	var terminals []terminal
	var tested, rightsurvive uint64
	for _, vi := range rightweak {
		vw := wordOf(vi, parent, last)
		vinv := inverse(vw)
		viperm := permutation(vinv)
		for ai := 0; ai < len(cosets); ai++ {
			tested++
			w := compose(&cperm[ai], els[vi])
			if !weak(w, 255) {
				continue
			}
			rightsurvive++
			inv := compose(&viperm, cosetinv[ai])
			if !weak(inv, 255) {
				continue
			}
			wd := append(append([]int(nil), cosetword[ai]...), vw...)
			must(element(wd) == w, "terminal word evaluates to w")
			must(len(wd) == int(cl[ai])+int(length[vi]), "terminal word length is coset length plus D7 length")
			terminals = append(terminals, terminal{w, inv, len(wd), wd})
		}
	}
	sort.Slice(terminals, func(i, j int) bool {
		a, b := &terminals[i], &terminals[j]
		if a.length != b.length {
			return a.length < b.length
		}
		return lessWord(a.wd, b.wd)
	})
	commuting := 0
	for i := range terminals {
		t := &terminals[i]
		support := supportOf(t.wd)
		if bits.OnesCount(support) == t.length && independent(support) {
			commuting++
		}
	}

	out := bufio.NewWriter(os.Stdout)
	defer out.Flush()
	fmt.Fprintf(out, "{\"rank\":8,\"parabolic\":\"D7\",\"parabolic_order\":%d,\"cosets\":%d,\"parabolic_right_terminals\":%d,"+
		"\"candidates_tested\":%d,\"full_right_terminals\":%d,\"terminal_count\":%d,\"commuting_terminals\":%d,"+
		"\"fc_count\":%d,\"complete\":true,\"parabolic_histogram\":[",
		len(els), len(cosets), len(rightweak), tested, rightsurvive, len(terminals), commuting, len(fcs))
	for k, v := range d7hist {
		if k > 0 {
			out.WriteByte(',')
		}
		fmt.Fprintf(out, "%d", v)
	}
	out.WriteString("],\"coset_histogram\":[")
	for k, v := range cosethist {
		if k > 0 {
			out.WriteByte(',')
		}
		fmt.Fprintf(out, "%d", v)
	}
	out.WriteString("],\"bad\":[")
	first := true
	for i := range terminals {
		t := &terminals[i]
		support := supportOf(t.wd)
		if bits.OnesCount(support) == t.length && independent(support) {
			continue
		}
		if !first {
			out.WriteByte(',')
		}
		first = false
		out.WriteString("\n{\"word\":")
		printWord(out, t.wd)
		rd, ld := desc(t.w), desc(t.inv)
		fmt.Fprintf(out, ",\"length\":%d,\"Rmask\":%d,\"Lmask\":%d,\"bottoms\":[", t.length, rd, ld)
		ff := true
		for x := 0; x < len(fcs); x++ {
			if desc(fcs[x])&rd == rd && fleft[x]&ld == ld && leq(fcs[x], int(fl[x]), t.w, t.length) {
				if !ff {
					out.WriteByte(',')
				}
				ff = false
				out.WriteString("{\"word\":")
				printWord(out, wordOf(uint32(x), fp, fg))
				fmt.Fprintf(out, ",\"length\":%d,\"rank\":%d}", int(fl[x]), t.length-int(fl[x]))
			}
		}
		out.WriteString("]}")
	}
	out.WriteString("\n]}\n")
}
