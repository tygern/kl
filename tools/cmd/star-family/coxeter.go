package main

import (
	"encoding/binary"
	"fmt"
	"math/bits"
)

// Exact simply-laced Coxeter group arithmetic on integral matrices (the
// geometric representation, stored by columns), following star()/Coxeter of
// research/broad_cells/targeted.py. Group elements are interned: equal
// matrices have the same identifier, so identity of identifiers is equality.

type elem struct {
	m      []int32 // column-major, n*n entries
	hash   uint64  // CPython hash of the tuple-of-column-tuples
	right  []int32 // right[s] = element after right multiplication by s, -1 if unknown
	desc   []int
	hasDsc bool
	length int // -1 if unknown
	inv    int32
	fc     int8 // 0 unknown, 1 true, 2 false
}

type pair struct{ a, b int32 }

type pairS struct {
	w int32
	s int
}

type correction struct {
	z     int32
	d, mu int
}

type coxeter struct {
	n       int
	edges   [][2]int
	adj     [][]int
	elems   []*elem
	index   map[string]int32
	e       int32
	klVals  map[pair][]int64
	klOrder []pair
	lowers  map[int32]*pySet
	corr    map[pairS][]correction
	keyBuf  []byte
}

func newCoxeter(n int, edges [][2]int) *coxeter {
	g := &coxeter{n: n, edges: edges, index: map[string]int32{}, klVals: map[pair][]int64{},
		lowers: map[int32]*pySet{}, corr: map[pairS][]correction{}}
	g.adj = make([][]int, n)
	for _, e := range edges {
		g.adj[e[0]] = append(g.adj[e[0]], e[1])
		g.adj[e[1]] = append(g.adj[e[1]], e[0])
	}
	id := make([]int32, n*n)
	for j := 0; j < n; j++ {
		id[j*n+j] = 1
	}
	g.e = g.intern(id)
	return g
}

func star(m int) *coxeter {
	edges := make([][2]int, 0, m)
	for i := 1; i <= m; i++ {
		edges = append(edges, [2]int{0, i})
	}
	return newCoxeter(m+1, edges)
}

// CPython hashes (64-bit build, Python >= 3.8).
const (
	xxPrime1 uint64 = 11400714785074694791
	xxPrime2 uint64 = 14029467366897019727
	xxPrime5 uint64 = 2870177450012600261
)

func tupleHash(lanes []uint64) uint64 {
	acc := xxPrime5
	for _, l := range lanes {
		acc += l * xxPrime2
		acc = bits.RotateLeft64(acc, 31)
		acc *= xxPrime1
	}
	acc += uint64(len(lanes)) ^ (xxPrime5 ^ 3527539)
	if acc == ^uint64(0) {
		return 1546275796
	}
	return acc
}

func intHash(v int64) uint64 {
	const modulus = (1 << 61) - 1
	neg := v < 0
	a := uint64(v)
	if neg {
		a = uint64(-v)
	}
	h := int64(a % modulus)
	if neg {
		h = -h
	}
	if h == -1 {
		h = -2
	}
	return uint64(h)
}

func (g *coxeter) matrixHash(m []int32) uint64 {
	n := g.n
	cols := make([]uint64, n)
	lanes := make([]uint64, n)
	for j := 0; j < n; j++ {
		for i := 0; i < n; i++ {
			lanes[i] = intHash(int64(m[j*n+i]))
		}
		cols[j] = tupleHash(lanes)
	}
	return tupleHash(cols)
}

func (g *coxeter) intern(m []int32) int32 {
	buf := g.keyBuf[:0]
	for _, v := range m {
		buf = binary.LittleEndian.AppendUint32(buf, uint32(v))
	}
	g.keyBuf = buf
	if id, ok := g.index[string(buf)]; ok {
		return id
	}
	e := &elem{m: append([]int32(nil), m...), hash: g.matrixHash(m), length: -1, inv: -1}
	e.right = make([]int32, g.n)
	for i := range e.right {
		e.right[i] = -1
	}
	id := int32(len(g.elems))
	g.elems = append(g.elems, e)
	g.index[string(buf)] = id
	return id
}

func checked32(v int64) int32 {
	if v > 1<<30 || v < -(1<<30) {
		panic(fmt.Sprintf("matrix entry %d outside the guarded range", v))
	}
	return int32(v)
}

func (g *coxeter) right(w int32, s int) int32 {
	e := g.elems[w]
	if r := e.right[s]; r >= 0 {
		return r
	}
	n := g.n
	out := append([]int32(nil), e.m...)
	for i := 0; i < n; i++ {
		out[s*n+i] = -e.m[s*n+i]
	}
	for _, t := range g.adj[s] {
		for i := 0; i < n; i++ {
			out[t*n+i] = checked32(int64(e.m[s*n+i]) + int64(e.m[t*n+i]))
		}
	}
	r := g.intern(out)
	g.elems[w].right[s] = r
	return r
}

// desc returns the right descent set in increasing order: the columns whose
// entries are all <= 0.
func (g *coxeter) desc(w int32) []int {
	e := g.elems[w]
	if e.hasDsc {
		return e.desc
	}
	var ds []int
	for s := 0; s < g.n; s++ {
		all := true
		for i := 0; i < g.n; i++ {
			if e.m[s*g.n+i] > 0 {
				all = false
				break
			}
		}
		if all {
			ds = append(ds, s)
		}
	}
	e.desc = ds
	e.hasDsc = true
	return ds
}

func contains(xs []int, v int) bool {
	for _, x := range xs {
		if x == v {
			return true
		}
	}
	return false
}

// peel returns the letters removed one at a time by right multiplication by the
// smallest descent, from w down to the identity. The reduced word of w is this
// sequence reversed.
func (g *coxeter) peel(w int32) []int {
	var out []int
	for w != g.e {
		ds := g.desc(w)
		if len(ds) == 0 {
			panic("non-identity element without right descent")
		}
		out = append(out, ds[0])
		w = g.right(w, ds[0])
	}
	return out
}

func (g *coxeter) word(w int32) []int {
	p := g.peel(w)
	for i, j := 0, len(p)-1; i < j; i, j = i+1, j-1 {
		p[i], p[j] = p[j], p[i]
	}
	return p
}

func (g *coxeter) length(w int32) int {
	e := g.elems[w]
	if e.length >= 0 {
		return e.length
	}
	if w == g.e {
		e.length = 0
		return 0
	}
	ds := g.desc(w)
	if len(ds) == 0 {
		panic("non-identity element without right descent")
	}
	l := g.length(g.right(w, ds[0])) + 1
	g.elems[w].length = l
	return l
}

func (g *coxeter) elt(word []int) int32 {
	w := g.e
	for _, s := range word {
		w = g.right(w, s)
	}
	return w
}

func (g *coxeter) inv(w int32) int32 {
	e := g.elems[w]
	if e.inv >= 0 {
		return e.inv
	}
	r := g.elt(g.peel(w)) // elt(reversed(word(w)))
	g.elems[w].inv = r
	return r
}

func (g *coxeter) fc(w int32) bool {
	e := g.elems[w]
	if e.fc != 0 {
		return e.fc == 1
	}
	ds := g.desc(w)
	ok := true
	for _, s := range ds {
		for _, t := range g.adj[s] {
			if contains(ds, t) {
				ok = false
			}
		}
	}
	if ok {
		for _, s := range ds {
			if !g.fc(g.right(w, s)) {
				ok = false
				break
			}
		}
	}
	if ok {
		g.elems[w].fc = 1
	} else {
		g.elems[w].fc = 2
	}
	return ok
}

func (g *coxeter) terminal(w int32) bool {
	for _, v := range []int32{w, g.inv(w)} {
		for _, s := range g.desc(v) {
			ds := g.desc(g.right(v, s))
			for _, t := range g.adj[s] {
				if contains(ds, t) {
					return false
				}
			}
		}
	}
	return true
}

// lower returns the Bruhat lower ideal as a frozenset-ordered set, built exactly
// as the original: base | {right(z, s) for z in base}.
func (g *coxeter) lower(w int32) *pySet {
	if s, ok := g.lowers[w]; ok {
		return s
	}
	var res *pySet
	if w == g.e {
		res = newPySet()
		res.add(g.e, g.elems[g.e].hash)
	} else {
		s := g.desc(w)[0]
		base := g.lower(g.right(w, s))
		comp := newPySet()
		for _, z := range base.order() {
			r := g.right(z, s)
			comp.add(r, g.elems[r].hash)
		}
		res = newPySet()
		res.merge(base)
		res.merge(comp)
	}
	g.lowers[w] = res
	return res
}

func (g *coxeter) leq(x, w int32) bool {
	for x != w {
		if g.length(x) >= g.length(w) {
			return false
		}
		s := g.desc(w)[0]
		if contains(g.desc(x), s) {
			x = g.right(x, s)
		}
		w = g.right(w, s)
	}
	return true
}

// eligible lists the fully commutative elements of the lower ideal of w whose
// left and right descent sets contain those of w.
func (g *coxeter) eligible(w int32) []int32 {
	rd := g.desc(w)
	ld := g.desc(g.inv(w))
	var out []int32
	for _, x := range g.lower(w).order() {
		if !g.fc(x) {
			continue
		}
		if !subset(rd, g.desc(x)) || !subset(ld, g.desc(g.inv(x))) {
			continue
		}
		out = append(out, x)
	}
	return out
}

func subset(a, b []int) bool {
	for _, x := range a {
		if !contains(b, x) {
			return false
		}
	}
	return true
}

// addPoly is add(a, b, shift, scale) of the original: a + scale*q^shift*b,
// trailing zeros removed.
func addPoly(a, b []int64, shift int, scale int64) []int64 {
	n := len(a)
	if len(b)+shift > n {
		n = len(b) + shift
	}
	c := make([]int64, n)
	copy(c, a)
	for i, v := range b {
		t := scale * v
		if v != 0 && t/v != scale {
			panic("polynomial coefficient overflow")
		}
		sum := c[i+shift] + t
		if (sum > c[i+shift]) != (t > 0) && t != 0 {
			panic("polynomial coefficient overflow")
		}
		c[i+shift] = sum
	}
	for len(c) > 0 && c[len(c)-1] == 0 {
		c = c[:len(c)-1]
	}
	return c
}

func (g *coxeter) corrections(w int32, s int) []correction {
	key := pairS{w, s}
	if c, ok := g.corr[key]; ok {
		return c
	}
	v := g.right(w, s)
	out := []correction{}
	lv := g.length(v)
	for _, z := range g.lower(v).order() {
		d := lv - g.length(z)
		if d <= 0 || d%2 == 0 || !contains(g.desc(z), s) {
			continue
		}
		p := g.kl(z, v)
		mu := int64(0)
		if len(p) > (d-1)/2 {
			mu = p[(d-1)/2]
		}
		if mu != 0 {
			out = append(out, correction{z, (d + 1) / 2, int(mu)})
		}
	}
	g.corr[key] = out
	return out
}

func (g *coxeter) kl(x, w int32) []int64 {
	k := pair{x, w}
	if p, ok := g.klVals[k]; ok {
		return p
	}
	var p []int64
	switch {
	case x == w:
		p = []int64{1}
	case !g.leq(x, w):
		p = []int64{}
	default:
		p = g.withDescent(x, w, g.desc(w)[0])
	}
	g.klVals[k] = p
	g.klOrder = append(g.klOrder, k)
	return p
}

func (g *coxeter) withDescent(x, w int32, s int) []int64 {
	v := g.right(w, s)
	xs := g.right(x, s)
	c := 0
	if contains(g.desc(x), s) {
		c = 1
	}
	p := addPoly(nil, g.kl(xs, v), 1-c, 1)
	p = addPoly(p, g.kl(x, v), c, 1)
	for _, t := range g.corrections(w, s) {
		if g.leq(x, t.z) {
			p = addPoly(p, g.kl(x, t.z), t.d, -int64(t.mu))
		}
	}
	return p
}
