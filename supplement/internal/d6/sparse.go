package d6

import "fmt"

// MaxLetters bounds the number of signed letters of an element.
const MaxLetters = 12

// Perm is a signed one-line permutation (comparable, usable as a map key).
// Unused positions are zero.
type Perm struct {
	N int8
	V [MaxLetters]int8
}

// NewPerm builds a Perm from one-line notation.
func NewPerm(w []int64) (Perm, error) {
	var p Perm
	if len(w) > MaxLetters {
		return p, fmt.Errorf("element with %d letters exceeds %d", len(w), MaxLetters)
	}
	p.N = int8(len(w))
	for i, a := range w {
		if a < -127 || a > 127 {
			return p, fmt.Errorf("letter %d out of range", a)
		}
		p.V[i] = int8(a)
	}
	return p, nil
}

// Slice returns the one-line notation.
func (p Perm) Slice() []int {
	out := make([]int, p.N)
	for i := range out {
		out[i] = int(p.V[i])
	}
	return out
}

// Less is the lexicographic order of Python tuples of ints.
func (p Perm) Less(q Perm) bool {
	m := p.N
	if q.N < m {
		m = q.N
	}
	for i := int8(0); i < m; i++ {
		if p.V[i] != q.V[i] {
			return p.V[i] < q.V[i]
		}
	}
	return p.N < q.N
}

type pair struct{ x, w Perm }

// SparseCoxeter is the targeted signed-permutation Coxeter group of type
// A, B or D with the independent length formula and lifting-property Bruhat
// comparison of computations/sparse_kl.py. Caches mirror functools.cache.
type SparseCoxeter struct {
	Kind byte
	Rank int

	lengthC   map[Perm]int
	descentsC map[Perm][]int
	leqC      map[pair]bool
	lowerC    map[Perm]map[Perm]struct{}
	rpolyC    map[pair]Poly
}

// NewSparseCoxeter returns the group of the given kind ('A','B','D') and rank.
func NewSparseCoxeter(kind string, rank int) (*SparseCoxeter, error) {
	if kind != "A" && kind != "B" && kind != "D" {
		return nil, fmt.Errorf("unsupported group type %q", kind)
	}
	if rank < 1 {
		return nil, fmt.Errorf("bad rank %d", rank)
	}
	return &SparseCoxeter{
		Kind: kind[0], Rank: rank,
		lengthC:   map[Perm]int{},
		descentsC: map[Perm][]int{},
		leqC:      map[pair]bool{},
		lowerC:    map[Perm]map[Perm]struct{}{},
		rpolyC:    map[pair]Poly{},
	}, nil
}

// MinLetters is the number of letters an element needs for every generator
// of this group to act.
func (g *SparseCoxeter) MinLetters() int {
	if g.Kind == 'A' {
		return g.Rank + 1
	}
	if g.Kind == 'D' && g.Rank < 2 {
		return 2
	}
	return g.Rank
}

// Length is the independent signed inversion length formula.
func (g *SparseCoxeter) Length(w Perm) int {
	if v, ok := g.lengthC[w]; ok {
		return v
	}
	n := int(w.N)
	inv, neg := 0, 0
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if w.V[i] > w.V[j] {
				inv++
			}
			if -w.V[i] > w.V[j] {
				neg++
			}
		}
	}
	var res int
	if g.Kind == 'A' {
		res = inv
	} else {
		res = inv + neg
		if g.Kind == 'B' {
			for i := 0; i < n; i++ {
				if w.V[i] < 0 {
					res++
				}
			}
		}
	}
	g.lengthC[w] = res
	return res
}

// Right is right multiplication by the simple reflection s.
func (g *SparseCoxeter) Right(w Perm, s int) Perm {
	a := w
	switch {
	case g.Kind == 'A':
		a.V[s], a.V[s+1] = a.V[s+1], a.V[s]
	case s != 0:
		a.V[s-1], a.V[s] = a.V[s], a.V[s-1]
	case g.Kind == 'B':
		a.V[0] = -a.V[0]
	default:
		a.V[0], a.V[1] = -a.V[1], -a.V[0]
	}
	return a
}

// Descents lists the right descents of w in increasing order.
func (g *SparseCoxeter) Descents(w Perm) []int {
	if v, ok := g.descentsC[w]; ok {
		return v
	}
	out := []int{}
	lw := g.Length(w)
	for s := 0; s < g.Rank; s++ {
		if g.Length(g.Right(w, s)) < lw {
			out = append(out, s)
		}
	}
	g.descentsC[w] = out
	return out
}

// HasDescent reports whether s is a right descent of w.
func (g *SparseCoxeter) HasDescent(w Perm, s int) bool {
	for _, d := range g.Descents(w) {
		if d == s {
			return true
		}
	}
	return false
}

// Leq is the Bruhat order by the lifting property.
func (g *SparseCoxeter) Leq(x, w Perm) bool {
	if x == w {
		return true
	}
	if g.Length(x) >= g.Length(w) {
		return false
	}
	k := pair{x, w}
	if v, ok := g.leqC[k]; ok {
		return v
	}
	d := g.Descents(w)
	if len(d) == 0 {
		panic("non-identity element without right descent")
	}
	s := d[0]
	xs, ws := g.Right(x, s), g.Right(w, s)
	var res bool
	if g.Length(xs) < g.Length(x) {
		res = g.Leq(xs, ws)
	} else {
		res = g.Leq(x, ws)
	}
	g.leqC[k] = res
	return res
}

// Lower is the Bruhat lower ideal {x : x <= w}. The returned set must not be
// modified.
func (g *SparseCoxeter) Lower(w Perm) map[Perm]struct{} {
	if v, ok := g.lowerC[w]; ok {
		return v
	}
	var res map[Perm]struct{}
	if g.Length(w) == 0 {
		res = map[Perm]struct{}{w: {}}
	} else {
		d := g.Descents(w)
		if len(d) == 0 {
			panic("non-identity element without right descent")
		}
		s := d[0]
		base := g.Lower(g.Right(w, s))
		res = make(map[Perm]struct{}, 2*len(base))
		for x := range base {
			res[x] = struct{}{}
			res[g.Right(x, s)] = struct{}{}
		}
	}
	g.lowerC[w] = res
	return res
}

// RPoly is the Kazhdan-Lusztig R-polynomial R_{x,w} (ascending powers of q).
func (g *SparseCoxeter) RPoly(x, w Perm) Poly {
	if x == w {
		return Poly{1}
	}
	if !g.Leq(x, w) {
		return Poly{}
	}
	k := pair{x, w}
	if v, ok := g.rpolyC[k]; ok {
		return v
	}
	s := g.Descents(w)[0]
	v, xs := g.Right(w, s), g.Right(x, s)
	var p Poly
	if g.Length(xs) < g.Length(x) {
		p = g.RPoly(xs, v)
	} else {
		// q*R(xs,v) + q*R(x,v) - R(x,v)
		p = Add(Poly{}, g.RPoly(xs, v), 1, 1)
		p = Add(p, g.RPoly(x, v), 1, 1)
		p = Add(p, g.RPoly(x, v), 0, -1)
	}
	g.rpolyC[k] = p
	return p
}
