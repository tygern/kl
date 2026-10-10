// Command e9-quotient-kl computes Kazhdan-Lusztig polynomials P_{x,w} for a
// fixed w in a simply laced Coxeter group, in the parabolic quotient W/W_I for
// a commuting set I subset R(w), in exact integer arithmetic with overflow
// checks.  It produces the affine E8 (= E9) odd-gap certificate of the
// manuscript: the length-33 reflection r_beta, beta = (1,2,3,3,2,2,1,1,2), and
// its two length-34 right extensions w s_0 and w s_1, whose eligible fully
// commutative lower endpoints have odd gaps. For these odd gaps, vanishing
// cannot be concluded from parity, and the polynomials must be computed.
//
// Port of research/review_checks/e9_quotient_kl.cpp (C++), keeping its
// structure, its assertions, its self-validation (random elements of A4, D4
// and E6 with random commuting I subset R(w), compared pair by pair with a
// naive ordinary Kazhdan-Lusztig recursion, using the same Mersenne
// twister mt19937 so the pair counts match) and its output text byte for byte.
//
// Method.  Elements are integer matrices on the span of the simple roots (row
// j is w(alpha_j) in the simple-root basis); right multiplication acts on
// rows, left multiplication on the coordinates.  For I subset R(w) write
// w = w' w_I with w' in W^I.  The lower interval [e,w] is a union of cosets
// x' W_I, and Pq(x,y) := P_{x w_I, y w_I} for x,y in W^I satisfies the
// recursion obtained by restricting the left recursion (Kazhdan-Lusztig 1979,
// (2.2.c)) to the left ideal H C'_{w_I}; see the original source for the
// formulas.  Every polynomial is checked for constant term 1, nonnegative
// coefficients and the degree bound.
//
// Imports: the Go standard library only.  This command is self-contained: it
// imports no other package of this module and shares no code with any other
// engine.
//
// Usage:
//
//	e9-quotient-kl -mode certify -out results/e9-odd-gap-certificate.json
//	e9-quotient-kl -mode validate -type E6 -seed 3 -count 120 -maxlen 14
//
// Hidden debugging mode (differential testing only, no certificate output):
// e9-quotient-kl -mode pairs -type T -pairs FILE reads lines "xword wword"
// (strings of generator labels, "e" for the identity) and prints one line
// "xword wword Pq Pord" per pair, where Pq is P_{x,w} computed by the
// quotient engine as in certify, with I the greedy maximal commuting subset
// of R(w) taken in ascending label order (certify uses I = R(w), which is
// commuting for the certified elements; x is replaced by its minimal coset
// representative, which leaves P_{x,w} unchanged since I subset R(w)), and
// Pord is the value of the naive ordinary recursion
// (Ord) when the lower ideal of w has at most 6000 elements, "-" otherwise.
// "[]" denotes the zero polynomial.  The certify and validate modes are not
// affected by this addition.
//
// E_n labelling: chain 0-1-...-(n-2), node n-1 attached to node 2.  A string
// of labels denotes the product of the simple reflections in the order
// written; W acts on the left and the matrix of w has columns w(alpha_j).
package main

import (
	"bufio"
	"flag"
	"fmt"
	"math"
	"math/bits"
	"os"
	"slices"
	"sort"
	"strconv"
	"strings"
)

// Mat is the matrix of an element: m[j*9+i] = coefficient of alpha_i in
// w(alpha_j).  The array is comparable and serves directly as a map key.
type Mat [81]int8

var N int
var ADJ [][]int

func fatal(format string, a ...any) {
	fmt.Fprintf(os.Stderr, format+"\n", a...)
	os.Exit(3)
}

func identity() Mat {
	var m Mat
	for j := 0; j < N; j++ {
		m[j*9+j] = 1
	}
	return m
}

func chk(v int) int8 {
	if v > 120 || v < -120 {
		fatal("int8 overflow")
	}
	return int8(v)
}

func rmul(w Mat, s int) Mat {
	r := w
	for i := 0; i < N; i++ {
		r[s*9+i] = chk(-int(w[s*9+i]))
	}
	for _, t := range ADJ[s] {
		for i := 0; i < N; i++ {
			r[t*9+i] = chk(int(w[t*9+i]) + int(w[s*9+i]))
		}
	}
	return r
}

func lmul(w Mat, s int) Mat {
	r := w
	for j := 0; j < N; j++ {
		p := 2 * int(w[j*9+s])
		for _, t := range ADJ[s] {
			p -= int(w[j*9+t])
		}
		r[j*9+s] = chk(int(w[j*9+s]) - p)
	}
	return r
}

func negcol(w *Mat, j int) bool {
	nz := false
	for i := 0; i < N; i++ {
		if w[j*9+i] > 0 {
			return false
		}
		if w[j*9+i] < 0 {
			nz = true
		}
	}
	if !nz {
		fatal("zero column")
	}
	return true
}

func rword(w Mat) []int {
	out := []int{}
	e := identity()
	for w != e {
		s := -1
		for j := 0; j < N; j++ {
			if negcol(&w, j) {
				s = j
				break
			}
		}
		if s < 0 {
			fatal("no descent")
		}
		w = rmul(w, s)
		out = append(out, s)
	}
	slices.Reverse(out)
	return out
}

func fromword(wd []int) Mat {
	w := identity()
	for _, s := range wd {
		w = rmul(w, s)
	}
	return w
}

func inverse(w Mat) Mat {
	wd := rword(w)
	slices.Reverse(wd)
	return fromword(wd)
}

func Rmask(w Mat) int {
	m := 0
	for j := 0; j < N; j++ {
		if negcol(&w, j) {
			m |= 1 << j
		}
	}
	return m
}

func Lmask(w Mat) int { return Rmask(inverse(w)) }

func pairing(a, b []int) int {
	p := 0
	for i := 0; i < N; i++ {
		p += 2 * a[i] * b[i]
		for _, t := range ADJ[i] {
			p -= a[i] * b[t]
		}
	}
	return p
}

func reflection(beta []int) Mat {
	m := identity()
	for j := 0; j < N; j++ {
		aj := make([]int, N)
		aj[j] = 1
		p := pairing(aj, beta)
		for i := 0; i < N; i++ {
			b := 0
			if i == j {
				b = 1
			}
			m[j*9+i] = chk(b - p*beta[i])
		}
	}
	return m
}

// Full commutativity: w is FC iff every element of its right weak lower
// interval has pairwise commuting right descents (Stembridge 1996, Prop. 2.1:
// a reduced word with a factor sts, m(s,t)=3, exists iff some such element has
// two adjacent right descents).
var FCMEMO = map[Mat]bool{}

func isFC(w Mat) bool {
	if v, ok := FCMEMO[w]; ok {
		return v
	}
	R := Rmask(w)
	ok := true
	for s := 0; s < N && ok; s++ {
		if R>>s&1 == 1 {
			for _, t := range ADJ[s] {
				if R>>t&1 == 1 {
					ok = false
				}
			}
			if ok && !isFC(rmul(w, s)) {
				ok = false
			}
		}
	}
	FCMEMO[w] = ok
	return ok
}

// Poly holds ascending coefficients.
type Poly = []int64

func trim(p Poly) Poly {
	for len(p) > 0 && p[len(p)-1] == 0 {
		p = p[:len(p)-1]
	}
	return p
}

func addchk(a, b int64) int64 {
	c := a + b
	if (a > 0 && b > 0 && c < 0) || (a < 0 && b < 0 && c >= 0) {
		fatal("overflow")
	}
	return c
}

func mulchk(a, b int64) int64 {
	if a == 0 || b == 0 {
		return 0
	}
	if (a == -1 && b == math.MinInt64) || (b == -1 && a == math.MinInt64) {
		fatal("overflow")
	}
	c := a * b
	if c/b != a {
		fatal("overflow")
	}
	return c
}

// addto adds scale * q^shift * b to a (growing a as needed) and returns it.
func addto(a Poly, b Poly, shift int, scale int64) Poly {
	if len(a) < len(b)+shift {
		a = append(a, make(Poly, len(b)+shift-len(a))...)
	}
	for i := range b {
		t := mulchk(b[i], scale)
		a[i+shift] = addchk(a[i+shift], t)
	}
	return a
}

func pjson(p Poly) string {
	var sb strings.Builder
	sb.WriteString("[")
	for i, v := range p {
		if i > 0 {
			sb.WriteString(",")
		}
		sb.WriteString(strconv.FormatInt(v, 10))
	}
	sb.WriteString("]")
	return sb.String()
}

func wjson(w []int) string {
	var sb strings.Builder
	sb.WriteString("[")
	for i, v := range w {
		if i > 0 {
			sb.WriteString(",")
		}
		sb.WriteString(strconv.Itoa(v))
	}
	sb.WriteString("]")
	return sb.String()
}

func wstr(w []int) string {
	var sb strings.Builder
	for _, x := range w {
		sb.WriteString(strconv.Itoa(x))
	}
	return sb.String()
}

func masklist(m int) []int {
	v := []int{}
	for i := 0; i < N; i++ {
		if m>>i&1 == 1 {
			v = append(v, i)
		}
	}
	return v
}

func ctz(x int) int { return bits.TrailingZeros(uint(x)) }

// Quot is the quotient W^I of minimal left coset representatives inside the
// lower interval of w w_I.
type Quot struct {
	I      int
	wI     Mat
	wIword []int
	lwI    int
	el     []Mat
	len    []int
	Lq     []int
	left   [][9]int // index of s x in the quotient, -1 if sx notin W^I (sx = x t), -2 if not in the ideal
	idx    map[Mat]int
	bucket [][]int // elements by left descent mask, sorted by length
}

func (q *Quot) inWI(x Mat) bool {
	for t := 0; t < N; t++ {
		if q.I>>t&1 == 1 {
			if negcol(&x, t) {
				return false
			}
		}
	}
	return true
}

// toFull returns x w_I for a quotient representative x.
func (q *Quot) toFull(x Mat) Mat {
	for _, t := range q.wIword {
		x = rmul(x, t)
	}
	return x
}

func (q *Quot) build(w Mat, Imask int) {
	q.I = Imask
	q.wI = identity()
	for {
		moved := false
		for t := 0; t < N; t++ {
			if q.I>>t&1 == 1 {
				if !negcol(&q.wI, t) {
					q.wI = rmul(q.wI, t)
					moved = true
				}
			}
		}
		if !moved {
			break
		}
	}
	q.wIword = rword(q.wI)
	q.lwI = len(q.wIword)
	if q.lwI != bits.OnesCount(uint(q.I)) {
		fatal("I is not a commuting set")
	}
	wp := q.toFull(w)
	if !q.inWI(wp) {
		fatal("I is not contained in R(w)")
	}
	wd := rword(wp)
	seen := map[Mat]int{}
	cur := []Mat{}
	e := identity()
	seen[e] = 0
	cur = append(cur, e)
	for k := len(wd) - 1; k >= 0; k-- {
		s := wd[k]
		sz := len(cur)
		for i := 0; i < sz; i++ {
			y := lmul(cur[i], s)
			if !q.inWI(y) {
				continue
			}
			if _, ok := seen[y]; ok {
				continue
			}
			seen[y] = len(cur)
			cur = append(cur, y)
		}
	}
	type pr struct{ l, i int }
	order := make([]pr, 0, len(cur))
	for i := range cur {
		order = append(order, pr{len(rword(cur[i])), i})
	}
	sort.Slice(order, func(a, b int) bool {
		if order[a].l != order[b].l {
			return order[a].l < order[b].l
		}
		return order[a].i < order[b].i
	})
	q.el = make([]Mat, len(cur))
	q.len = make([]int, len(cur))
	q.idx = make(map[Mat]int, len(cur))
	for i := range order {
		q.el[i] = cur[order[i].i]
		q.len[i] = order[i].l
		q.idx[q.el[i]] = i
	}
	q.left = make([][9]int, len(q.el))
	q.Lq = make([]int, len(q.el))
	for i := range q.el {
		for s := 0; s < N; s++ {
			y := lmul(q.el[i], s)
			if !q.inWI(y) {
				q.left[i][s] = -1
				q.Lq[i] |= 1 << s
				continue
			}
			j, ok := q.idx[y]
			if !ok {
				q.left[i][s] = -2
				continue
			}
			q.left[i][s] = j
			if q.len[j] < q.len[i] {
				q.Lq[i] |= 1 << s
			} else if q.len[j] != q.len[i]+1 {
				fatal("length mismatch")
			}
		}
	}
	q.bucket = make([][]int, 1<<N)
	for i := range q.el {
		q.bucket[q.Lq[i]] = append(q.bucket[q.Lq[i]], i)
	}
}

type lpair struct{ x, id int }

type corrTerm struct {
	z     int
	coef  int64
	shift int
}

// KLQ computes the quotient polynomials Pq(x, y).
type KLQ struct {
	Q        *Quot
	polys    []Poly
	pid      map[string]int
	lst      [][]lpair
	done     []bool
	pairs    int
	computed int
}

func newKLQ(q *Quot) *KLQ {
	return &KLQ{Q: q, pid: map[string]int{}, lst: make([][]lpair, len(q.el)), done: make([]bool, len(q.el))}
}

func (k *KLQ) id(p Poly) int {
	buf := make([]byte, 0, len(p)*8)
	for _, v := range p {
		u := uint64(v)
		for b := 0; b < 8; b++ {
			buf = append(buf, byte(u>>(8*b)))
		}
	}
	key := string(buf)
	if i, ok := k.pid[key]; ok {
		return i
	}
	i := len(k.polys)
	k.polys = append(k.polys, append(Poly(nil), p...))
	k.pid[key] = i
	return i
}

// get returns Pq(x, z); a nil slice stands for the zero polynomial.
func (k *KLQ) get(x, z int) Poly {
	Q := k.Q
	if Q.len[x] > Q.len[z] {
		return nil
	}
	Lz := Q.Lq[z]
	for {
		d := Lz &^ Q.Lq[x]
		if d == 0 {
			break
		}
		t := ctz(d)
		nx := Q.left[x][t]
		if nx < 0 {
			if nx == -2 {
				return nil
			}
			fatal("impossible reduction")
		}
		x = nx
		if Q.len[x] > Q.len[z] {
			return nil
		}
	}
	v := k.lst[z]
	i := sort.Search(len(v), func(i int) bool { return v[i].x >= x })
	if i == len(v) || v[i].x != x {
		return nil
	}
	return k.polys[v[i].id]
}

func (k *KLQ) compute(y int) {
	Q := k.Q
	if k.done[y] {
		return
	}
	if Q.len[y] == 0 {
		k.lst[y] = []lpair{{y, k.id(Poly{1})}}
		k.done[y] = true
		return
	}
	s := -1
	for c := 0; c < N; c++ {
		if (Q.Lq[y]>>c&1) == 1 && Q.left[y][c] >= 0 {
			if s < 0 {
				s = c
			}
			if k.done[Q.left[y][c]] {
				s = c
				break
			}
		}
	}
	if s < 0 {
		fatal("no quotient descent")
	}
	u := Q.left[y][s]
	k.compute(u)
	corr := []corrTerm{}
	for _, pr := range k.lst[u] {
		z := pr.x
		if z == u || (Q.Lq[z]>>s&1) == 0 {
			continue
		}
		gap := Q.len[u] - Q.len[z]
		if gap%2 == 0 {
			continue
		}
		p := k.polys[pr.id]
		top := (gap - 1) / 2
		if len(p) > top+1 {
			fatal("degree violation")
		}
		if len(p) == top+1 && p[top] != 0 {
			corr = append(corr, corrTerm{z, p[top], (Q.len[y] - Q.len[z]) / 2})
		}
	}
	for t := 0; t < N; t++ {
		if (Q.Lq[u]>>t&1) == 1 && Q.left[u][t] >= 0 {
			z := Q.left[u][t]
			if (Q.Lq[z] >> s & 1) == 1 {
				corr = append(corr, corrTerm{z, 1, 1})
			}
		}
	}
	for _, c := range corr {
		k.compute(c.z)
	}
	out := []lpair{}
	Ly := Q.Lq[y]
	full := (1 << N) - 1
	freebits := full &^ Ly
	for sub := freebits; ; sub = (sub - 1) & freebits {
		D := Ly | sub
		b := Q.bucket[D]
		for _, x := range b {
			if Q.len[x] >= Q.len[y] {
				break
			}
			var p Poly
			sx := Q.left[x][s]
			if sx == -1 {
				a := k.get(x, u)
				p = addto(p, a, 0, 1)
				p = addto(p, a, 1, 1)
			} else {
				if sx < 0 || Q.len[sx] >= Q.len[x] {
					fatal("candidate without s-descent")
				}
				p = addto(p, k.get(sx, u), 0, 1)
				p = addto(p, k.get(x, u), 1, 1)
			}
			for _, c := range corr {
				z := c.z
				if Q.len[x] > Q.len[z] {
					continue
				}
				pz := k.get(x, z)
				if len(pz) != 0 {
					p = addto(p, pz, c.shift, -c.coef)
				}
			}
			p = trim(p)
			if len(p) == 0 {
				continue
			}
			gap := Q.len[y] - Q.len[x]
			if len(p)-1 > (gap-1)/2 {
				fatal("degree bound violated")
			}
			for _, c := range p {
				if c < 0 {
					fatal("negative coefficient")
				}
			}
			if p[0] != 1 {
				fatal("constant term != 1")
			}
			out = append(out, lpair{x, k.id(p)})
		}
		if sub == 0 {
			break
		}
	}
	out = append(out, lpair{y, k.id(Poly{1})})
	sort.Slice(out, func(a, b int) bool {
		if out[a].x != out[b].x {
			return out[a].x < out[b].x
		}
		return out[a].id < out[b].id
	})
	k.lst[y] = out
	k.done[y] = true
	k.pairs += len(out)
	k.computed++
}

// Ord is the naive ordinary Kazhdan-Lusztig recursion on the lower interval of w, for validation.
type Ord struct {
	el   []Mat
	len  []int
	idx  map[Mat]int
	left [][9]int
	L    []int
	memo map[[2]int]Poly
}

func (o *Ord) build(w Mat) {
	wd := rword(w)
	seen := map[Mat]int{}
	cur := []Mat{}
	e := identity()
	seen[e] = 0
	cur = append(cur, e)
	for _, s := range wd {
		sz := len(cur)
		for i := 0; i < sz; i++ {
			y := rmul(cur[i], s)
			if _, ok := seen[y]; ok {
				continue
			}
			seen[y] = len(cur)
			cur = append(cur, y)
		}
	}
	type pr struct{ l, i int }
	order := []pr{}
	for i := range cur {
		order = append(order, pr{len(rword(cur[i])), i})
	}
	sort.Slice(order, func(a, b int) bool {
		if order[a].l != order[b].l {
			return order[a].l < order[b].l
		}
		return order[a].i < order[b].i
	})
	o.el = make([]Mat, len(cur))
	o.len = make([]int, len(cur))
	o.idx = map[Mat]int{}
	for i := range order {
		o.el[i] = cur[order[i].i]
		o.len[i] = order[i].l
		o.idx[o.el[i]] = i
	}
	o.left = make([][9]int, len(o.el))
	o.L = make([]int, len(o.el))
	o.memo = map[[2]int]Poly{}
	for i := range o.el {
		for s := 0; s < N; s++ {
			y := lmul(o.el[i], s)
			j, ok := o.idx[y]
			if !ok {
				o.left[i][s] = -2
				continue
			}
			o.left[i][s] = j
			if o.len[j] < o.len[i] {
				o.L[i] |= 1 << s
			}
		}
	}
}

func (o *Ord) leq(x, y int) bool {
	for x != y {
		if o.len[x] >= o.len[y] {
			return false
		}
		s := ctz(o.L[y])
		if o.L[x]>>s&1 == 1 {
			x = o.left[x][s]
		}
		y = o.left[y][s]
	}
	return true
}

func (o *Ord) P(x, y int) Poly {
	if x == y {
		return Poly{1}
	}
	if !o.leq(x, y) {
		return Poly{}
	}
	key := [2]int{x, y}
	if v, ok := o.memo[key]; ok {
		return v
	}
	s := ctz(o.L[y])
	u := o.left[y][s]
	var p Poly
	c := o.L[x] >> s & 1
	sx := o.left[x][s]
	if sx == -2 {
		fatal("ord: sx missing")
	}
	p = addto(p, o.P(sx, u), 1-c, 1)
	p = addto(p, o.P(x, u), c, 1)
	for z := 0; z < len(o.el); z++ {
		if o.len[z] >= o.len[u] {
			break
		}
		if o.L[z]>>s&1 == 0 {
			continue
		}
		gap := o.len[u] - o.len[z]
		if gap%2 == 0 || !o.leq(x, z) {
			continue
		}
		pz := o.P(z, u)
		top := (gap - 1) / 2
		if len(pz) == top+1 && pz[top] != 0 {
			p = addto(p, o.P(x, z), (o.len[y]-o.len[z])/2, -pz[top])
		}
	}
	p = trim(p)
	o.memo[key] = p
	return p
}

func setup(typ string) {
	ADJ = nil
	var add func(a, b int)
	add = func(a, b int) {
		ADJ[a] = append(ADJ[a], b)
		ADJ[b] = append(ADJ[b], a)
	}
	n, err := strconv.Atoi(typ[1:])
	if err != nil {
		fmt.Fprintln(os.Stderr, "unknown type")
		os.Exit(2)
	}
	N = n
	switch typ[0] {
	case 'E':
		ADJ = make([][]int, N)
		for i := 0; i+1 < N-1; i++ {
			add(i, i+1)
		}
		add(2, N-1)
	case 'D':
		ADJ = make([][]int, N)
		for i := 0; i+1 < N-1; i++ {
			add(i, i+1)
		}
		add(N-3, N-1)
	case 'A':
		ADJ = make([][]int, N)
		for i := 0; i+1 < N; i++ {
			add(i, i+1)
		}
	default:
		fmt.Fprintln(os.Stderr, "unknown type")
		os.Exit(2)
	}
	if N > 9 {
		fmt.Fprintln(os.Stderr, "rank at most 9")
		os.Exit(2)
	}
}

// mt19937 is the 32-bit Mersenne twister (std::mt19937), so that the random
// validation samples coincide with those of the original program.
type mt19937 struct {
	s [624]uint32
	i int
}

func newMT(seed uint32) *mt19937 {
	m := &mt19937{}
	m.s[0] = seed
	for i := 1; i < 624; i++ {
		m.s[i] = 1812433253*(m.s[i-1]^(m.s[i-1]>>30)) + uint32(i)
	}
	m.i = 624
	return m
}

func (m *mt19937) next() uint32 {
	if m.i >= 624 {
		for k := 0; k < 624; k++ {
			y := (m.s[k] & 0x80000000) | (m.s[(k+1)%624] & 0x7fffffff)
			v := m.s[(k+397)%624] ^ (y >> 1)
			if y&1 == 1 {
				v ^= 0x9908b0df
			}
			m.s[k] = v
		}
		m.i = 0
	}
	y := m.s[m.i]
	m.i++
	y ^= y >> 11
	y ^= (y << 7) & 0x9d2c5680
	y ^= (y << 15) & 0xefc60000
	y ^= y >> 18
	return y
}

// validate compares polynomials for every quotient pair, using random
// elements w and commuting I subset R(w), with the naive recursion.
// It returns the number of pairs compared.
func validate(typ string, seed uint32, count, maxlen int) int64 {
	setup(typ)
	rng := newMT(seed)
	var pairs int64
	tested := 0
	for it := 0; it < count; it++ {
		L := 1 + int(rng.next()%uint32(maxlen))
		w := identity()
		for i := 0; i < L; i++ {
			s := int(rng.next() % uint32(N))
			if !negcol(&w, s) {
				w = rmul(w, s)
			}
		}
		R := Rmask(w)
		if R == 0 {
			continue
		}
		I := 0
		for s := 0; s < N; s++ {
			if R>>s&1 == 1 && rng.next()%2 == 1 {
				adj := false
				for _, t := range ADJ[s] {
					if I>>t&1 == 1 {
						adj = true
					}
				}
				if !adj {
					I |= 1 << s
				}
			}
		}
		if I == 0 {
			I = 1 << ctz(R)
		}
		var o Ord
		o.build(w)
		if len(o.el) > 6000 {
			continue
		}
		var q Quot
		q.build(w, I)
		k := newKLQ(&q)
		wp := q.toFull(w)
		yq := q.idx[wp]
		k.compute(yq)
		ow := o.idx[w]
		for xi := range q.el {
			xm := q.toFull(q.el[xi])
			var ref Poly
			if f, ok := o.idx[xm]; ok {
				ref = o.P(f, ow)
			}
			if !slices.Equal(ref, k.get(xi, yq)) {
				fmt.Fprintf(os.Stderr, "validation mismatch in %s\n", typ)
				os.Exit(1)
			}
			pairs++
		}
		tested++
	}
	fmt.Fprintf(os.Stderr, "validated %s: %d elements, polynomials match the naive recursion on %d quotient pairs\n", typ, tested, pairs)
	return pairs
}

// parseWord reads a word for the -pairs debugging mode: a string of generator
// labels 0..N-1, or "e" for the identity.
func parseWord(w string) []int {
	if w == "e" {
		return nil
	}
	out := make([]int, 0, len(w))
	for _, c := range []byte(w) {
		if c < '0' || int(c-'0') >= N {
			fmt.Fprintln(os.Stderr, "bad generator label in word", w)
			os.Exit(2)
		}
		out = append(out, int(c-'0'))
	}
	return out
}

// pairTop caches the quotient and the naive engine for one top element w.
type pairTop struct {
	q  *Quot
	k  *KLQ
	yq int
	o  *Ord // nil when the lower ideal of w exceeds 6000 elements
	ow int
}

// runPairs implements the hidden -mode pairs (see the file header): for
// every line "xword wword" it prints P_{x,w} from the quotient engine with
// I = R(w) and, for small ideals, from the naive recursion.  It is used only
// for differential testing against other implementations.
func runPairs(typ, path string) {
	setup(typ)
	f, err := os.Open(path)
	if err != nil {
		fatal("cannot open pairs file %s", path)
	}
	defer f.Close()
	tops := map[string]*pairTop{}
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		fields := strings.Fields(sc.Text())
		if len(fields) == 0 || strings.HasPrefix(fields[0], "#") {
			continue
		}
		if len(fields) < 2 {
			fatal("pairs line needs two words")
		}
		t, ok := tops[fields[1]]
		if !ok {
			w := fromword(parseWord(fields[1]))
			// I = greedy maximal commuting subset of R(w) (ascending labels).
			R, I := Rmask(w), 0
			for s := 0; s < N; s++ {
				if R>>s&1 == 1 {
					adjacent := false
					for _, u := range ADJ[s] {
						if I>>u&1 == 1 {
							adjacent = true
						}
					}
					if !adjacent {
						I |= 1 << s
					}
				}
			}
			t = &pairTop{q: &Quot{}}
			t.q.build(w, I)
			t.k = newKLQ(t.q)
			t.yq = t.q.idx[t.q.toFull(w)]
			t.k.compute(t.yq)
			var o Ord
			o.build(w)
			if len(o.el) <= 6000 {
				t.o = &o
				t.ow = o.idx[w]
			}
			tops[fields[1]] = t
		}
		x := fromword(parseWord(fields[0]))
		// Minimal coset representative of x W_I: strip right descents in I.
		xq := x
		for moved := true; moved; {
			moved = false
			for s := 0; s < N; s++ {
				if t.q.I>>s&1 == 1 && negcol(&xq, s) {
					xq = rmul(xq, s)
					moved = true
				}
			}
		}
		var pq Poly
		if xi, ok := t.q.idx[xq]; ok {
			pq = t.k.get(xi, t.yq)
		}
		ord := "-"
		if t.o != nil {
			var po Poly
			if xi, ok := t.o.idx[x]; ok {
				po = t.o.P(xi, t.ow)
			}
			ord = pjson(po)
		}
		fmt.Printf("%s %s %s %s\n", fields[0], fields[1], pjson(pq), ord)
	}
	if err := sc.Err(); err != nil {
		fatal("reading pairs file: %v", err)
	}
}

type target struct {
	name string
	ext  []int
}

var failures = 0

func expect(cond bool, what string) {
	if !cond {
		fmt.Fprintf(os.Stderr, "EXPECTATION FAILED: %s\n", what)
		failures++
	}
}

func intsEqual(a, b []int) bool { return slices.Equal(a, b) }

type row struct {
	xs  string
	xl  int
	gap int
	mu  int64
	p   Poly
}

func certifyOne(tg target, beta, baseword []int) string {
	setup("E9")
	w := reflection(beta)
	expect(w == fromword(baseword), "r_beta equals the manuscript's reduced word")
	for _, s := range tg.ext {
		if negcol(&w, s) {
			fatal("extension letter is a descent")
		}
		w = rmul(w, s)
	}
	wd := rword(w)
	lw := len(wd)
	R, L := Rmask(w), Lmask(w)
	I := R
	winv := inverse(w)
	invol := winv == w
	term := true
	for _, v := range []Mat{w, winv} {
		for s := 0; s < N; s++ {
			if negcol(&v, s) {
				vs := rmul(v, s)
				for _, t := range ADJ[s] {
					if negcol(&vs, t) {
						term = false
					}
				}
			}
		}
	}
	var q Quot
	q.build(w, I)
	qsize := len(q.el)
	popI := bits.OnesCount(uint(I))
	ideal := int64(qsize) * (int64(1) << popI)
	// Eligible fully commutative lower endpoints: x <= w, L(w) subset L(x), R(w) subset R(x), x FC.
	elig := []int{}
	for xi := 0; xi < qsize; xi++ {
		xm := q.toFull(q.el[xi])
		if (Lmask(xm)&L) != L || (Rmask(xm)&R) != R {
			continue
		}
		if !isFC(xm) {
			continue
		}
		elig = append(elig, xi)
	}
	wp := q.toFull(w)
	yq := q.idx[wp]
	k := newKLQ(&q)
	k.compute(yq)
	var js strings.Builder
	bs := func(b bool) string {
		if b {
			return "true"
		}
		return "false"
	}
	fmt.Fprintf(&js, "    {\"name\": \"%s\", \"word\": %s, \"word_string\": \"%s\", \"length\": %d", tg.name, wjson(wd), wstr(wd), lw)
	fmt.Fprintf(&js, ", \"L\": %s, \"R\": %s, \"involution\": %s, \"terminal\": %s", wjson(masklist(L)), wjson(masklist(R)), bs(invol), bs(term))
	fmt.Fprintf(&js, ", \"I\": %s, \"W_I_order\": %d, \"quotient_size\": %d, \"lower_ideal_size\": %d", wjson(masklist(I)), 1<<popI, qsize, ideal)
	fmt.Fprintf(&js, ", \"quotient_elements_computed\": %d, \"stored_pairs\": %d, \"distinct_polynomials\": %d", k.computed, k.pairs, len(k.polys))
	js.WriteString(",\n     \"eligible_fully_commutative_bottoms\": [")
	rows := []row{}
	for _, xi := range elig {
		xm := q.toFull(q.el[xi])
		xw := rword(xm)
		p := k.get(xi, yq)
		gap := lw - len(xw)
		var mu int64
		if gap%2 == 1 {
			top := (gap - 1) / 2
			if len(p) == top+1 {
				mu = p[top]
			}
		}
		rows = append(rows, row{wstr(xw), len(xw), gap, mu, p})
	}
	sort.Slice(rows, func(a, b int) bool { return rows[a].xs < rows[b].xs })
	for i, r := range rows {
		if i > 0 {
			js.WriteString(",\n       ")
		} else {
			js.WriteString("\n       ")
		}
		fmt.Fprintf(&js, "{\"x\": \"%s\", \"length\": %d, \"gap\": %d, \"P_ascending\": %s, \"degree\": %d, \"mu\": %d}", r.xs, r.xl, r.gap, pjson(r.p), len(r.p)-1, r.mu)
	}
	js.WriteString("],\n")
	// Histogram of nonzero mu over all extremal lower endpoints in the quotient.
	var maxmu int64
	hist := map[int64]int{}
	for _, pr := range k.lst[yq] {
		x := pr.x
		if x == yq {
			continue
		}
		gap := q.len[yq] - q.len[x]
		if gap%2 == 0 {
			continue
		}
		p := k.polys[pr.id]
		top := (gap - 1) / 2
		if len(p) == top+1 && p[top] != 0 {
			hist[p[top]]++
			if p[top] > maxmu {
				maxmu = p[top]
			}
		}
	}
	js.WriteString("     \"nonzero_mu_histogram_over_quotient_representatives\": {")
	keys := make([]int64, 0, len(hist))
	for h := range hist {
		keys = append(keys, h)
	}
	sort.Slice(keys, func(a, b int) bool { return keys[a] < keys[b] })
	for i, h := range keys {
		if i > 0 {
			js.WriteString(", ")
		}
		fmt.Fprintf(&js, "\"%d\": %d", h, hist[h])
	}
	js.WriteString("}")
	fmt.Fprintf(&js, ", \"max_mu_over_quotient_representatives\": %d", maxmu)
	checked, bad, missing := 0, 0, 0
	if invol {
		for _, pr := range k.lst[yq] {
			X := q.toFull(q.el[pr.x])
			Xi := q.toFull(inverse(X))
			it, ok := q.idx[Xi]
			if !ok {
				missing++
				continue
			}
			if !slices.Equal(k.get(it, yq), k.polys[pr.id]) {
				bad++
			}
			checked++
		}
	}
	fmt.Fprintf(&js, ", \"inversion_symmetry\": {\"pairs_checked\": %d, \"mismatches\": %d, \"inverse_outside_quotient\": %d}", checked, bad, missing)
	var maxc int64
	maxdeg := 0
	for _, p := range k.polys {
		for _, c := range p {
			if c > maxc {
				maxc = c
			}
		}
		if len(p) > maxdeg {
			maxdeg = len(p)
		}
	}
	md := 0
	if maxdeg > 0 {
		md = maxdeg - 1
	}
	fmt.Fprintf(&js, ", \"max_coefficient\": %d, \"max_degree\": %d}", maxc, md)
	expect(bad == 0, "P(x,w) = P(x^-1,w) for the involution w")
	// Expected values (reviewer computation of 7 October 2026).
	expectRow := func(xs string, gap int, P Poly, mu int64) {
		found := false
		for _, r := range rows {
			if r.xs == xs {
				found = true
				expect(r.gap == gap, "gap of "+xs)
				expect(slices.Equal(r.p, P), "polynomial of "+xs)
				expect(r.mu == mu, "mu of "+xs)
			}
		}
		expect(found, "eligible bottom "+xs+" present")
	}
	switch tg.name {
	case "w33":
		expect(lw == 33 && invol && term, "length 33 terminal involution")
		expect(intsEqual(masklist(R), []int{3, 5, 7, 8}) && intsEqual(masklist(L), []int{3, 5, 7, 8}), "L = R = {3,5,7,8}")
		expect(qsize == 364156 && ideal == 5826496, "quotient 364,156 cosets, ideal 5,826,496 elements")
		expect(len(rows) == 5, "exactly five eligible FC bottoms")
		expectRow("8753", 29, Poly{1, 23, 236, 1391, 5298, 13861, 25666, 33996, 31954, 20820, 8985, 2380, 346, 23}, 0)
		expectRow("87531", 28, Poly{1, 23, 236, 1389, 5262, 13609, 24701, 31740, 28612, 17622, 7052, 1654, 194, 9}, 0)
		expectRow("87530", 28, Poly{1, 23, 235, 1370, 5128, 13088, 23433, 29667, 26255, 15772, 6082, 1350, 144, 5}, 0)
		expectRow("875301", 27, Poly{1, 23, 234, 1348, 4931, 12136, 20641, 24398, 19796, 10646, 3551, 621, 40, 1}, 1)
		expectRow("875310", 27, Poly{1, 23, 234, 1348, 4931, 12136, 20641, 24398, 19796, 10646, 3551, 621, 40, 1}, 1)
		distinct := map[int]bool{}
		for _, s := range wd {
			distinct[s] = true
		}
		expect(len(distinct) == 9, "full support")
	case "w34_s0":
		expect(lw == 34 && invol && term, "length 34 terminal involution")
		expect(intsEqual(masklist(R), []int{0, 3, 5, 7, 8}) && intsEqual(masklist(L), []int{0, 3, 5, 7, 8}), "L = R = {0,3,5,7,8}")
		expect(qsize == 216990 && ideal == 6943680, "quotient 216,990 cosets, ideal 6,943,680 elements")
		expect(len(rows) == 1, "unique eligible FC bottom")
		expectRow("87530", 29, Poly{1, 18, 146, 712, 2342, 5456, 9200, 11264, 9947, 6220, 2649, 721, 116, 10, 1}, 1)
	case "w34_s1":
		expect(lw == 34 && invol && term, "length 34 terminal involution")
		expect(intsEqual(masklist(R), []int{1, 3, 5, 7, 8}) && intsEqual(masklist(L), []int{1, 3, 5, 7, 8}), "L = R = {1,3,5,7,8}")
		expect(qsize == 207866 && ideal == 6651712, "quotient 207,866 cosets, ideal 6,651,712 elements")
		expect(len(rows) == 1, "unique eligible FC bottom")
		expectRow("87531", 29, Poly{1, 18, 140, 639, 1957, 4280, 6858, 8114, 7084, 4497, 1996, 579, 100, 10, 1}, 1)
	}
	fmt.Fprintf(os.Stderr, "%s: length %d, quotient %d, ideal %d, %d eligible FC bottoms\n", tg.name, lw, qsize, ideal, len(rows))
	FCMEMO = map[Mat]bool{}
	return js.String()
}

func main() {
	mode := flag.String("mode", "", "certify or validate")
	out := flag.String("out", "", "output JSON path (certify)")
	typ := flag.String("type", "", "Coxeter type, e.g. E6 (validate)")
	seed := flag.Uint("seed", 0, "random seed (validate)")
	count := flag.Int("count", 0, "number of random elements (validate)")
	maxlen := flag.Int("maxlen", 0, "maximal random word length (validate)")
	pairs := flag.String("pairs", "", "debugging: file of \"xword wword\" lines (mode pairs)")
	flag.Parse()
	switch *mode {
	case "pairs":
		if *typ == "" || *pairs == "" {
			fmt.Fprintln(os.Stderr, "usage: e9-quotient-kl -mode pairs -type T -pairs FILE")
			os.Exit(2)
		}
		runPairs(*typ, *pairs)
		return
	case "validate":
		if *typ == "" || *count <= 0 || *maxlen <= 0 {
			fmt.Fprintln(os.Stderr, "usage: e9-quotient-kl -mode validate -type T -seed S -count C -maxlen M")
			os.Exit(2)
		}
		p := validate(*typ, uint32(*seed), *count, *maxlen)
		fmt.Printf("validated: %d pairs\n", p)
		return
	case "certify":
		if *out == "" {
			fmt.Fprintln(os.Stderr, "usage: e9-quotient-kl -mode certify -out <output.json>")
			os.Exit(2)
		}
	default:
		fmt.Fprintln(os.Stderr, "usage: e9-quotient-kl -mode certify -out <output.json> | -mode validate -type T -seed S -count C -maxlen M")
		os.Exit(2)
	}
	va := validate("A4", 1, 200, 10)
	vd := validate("D4", 2, 200, 12)
	ve := validate("E6", 3, 120, 14)
	beta := []int{1, 2, 3, 3, 2, 2, 1, 1, 2}
	baseword := []int{8, 7, 5, 6, 3, 4, 5, 2, 3, 4, 1, 2, 3, 0, 1, 2, 8, 2, 3, 4, 1, 2, 3, 0, 1, 2, 8, 5, 6, 7, 4, 5, 3}
	targets := []target{{"w33", nil}, {"w34_s0", []int{0}}, {"w34_s1", []int{1}}}
	parts := []string{}
	for _, t := range targets {
		parts = append(parts, certifyOne(t, beta, baseword))
	}
	var sb strings.Builder
	sb.WriteString("{\n  \"group\": \"affine E8 = E9, chain 0-1-...-7, node 8 attached to node 2\",\n")
	sb.WriteString("  \"method\": \"Kazhdan-Lusztig recursion in the parabolic quotient W/W_I, I = R(w); exact integers with overflow checks\",\n")
	sb.WriteString("  \"beta\": [1,2,3,3,2,2,1,1,2],\n  \"base_element\": \"r_beta, the manuscript's length-33 reflection\",\n")
	fmt.Fprintf(&sb, "  \"engine_validation\": {\"A4_pairs\": %d, \"D4_pairs\": %d, \"E6_pairs\": %d, \"all_equal_to_naive_recursion\": true},\n", va, vd, ve)
	sb.WriteString("  \"elements\": [\n")
	for i, p := range parts {
		sb.WriteString(p)
		if i+1 < len(parts) {
			sb.WriteString(",\n")
		} else {
			sb.WriteString("\n")
		}
	}
	status := "passed"
	if failures > 0 {
		status = "FAILED"
	}
	fmt.Fprintf(&sb, "  ],\n  \"all_eligible_mu_values_in\": [0, 1],\n  \"status\": \"%s\"\n}\n", status)
	if err := os.WriteFile(*out, []byte(sb.String()), 0o644); err != nil {
		fmt.Fprintln(os.Stderr, "cannot open output")
		os.Exit(2)
	}
	fmt.Printf("E9 odd-gap certificate: %s\n", status)
	if failures > 0 {
		os.Exit(1)
	}
}
