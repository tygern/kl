package main

import (
	"math/big"
	"sort"
)

// vec is a vector in the simple-root basis; mat is a matrix stored as its
// columns (mat[j] is the image of alpha_j), as in the original.
type vec = []int
type mat = []vec

func keyVec(v vec) string {
	b := make([]byte, 0, 2*len(v))
	for _, c := range v {
		if c < -32768 || c > 32767 {
			panic("coefficient out of int16 range")
		}
		b = append(b, byte(uint16(int16(c))>>8), byte(uint16(int16(c))))
	}
	return string(b)
}

func keyMat(M mat) string {
	b := make([]byte, 0, 2*len(M)*len(M))
	for _, c := range M {
		b = append(b, keyVec(c)...)
	}
	return string(b)
}

func equalVec(a, b vec) bool {
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

func equalMat(a, b mat) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !equalVec(a[i], b[i]) {
			return false
		}
	}
	return true
}

func cmpVec(a, b vec) int {
	for i := range a {
		if a[i] != b[i] {
			if a[i] < b[i] {
				return -1
			}
			return 1
		}
	}
	return 0
}

// bits is a subset of the (at most 128) positive roots.
type bits [2]uint64

func (b *bits) set(i int)           { b[i/64] |= 1 << uint(i%64) }
func (b bits) subsetOf(c bits) bool { return b[0]&^c[0] == 0 && b[1]&^c[1] == 0 }

// rootSystem is a simply laced root system given by Dynkin edges.
type rootSystem struct {
	n        int
	edges    [][2]int
	adj      [][]int
	adjSet   []map[int]bool
	identity mat
	roots    []vec
}

func newRootSystem(n int, edges [][2]int) *rootSystem {
	rs := &rootSystem{n: n, edges: edges}
	rs.adj = make([][]int, n)
	rs.adjSet = make([]map[int]bool, n)
	for i := range rs.adjSet {
		rs.adjSet[i] = map[int]bool{}
	}
	for _, e := range edges {
		a, b := e[0], e[1]
		rs.adjSet[a][b] = true
		rs.adjSet[b][a] = true
	}
	for i := 0; i < n; i++ {
		for t := range rs.adjSet[i] {
			rs.adj[i] = append(rs.adj[i], t)
		}
		sort.Ints(rs.adj[i])
	}
	rs.identity = rs.unit()
	return rs
}

func (rs *rootSystem) unit() mat {
	M := make(mat, rs.n)
	for j := range M {
		M[j] = make(vec, rs.n)
		M[j][j] = 1
	}
	return M
}

func (rs *rootSystem) pairing(u, v vec) int {
	total := 0
	for i := 0; i < rs.n; i++ {
		sv := 0
		for _, t := range rs.adj[i] {
			sv += v[t]
		}
		total += 2*u[i]*v[i] - u[i]*sv
	}
	return total
}

func (rs *rootSystem) simpleReflect(s int, v vec) vec {
	c := 2 * v[s]
	for _, t := range rs.adj[s] {
		c -= v[t]
	}
	w := append(vec(nil), v...)
	w[s] -= c
	return w
}

func isNegative(v vec) bool {
	any := false
	for _, c := range v {
		if c > 0 {
			return false
		}
		if c < 0 {
			any = true
		}
	}
	return any
}

func isPositive(v vec) bool {
	any := false
	for _, c := range v {
		if c < 0 {
			return false
		}
		if c > 0 {
			any = true
		}
	}
	return any
}

// positiveRoots returns the positive roots sorted by (height, coefficients).
func (rs *rootSystem) positiveRoots() []vec {
	if rs.roots != nil {
		return rs.roots
	}
	seen := map[string]vec{}
	var frontier []vec
	for j := 0; j < rs.n; j++ {
		e := make(vec, rs.n)
		e[j] = 1
		seen[keyVec(e)] = e
		frontier = append(frontier, e)
	}
	for len(frontier) > 0 {
		var next []vec
		for _, r := range frontier {
			for s := 0; s < rs.n; s++ {
				q := rs.simpleReflect(s, r)
				ok := true
				for _, c := range q {
					if c < 0 {
						ok = false
					}
				}
				if !ok {
					continue
				}
				k := keyVec(q)
				if _, dup := seen[k]; !dup {
					seen[k] = q
					next = append(next, q)
				}
			}
		}
		frontier = next
	}
	roots := make([]vec, 0, len(seen))
	for _, r := range seen {
		roots = append(roots, r)
	}
	height := func(r vec) int {
		h := 0
		for _, c := range r {
			h += c
		}
		return h
	}
	sort.Slice(roots, func(a, b int) bool {
		ha, hb := height(roots[a]), height(roots[b])
		if ha != hb {
			return ha < hb
		}
		return cmpVec(roots[a], roots[b]) < 0
	})
	if len(roots) > 128 {
		panic("too many positive roots for bitset")
	}
	rs.roots = roots
	return roots
}

func (rs *rootSystem) rightMult(M mat, s int) mat {
	cols := make(mat, rs.n)
	for j := range M {
		cols[j] = append(vec(nil), M[j]...)
	}
	old := M[s]
	neg := make(vec, rs.n)
	for i, v := range old {
		neg[i] = -v
	}
	cols[s] = neg
	for _, t := range rs.adj[s] {
		for i := range cols[t] {
			cols[t][i] += old[i]
		}
	}
	return cols
}

func (rs *rootSystem) leftMult(s int, M mat) mat {
	out := make(mat, rs.n)
	for j, c := range M {
		out[j] = rs.simpleReflect(s, c)
	}
	return out
}

func (rs *rootSystem) wordMatrix(word []int) mat {
	M := rs.identity
	for _, s := range word {
		M = rs.rightMult(M, s)
	}
	return M
}

func (rs *rootSystem) apply(M mat, v vec) vec {
	out := make(vec, rs.n)
	for j := 0; j < rs.n; j++ {
		if v[j] != 0 {
			for i := 0; i < rs.n; i++ {
				out[i] += M[j][i] * v[j]
			}
		}
	}
	return out
}

func (rs *rootSystem) mult(A, B mat) mat {
	out := make(mat, rs.n)
	for j := 0; j < rs.n; j++ {
		out[j] = rs.apply(A, B[j])
	}
	return out
}

func reverseInts(w []int) []int {
	out := make([]int, len(w))
	for i, x := range w {
		out[len(w)-1-i] = x
	}
	return out
}

func (rs *rootSystem) inverse(M mat) mat {
	word := rs.reducedWord(M)
	return rs.wordMatrix(reverseInts(word))
}

// rightDescents returns the sorted right descent set.
func (rs *rootSystem) rightDescents(M mat) []int {
	var out []int
	for s := 0; s < rs.n; s++ {
		if isNegative(M[s]) {
			out = append(out, s)
		}
	}
	return out
}

func (rs *rootSystem) leftDescents(M mat) []int {
	return rs.rightDescents(rs.inverse(M))
}

func (rs *rootSystem) length(M mat) int {
	c := 0
	for _, r := range rs.positiveRoots() {
		if isNegative(rs.apply(M, r)) {
			c++
		}
	}
	return c
}

func (rs *rootSystem) inversions(M mat) bits {
	var b bits
	for i, r := range rs.positiveRoots() {
		if isNegative(rs.apply(M, r)) {
			b.set(i)
		}
	}
	return b
}

func (rs *rootSystem) rightTerminal(M mat) bool {
	for _, s := range rs.rightDescents(M) {
		for _, t := range rs.adj[s] {
			sum := make(vec, rs.n)
			for i := range sum {
				sum[i] = M[s][i] + M[t][i]
			}
			if !isPositive(sum) {
				return false
			}
		}
	}
	return true
}

func (rs *rootSystem) terminal(M mat) bool {
	return rs.rightTerminal(M) && rs.rightTerminal(rs.inverse(M))
}

// reducedWord strips right descents (smallest first) and returns the word of M.
func (rs *rootSystem) reducedWord(M mat) []int {
	var word []int
	for {
		R := rs.rightDescents(M)
		if len(R) == 0 {
			break
		}
		s := R[0]
		word = append(word, s)
		M = rs.rightMult(M, s)
	}
	return reverseInts(word)
}

func (rs *rootSystem) reflection(beta vec) mat {
	cols := make(mat, rs.n)
	for j := 0; j < rs.n; j++ {
		e := make(vec, rs.n)
		e[j] = 1
		c := rs.pairing(e, beta)
		col := make(vec, rs.n)
		for i := range col {
			col[i] = e[i] - c*beta[i]
		}
		cols[j] = col
	}
	return cols
}

func (rs *rootSystem) longestElement(J []int) mat {
	M := rs.identity
	for {
		asc := -1
		for _, s := range J {
			if !isNegative(M[s]) {
				asc = s
				break
			}
		}
		if asc < 0 {
			return M
		}
		M = rs.rightMult(M, asc)
	}
}

// isPrefix: U <= W in right weak order iff N(U^-1) is contained in N(W^-1).
func (rs *rootSystem) isPrefix(U, W mat) bool {
	return rs.inversions(rs.inverse(U)).subsetOf(rs.inversions(rs.inverse(W)))
}

func (rs *rootSystem) isSuffix(V, W mat) bool {
	return rs.inversions(V).subsetOf(rs.inversions(W))
}

// suffixList is the set of all V with W = U V length-additively, in the
// breadth-first insertion order of the original (left descents ascending).
type suffixList struct {
	order []mat
	u     map[string]mat
}

func (rs *rootSystem) suffixes(W mat) *suffixList {
	out := &suffixList{u: map[string]mat{keyMat(W): rs.identity}}
	out.order = []mat{W}
	for head := 0; head < len(out.order); head++ {
		V := out.order[head]
		for _, s := range rs.leftDescents(V) {
			Y := rs.leftMult(s, V)
			k := keyMat(Y)
			if _, ok := out.u[k]; !ok {
				out.u[k] = rs.mult(out.u[keyMat(V)], rs.wordMatrix([]int{s}))
				out.order = append(out.order, Y)
			}
		}
	}
	return out
}

func en(n int) *rootSystem {
	var edges [][2]int
	for i := 0; i < n-2; i++ {
		edges = append(edges, [2]int{i, i + 1})
	}
	edges = append(edges, [2]int{2, n - 1})
	return newRootSystem(n, edges)
}

// rankOf is the exact rank over the rationals.
func rankOf(rows [][]int) int {
	if len(rows) == 0 {
		return 0
	}
	A := make([][]*big.Rat, len(rows))
	for i, r := range rows {
		A[i] = make([]*big.Rat, len(r))
		for j, x := range r {
			A[i][j] = new(big.Rat).SetInt64(int64(x))
		}
	}
	rk := 0
	for c := 0; c < len(A[0]); c++ {
		p := -1
		for r := rk; r < len(A); r++ {
			if A[r][c].Sign() != 0 {
				p = r
				break
			}
		}
		if p < 0 {
			continue
		}
		A[rk], A[p] = A[p], A[rk]
		piv := new(big.Rat).Set(A[rk][c])
		for j := range A[rk] {
			A[rk][j] = new(big.Rat).Quo(A[rk][j], piv)
		}
		for r := range A {
			if r != rk && A[r][c].Sign() != 0 {
				f := new(big.Rat).Set(A[r][c])
				for j := range A[r] {
					A[r][j] = new(big.Rat).Sub(A[r][j], new(big.Rat).Mul(f, A[rk][j]))
				}
			}
		}
		rk++
	}
	return rk
}
