// Package indefinite is the exact integer arithmetic for generalized
// E-type (T-shaped) Coxeter diagrams of arbitrary rank, used for the all-k
// audits of full-support reflection families.
//
// It ports the library part of research/en_affine_referee/verify_indefinite.py
// (the E13 main of that file is deliberately not ported). It imports only the
// standard library and shares no code with any family engine.
//
// All arithmetic is int64 with explicit overflow guards (panics). Failed
// assertions panic with an *AssertionError; commands turn those into a clear
// message and a non-zero exit status.
package indefinite

import (
	"fmt"
	"sort"
)

// Vec is an integer vector; Mat is a row-major integer matrix.
type Vec = []int64
type Mat = [][]int64

// AssertionError is the panic value used for every failed assertion.
type AssertionError struct{ Msg string }

func (e *AssertionError) Error() string { return "assertion failed: " + e.Msg }

// Assert panics with an AssertionError when cond is false.
func Assert(cond bool, format string, args ...any) {
	if !cond {
		panic(&AssertionError{Msg: fmt.Sprintf(format, args...)})
	}
}

func add64(a, b int64) int64 {
	c := a + b
	if (c > a) != (b > 0) {
		panic("int64 overflow in addition")
	}
	return c
}

func mul64(a, b int64) int64 {
	if a == 0 || b == 0 {
		return 0
	}
	c := a * b
	if c/b != a || (a == -1 && b == -1<<63) || (b == -1 && a == -1<<63) {
		panic("int64 overflow in multiplication")
	}
	return c
}

// Ctx holds the diagram for one rank: path 0-1-...-(N-2) with an extra node
// N-1 attached to node 2.
type Ctx struct {
	N     int
	Edges [][2]int
	Adj   [][]int // sorted adjacency lists
	E     Mat     // identity
	Zero  Mat
}

// Configure selects the generalized E diagram of the given rank
// (python: configure(rank)).
func Configure(rank int) *Ctx {
	c := &Ctx{N: rank}
	for i := 0; i < rank-2; i++ {
		c.Edges = append(c.Edges, [2]int{i, i + 1})
	}
	c.Edges = append(c.Edges, [2]int{2, rank - 1})
	c.Adj = make([][]int, rank)
	for _, e := range c.Edges {
		c.Adj[e[0]] = append(c.Adj[e[0]], e[1])
		c.Adj[e[1]] = append(c.Adj[e[1]], e[0])
	}
	for i := range c.Adj {
		sort.Ints(c.Adj[i])
	}
	c.E = make(Mat, rank)
	c.Zero = make(Mat, rank)
	for i := 0; i < rank; i++ {
		c.E[i] = make([]int64, rank)
		c.Zero[i] = make([]int64, rank)
		c.E[i][i] = 1
	}
	return c
}

// VecEq reports whether two vectors are equal.
func VecEq(a, b Vec) bool {
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

// MatEq reports whether two matrices are equal.
func MatEq(a, b Mat) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !VecEq(a[i], b[i]) {
			return false
		}
	}
	return true
}

// Pair is the symmetric bilinear form sum_i a_i (2 b_i - sum_{j~i} b_j).
func (c *Ctx) Pair(a, b Vec) int64 {
	var total int64
	for i := 0; i < c.N; i++ {
		s := mul64(2, b[i])
		for _, j := range c.Adj[i] {
			s = add64(s, -b[j])
		}
		total = add64(total, mul64(a[i], s))
	}
	return total
}

// PairingVector returns (<alpha_s, a>)_s.
func (c *Ctx) PairingVector(a Vec) Vec {
	out := make(Vec, c.N)
	for s := 0; s < c.N; s++ {
		out[s] = c.Pair(c.E[s], a)
	}
	return out
}

// MV is the matrix-vector product.
func (c *Ctx) MV(a Mat, v Vec) Vec {
	out := make(Vec, len(a))
	for i, row := range a {
		var s int64
		for j := range row {
			s = add64(s, mul64(row[j], v[j]))
		}
		out[i] = s
	}
	return out
}

// MM is the N x N matrix product.
func (c *Ctx) MM(a, b Mat) Mat {
	out := make(Mat, c.N)
	for i := 0; i < c.N; i++ {
		out[i] = make([]int64, c.N)
		for j := 0; j < c.N; j++ {
			var s int64
			for k := 0; k < c.N; k++ {
				s = add64(s, mul64(a[i][k], b[k][j]))
			}
			out[i][j] = s
		}
	}
	return out
}

// Add returns a + k*b.
func Add(a, b Vec, k int64) Vec {
	out := make(Vec, len(a))
	for i := range a {
		out[i] = add64(a[i], mul64(k, b[i]))
	}
	return out
}

// Simple applies the simple reflection s to the root vector a.
func (c *Ctx) Simple(a Vec, s int) Vec {
	b := append(Vec(nil), a...)
	v := -a[s]
	for _, t := range c.Adj[s] {
		v = add64(v, a[t])
	}
	b[s] = v
	return b
}

// Witness is the lowering sequence certifying that a root is real.
type Witness struct {
	Height             int64 `json:"height"`
	LoweringGenerators []int `json:"lowering_generators"`
	SimpleRoot         int   `json:"simple_root"`
}

func minOf(a Vec) int64 {
	m := a[0]
	for _, x := range a {
		if x < m {
			m = x
		}
	}
	return m
}

// Min returns the minimum entry of a.
func Min(a Vec) int64 { return minOf(a) }

func sum(a Vec) int64 {
	var s int64
	for _, x := range a {
		s = add64(s, x)
	}
	return s
}

func (c *Ctx) simpleIndex(a Vec) int {
	for i := 0; i < c.N; i++ {
		if VecEq(a, c.E[i]) {
			return i
		}
	}
	return -1
}

// RealRootWitness returns an actual sequence of simple reflections reducing a
// to a simple root.
func (c *Ctx) RealRootWitness(a Vec) Witness {
	Assert(minOf(a) >= 0, "real_root_witness: negative entry")
	original := a
	steps := []int{}
	for c.simpleIndex(a) < 0 {
		choice := -1
		for s := 0; s < c.N; s++ {
			if c.Pair(c.E[s], a) > 0 {
				choice = s
				break
			}
		}
		Assert(choice >= 0, "not certified real %v", a)
		b := c.Simple(a, choice)
		Assert(minOf(b) >= 0 && sum(b) < sum(a), "lowering step not decreasing")
		steps = append(steps, choice)
		a = b
	}
	Assert(c.Pair(original, original) == 2, "root norm is not 2")
	return Witness{Height: sum(original), LoweringGenerators: steps, SimpleRoot: c.simpleIndex(a)}
}

// Reflection returns the matrix of the reflection in the root a.
func (c *Ctx) Reflection(a Vec) Mat {
	p := c.PairingVector(a)
	out := make(Mat, c.N)
	for i := 0; i < c.N; i++ {
		out[i] = make([]int64, c.N)
		for j := 0; j < c.N; j++ {
			out[i][j] = add64(c.E[i][j], -mul64(a[i], p[j]))
		}
	}
	return out
}

// Right multiplies the matrix a on the right by the simple reflection s.
func (c *Ctx) Right(a Mat, s int) Mat {
	b := make(Mat, c.N)
	for i := 0; i < c.N; i++ {
		b[i] = append([]int64(nil), a[i]...)
		b[i][s] = -a[i][s]
		for _, t := range c.Adj[s] {
			b[i][t] = add64(a[i][t], a[i][s])
		}
	}
	return b
}

// ReducedMatrix multiplies out a word, asserting that it is reduced.
func (c *Ctx) ReducedMatrix(word []int) Mat {
	a := c.E
	for _, s := range word {
		for i := 0; i < c.N; i++ {
			Assert(a[i][s] >= 0, "nonreduced supplied word")
		}
		a = c.Right(a, s)
	}
	return a
}
