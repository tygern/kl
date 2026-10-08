// Package parabolic is the shared exact-arithmetic library of the two
// parabolic terminal enumerators for the finite Weyl groups E6, E7 and E8.
//
// It ports the library part of research/en_e8/parabolic_terminals.cpp (the
// part that research/en_e8/recursive_terminals.cpp includes as a header):
// root arithmetic, the packed one-line state of a Weyl group element (the
// images of the simple roots, one byte each), right and left multiplication
// by simple reflections, descent sets, the one-sided terminal enumeration of
// a standard parabolic subgroup, minimal coset representatives, and the
// fully commutative elements with their star-operation closure check.
//
// The package imports only the Go standard library. It deliberately knows
// nothing about the other engines of the supplement.
package parabolic

import (
	"fmt"
	"math/bits"
	"os"
	"sort"
)

// State packs the images of the simple roots a_0..a_{n-1} (root ids, one byte
// each) of a Weyl group element w: byte s is the id of w(a_s).
type State = uint64

// Element is a group element with its inverse and a reduced word.
type Element struct {
	W    State
	Inv  State
	Word []int
}

// Failure is the panic value used for every failed assertion.
type Failure struct{ Msg string }

func (f Failure) Error() string { return "assertion failed: " + f.Msg }

// Assert panics with a Failure when cond is false.
func Assert(cond bool, format string, args ...any) {
	if !cond {
		panic(Failure{fmt.Sprintf(format, args...)})
	}
}

// ExitOnFailure must be deferred directly from main; it turns a Failure
// panic into a clear message on stderr and exit status 1.
func ExitOnFailure() {
	if r := recover(); r != nil {
		if f, ok := r.(Failure); ok {
			fmt.Fprintln(os.Stderr, f.Error())
		} else {
			fmt.Fprintln(os.Stderr, "fatal error:", r)
		}
		os.Exit(1)
	}
}

// MulChecked multiplies non-negative int64 values, panicking on overflow.
func MulChecked(a, b int64) int64 {
	if a < 0 || b < 0 {
		panic(Failure{"negative operand in MulChecked"})
	}
	if a != 0 && b > (1<<62)/a {
		panic(Failure{"int64 overflow in multiplication"})
	}
	return a * b
}

// Group holds the root system data of E_n (n = 6, 7, 8) in the Bourbaki-like
// labelling used by the original: the chain 0-1-...-(n-2) with the extra
// node n-1 attached to node 2.
type Group struct {
	N        int
	Omitted  int
	J        uint32 // parabolic generator mask: all generators except Omitted
	Adj      [8][]int
	Edges    [][2]int
	Roots    [][8]int
	rid      map[[8]int]int
	adds     []int16 // adds[i*R+j] = id of Roots[i]+Roots[j], or -1
	neg      []uint8
	reflect  [8][]uint8 // reflect[s][r] = id of s_s(Roots[r])
	positive []bool
	Identity State

	// FCStarChecks is the number of star-operation checks performed by the
	// last call of FullyCommutative.
	FCStarChecks uint64
}

// Image returns the root id of w(a_s).
func Image(w State, s int) int { return int(w>>(8*uint(s))) & 255 }

// Replace returns w with the image of a_s replaced by value.
func Replace(w State, s int, value int) State {
	Assert(value >= 0 && value < 256, "root id %d out of range", value)
	return (w &^ (State(255) << (8 * uint(s)))) | (State(value) << (8 * uint(s)))
}

// New builds the root system data for E_n, asserting the root count.
func New(n int) *Group {
	Assert(n >= 6 && n <= 8, "rank must be 6, 7 or 8")
	g := &Group{N: n}
	for s := 0; s < n-2; s++ {
		g.Edges = append(g.Edges, [2]int{s, s + 1})
	}
	g.Edges = append(g.Edges, [2]int{2, n - 1})
	for _, e := range g.Edges {
		g.Adj[e[0]] = append(g.Adj[e[0]], e[1])
		g.Adj[e[1]] = append(g.Adj[e[1]], e[0])
	}
	g.Omitted = n - 2
	g.J = (uint32(1)<<uint(n) - 1) ^ (uint32(1) << uint(g.Omitted))
	g.rid = map[[8]int]int{}
	for s := 0; s < n; s++ {
		var a [8]int
		a[s] = 1
		g.rid[a] = len(g.Roots)
		g.Roots = append(g.Roots, a)
	}
	for i := 0; i < len(g.Roots); i++ {
		for s := 0; s < n; s++ {
			b := g.Roots[i]
			b[s] = -b[s]
			for _, t := range g.Adj[s] {
				b[s] += g.Roots[i][t]
			}
			if _, ok := g.rid[b]; !ok {
				g.rid[b] = len(g.Roots)
				g.Roots = append(g.Roots, b)
			}
		}
	}
	want := 240
	if n == 6 {
		want = 72
	} else if n == 7 {
		want = 126
	}
	Assert(len(g.Roots) == want, "root count %d, expected %d", len(g.Roots), want)
	R := len(g.Roots)
	g.positive = make([]bool, R)
	g.neg = make([]uint8, R)
	g.adds = make([]int16, R*R)
	for s := 0; s < n; s++ {
		g.reflect[s] = make([]uint8, R)
	}
	for i := 0; i < R; i++ {
		b := g.Roots[i]
		g.positive[i] = true
		for s := 0; s < n; s++ {
			g.positive[i] = g.positive[i] && b[s] >= 0
			b[s] = -b[s]
		}
		id, ok := g.rid[b]
		Assert(ok, "negative root missing")
		g.neg[i] = uint8(id)
		for s := 0; s < n; s++ {
			b = g.Roots[i]
			b[s] = -b[s]
			for _, t := range g.Adj[s] {
				b[s] += g.Roots[i][t]
			}
			id, ok := g.rid[b]
			Assert(ok, "reflected root missing")
			g.reflect[s][i] = uint8(id)
		}
		for j := 0; j < R; j++ {
			var c [8]int
			for s := 0; s < n; s++ {
				c[s] = g.Roots[i][s] + g.Roots[j][s]
			}
			if id, ok := g.rid[c]; ok {
				g.adds[i*R+j] = int16(id)
			} else {
				g.adds[i*R+j] = -1
			}
		}
	}
	g.Identity = 0
	for s := 0; s < n; s++ {
		g.Identity = Replace(g.Identity, s, s)
	}
	return g
}

// NumRoots returns the number of roots.
func (g *Group) NumRoots() int { return len(g.Roots) }

// Positive reports whether root id r is a positive root.
func (g *Group) Positive(r int) bool { return g.positive[r] }

// ReflectRoot returns the id of s_s(Roots[r]).
func (g *Group) ReflectRoot(s, r int) int { return int(g.reflect[s][r]) }

// Right returns w*s (right multiplication by the simple reflection s).
func (g *Group) Right(w State, s int) State {
	a := Image(w, s)
	v := Replace(w, s, int(g.neg[a]))
	R := len(g.Roots)
	for _, t := range g.Adj[s] {
		v = Replace(v, t, int(g.adds[a*R+Image(w, t)]))
	}
	return v
}

// Left returns s*w.
func (g *Group) Left(w State, s int) State {
	var v State
	for t := 0; t < g.N; t++ {
		v = Replace(v, t, int(g.reflect[s][Image(w, t)]))
	}
	return v
}

// Desc returns the right descent mask of w.
func (g *Group) Desc(w State) uint32 {
	var result uint32
	for s := 0; s < g.N; s++ {
		if !g.positive[Image(w, s)] {
			result |= 1 << uint(s)
		}
	}
	return result
}

// Commuting reports whether the generators in mask pairwise commute.
func (g *Group) Commuting(mask uint32) bool {
	for _, e := range g.Edges {
		if mask&(1<<uint(e[0])) != 0 && mask&(1<<uint(e[1])) != 0 {
			return false
		}
	}
	return true
}

// WeakRight reports whether w is right terminal with respect to the
// generators in mask.
func (g *Group) WeakRight(w State, mask uint32) bool {
	d := g.Desc(w) & mask
	R := len(g.Roots)
	for s := 0; s < g.N; s++ {
		if d&(1<<uint(s)) == 0 {
			continue
		}
		for _, t := range g.Adj[s] {
			if mask&(1<<uint(t)) != 0 {
				z := int(g.adds[Image(w, s)*R+Image(w, t)])
				Assert(z >= 0, "weakright: sum of roots missing")
				if !g.positive[z] {
					return false
				}
			}
		}
	}
	return true
}

// ReducedWord returns the reduced word of w obtained by repeatedly removing
// the lowest right descent.
func (g *Group) ReducedWord(w State) []int {
	var reverse []int
	for w != g.Identity {
		d := g.Desc(w)
		Assert(d != 0, "reduced_word: no descent")
		s := bits.TrailingZeros32(d)
		reverse = append(reverse, s)
		w = g.Right(w, s)
		Assert(len(reverse) <= 120, "reduced_word: word too long")
	}
	out := make([]int, len(reverse))
	for i, s := range reverse {
		out[len(reverse)-1-i] = s
	}
	return out
}

// Inverse returns the element of the inverse of the product of word.
func (g *Group) Inverse(word []int) State {
	v := g.Identity
	for i := len(word) - 1; i >= 0; i-- {
		v = g.Right(v, word[i])
	}
	return v
}

// Leq is the Bruhat order test x <= w (lx, lw are the lengths).
func (g *Group) Leq(x State, lx int, w State, lw int) bool {
	for x != w {
		if lx >= lw {
			return false
		}
		d := g.Desc(w)
		Assert(d != 0, "leq: no descent")
		s := bits.TrailingZeros32(d)
		if !g.positive[Image(x, s)] {
			x = g.Right(x, s)
			lx--
		}
		w = g.Right(w, s)
		lw--
	}
	return true
}

// ParabolicOneSided enumerates the standard parabolic subgroup W_J by
// breadth-first search and returns its right-terminal elements (those w with
// WeakRight(w, J)) together with the order of the subgroup.
func (g *Group) ParabolicOneSided() ([]Element, uint64) {
	n := g.N
	expected := 2903040
	if n == 6 {
		expected = 1920
	} else if n == 7 {
		expected = 51840
	}
	group := make([]State, 1, expected)
	inverses := make([]State, 1, expected)
	group[0], inverses[0] = g.Identity, g.Identity
	parents := make([]uint32, 1, expected)
	last := make([]uint8, 1, expected)
	lengths := make([]uint8, 1, expected)
	ids := make(map[State]uint32, expected)
	ids[g.Identity] = 0
	var accepted []Element
	for w := 0; w < len(group); w++ {
		for s := 0; s < n; s++ {
			if g.J&(1<<uint(s)) == 0 {
				continue
			}
			v := g.Right(group[w], s)
			if _, ok := ids[v]; !ok {
				id := uint32(len(group))
				ids[v] = id
				group = append(group, v)
				inverses = append(inverses, g.Left(inverses[w], s))
				parents = append(parents, uint32(w))
				last = append(last, uint8(s))
				lengths = append(lengths, lengths[w]+1)
			}
		}
		if g.WeakRight(group[w], g.J) {
			var word []int
			for p := uint32(w); p != 0; p = parents[p] {
				word = append(word, int(last[p]))
			}
			for i, j := 0, len(word)-1; i < j; i, j = i+1, j-1 {
				word[i], word[j] = word[j], word[i]
			}
			Assert(len(word) == int(lengths[w]), "parabolic word length mismatch")
			accepted = append(accepted, Element{group[w], inverses[w], word})
		}
	}
	Assert(len(group) == expected, "parabolic order %d, expected %d", len(group), expected)
	fmt.Fprintf(os.Stderr, "Parabolic order %d, one-sided terminals %d\n", len(group), len(accepted))
	return accepted, uint64(len(group))
}

// MinimalCosets returns the minimal representatives of the cosets x*W_H for
// x in W_K (H a subset of K), generated by left multiplication by generators
// of K and reduction by right descents in H. expected is the asserted count.
func (g *Group) MinimalCosets(K, H uint32, expected int) []Element {
	reps := []State{g.Identity}
	ids := map[State]int{g.Identity: 0}
	for i := 0; i < len(reps); i++ {
		for s := 0; s < g.N; s++ {
			if K&(1<<uint(s)) == 0 {
				continue
			}
			a := g.Left(reps[i], s)
			for g.Desc(a)&H != 0 {
				a = g.Right(a, bits.TrailingZeros32(g.Desc(a)&H))
			}
			if _, ok := ids[a]; !ok {
				ids[a] = len(reps)
				reps = append(reps, a)
			}
		}
	}
	Assert(len(reps) == expected, "coset count %d, expected %d", len(reps), expected)
	var result []Element
	for _, a := range reps {
		Assert(g.Desc(a)&H == 0, "representative is not minimal")
		word := g.ReducedWord(a)
		for _, s := range word {
			Assert(K&(1<<uint(s)) != 0, "representative word leaves W_K")
		}
		result = append(result, Element{a, g.Inverse(word), word})
	}
	return result
}

// Cosets returns the minimal right coset representatives of W_J in W.
func (g *Group) Cosets() []Element {
	expected := 240
	if g.N == 6 {
		expected = 27
	} else if g.N == 7 {
		expected = 56
	}
	all := uint32(1)<<uint(g.N) - 1
	result := g.MinimalCosets(all, g.J, expected)
	fmt.Fprintf(os.Stderr, "Minimal right cosets %d\n", len(result))
	return result
}

// RootPermutation returns the permutation of root ids induced by the element
// whose reduced word is word (letters applied from the right end), asserting
// that it agrees with the images of the simple roots stored in a.
func (g *Group) RootPermutation(a Element) []uint8 {
	R := len(g.Roots)
	permutation := make([]uint8, R)
	for r := 0; r < R; r++ {
		permutation[r] = uint8(r)
	}
	for i := len(a.Word) - 1; i >= 0; i-- {
		for r := 0; r < R; r++ {
			permutation[r] = g.reflect[a.Word[i]][permutation[r]]
		}
	}
	for s := 0; s < g.N; s++ {
		Assert(int(permutation[s]) == Image(a.W, s), "coset permutation disagrees with representative")
	}
	return permutation
}

// Compose applies a root permutation to every simple-root image of v.
func (g *Group) Compose(permutation []uint8, v State) State {
	var w State
	for s := 0; s < g.N; s++ {
		w = Replace(w, s, int(permutation[Image(v, s)]))
	}
	return w
}

// FullyCommutative enumerates all fully commutative elements and checks
// closure under the star operations. It records the number of star checks in
// g.FCStarChecks.
func (g *Group) FullyCommutative() []Element {
	n := g.N
	fc := []Element{{g.Identity, g.Identity, nil}}
	ids := map[State]int{g.Identity: 0}
	for i := 0; i < len(fc); i++ {
		for s := 0; s < n; s++ {
			// Appending a descent cannot discover a longer FC element.
			if !g.positive[Image(fc[i].W, s)] {
				continue
			}
			w := g.Right(fc[i].W, s)
			if _, ok := ids[w]; ok {
				continue
			}
			d := g.Desc(w)
			if !g.Commuting(d) {
				continue
			}
			ok := true
			for t := 0; t < n; t++ {
				if d&(1<<uint(t)) == 0 {
					continue
				}
				found, present := ids[g.Right(w, t)]
				if !present {
					ok = false
				} else {
					Assert(len(fc[found].Word) == len(fc[i].Word), "fc length mismatch")
				}
			}
			if !ok {
				continue
			}
			word := append(append([]int(nil), fc[i].Word...), s)
			inv := g.Left(fc[i].Inv, s)
			ids[w] = len(fc)
			fc = append(fc, Element{w, inv, word})
		}
	}
	want := 10846
	if n == 6 {
		want = 662
	} else if n == 7 {
		want = 2670
	}
	Assert(len(fc) == want, "fc count %d, expected %d", len(fc), want)
	var starChecks uint64
	for _, x := range fc {
		for side := 0; side < 2; side++ {
			for _, e := range g.Edges {
				s, t := e[0], e[1]
				w := x.W
				if side == 1 {
					w = x.Inv
				}
				mask := uint32(1)<<uint(s) | uint32(1)<<uint(t)
				if bits.OnesCount32(g.Desc(w)&mask) != 1 {
					continue
				}
				valid := 0
				for _, u := range [2]int{s, t} {
					v := g.Right(w, u)
					if bits.OnesCount32(g.Desc(v)&mask) == 1 {
						valid++
						_, present := ids[v]
						Assert(present, "star operation leaves the fully commutative set")
					}
				}
				Assert(valid == 1, "star operation not unique")
				starChecks++
			}
		}
	}
	fmt.Fprintf(os.Stderr, "Fully commutative elements %d, star checks %d\n", len(fc), starChecks)
	g.FCStarChecks = starChecks
	return fc
}

// FCIndex maps the State of every fully commutative element to its index.
func FCIndex(fc []Element) map[State]int {
	m := make(map[State]int, len(fc))
	for i := range fc {
		if _, ok := m[fc[i].W]; !ok {
			m[fc[i].W] = i
		}
	}
	return m
}

// SupportCommuting asserts the support condition for a commuting terminal
// word: its letters are distinct and pairwise commuting. Returns the support.
func (g *Group) SupportCommuting(word []int) uint32 {
	var support uint32
	for _, s := range word {
		support |= 1 << uint(s)
	}
	Assert(bits.OnesCount32(support) == len(word) && g.Commuting(support), "commuting terminal has wrong support")
	return support
}

// SortByLengthThenWord sorts elements by word length, then lexicographically.
func SortByLengthThenWord(es []Element) {
	sort.Slice(es, func(i, j int) bool {
		a, b := es[i].Word, es[j].Word
		if len(a) != len(b) {
			return len(a) < len(b)
		}
		for k := range a {
			if a[k] != b[k] {
				return a[k] < b[k]
			}
		}
		return false
	})
}

// Bottom is a fully commutative element below a non-commuting terminal.
type Bottom struct {
	Word   []int `json:"word"`
	Length int   `json:"length"`
	Rank   int   `json:"rank"`
}

// Bad is the certificate record of a non-commuting terminal element.
type Bad struct {
	Word    []int    `json:"word"`
	Length  int      `json:"length"`
	Rmask   uint32   `json:"Rmask"`
	Lmask   uint32   `json:"Lmask"`
	Bottoms []Bottom `json:"bottoms"`
}

// BadRecords builds the certificate records of the non-commuting terminals,
// listing for each the fully commutative elements below it with the same
// right and left descent containment (the "bottoms").
func (g *Group) BadRecords(bads, fc []Element) []Bad {
	out := make([]Bad, 0, len(bads))
	for _, b := range bads {
		rd, ld := g.Desc(b.W), g.Desc(b.Inv)
		rec := Bad{Word: append([]int{}, b.Word...), Length: len(b.Word), Rmask: rd, Lmask: ld, Bottoms: []Bottom{}}
		for _, x := range fc {
			if g.Desc(x.W)&rd == rd && g.Desc(x.Inv)&ld == ld && g.Leq(x.W, len(x.Word), b.W, len(b.Word)) {
				w := append([]int{}, x.Word...)
				rec.Bottoms = append(rec.Bottoms, Bottom{Word: w, Length: len(x.Word), Rank: len(b.Word) - len(x.Word)})
			}
		}
		out = append(out, rec)
	}
	return out
}
