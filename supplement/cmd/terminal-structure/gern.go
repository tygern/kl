package main

import (
	"fmt"
	"sort"
	"strconv"
	"strings"
)

// ---------- Gern's D_m ----------

// gernInterval is Gern Definition 2.3.1.
func gernInterval(i, j int) []int {
	if 0 <= j && j < i && i >= 2 {
		return reverseInts(gernInterval(j, i))
	}
	if i == 1 && j >= 3 {
		out := []int{1}
		for x := 3; x <= j; x++ {
			out = append(out, x)
		}
		return out
	}
	if i == 0 && j >= 2 {
		var out []int
		for x := 1; x <= j; x++ {
			out = append(out, x)
		}
		return out
	}
	if 2 <= i && i <= j {
		var out []int
		for x := i; x <= j; x++ {
			out = append(out, x)
		}
		return out
	}
	panic(fmt.Sprintf("gern_interval: invalid arguments (%d, %d)", i, j))
}

// gernWnWord is the bracket form of w_n, Gern Lemma 2.3.4 (n even).
func gernWnWord(n int) []int {
	if n%2 != 0 {
		panic("gern_wn_word: n must be even")
	}
	k := n/2 - 2
	var word []int
	for i := 2; i <= n; i += 2 {
		word = append(word, gernInterval(i, 0)...)
	}
	for i := 0; i <= k; i++ {
		word = append(word, gernInterval(n-k+i, n-2*k+2*i)...)
	}
	return word
}

func gernToEn(label, n int) int {
	if label == 1 {
		return 1
	}
	if label == 2 {
		return n - 1
	}
	return label - 1
}

func absInt(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func sgn(x int) int {
	if x > 0 {
		return 1
	}
	return -1
}

// signedPerm is the one-line signed permutation of a Gern word (left action,
// composed left to right).
func signedPerm(word []int, m int) []int {
	gen := func(label int) []int {
		p := make([]int, m)
		for i := range p {
			p[i] = i + 1
		}
		if label == 1 {
			p[0], p[1] = -2, -1
		} else {
			p[label-2], p[label-1] = label, label-1
		}
		return p
	}
	w := make([]int, m)
	for i := range w {
		w[i] = i + 1
	}
	for _, s := range word {
		g := gen(s)
		nw := make([]int, m)
		for i, a := range g {
			nw[i] = sgn(a) * w[absInt(a)-1]
		}
		w = nw
	}
	return w
}

// gernCor2219 is the signed permutation of w_n inside D_m, Gern Corollary 2.2.19.
func gernCor2219(n, m int) []int {
	w := make([]int, 0, m)
	for i := 1; i <= m; i++ {
		switch {
		case i > n:
			w = append(w, i)
		case i == 1:
			if (n/2)%2 == 0 {
				w = append(w, 1)
			} else {
				w = append(w, -1)
			}
		case i%2 == 1:
			w = append(w, i)
		default:
			w = append(w, -(n + 2 - i))
		}
	}
	return w
}

func keyInts(w []int) string {
	parts := make([]string, len(w))
	for i, x := range w {
		parts[i] = strconv.Itoa(x)
	}
	return strings.Join(parts, ",")
}

// dTerm / dRoot: a root of D_m is c1 e_i1 + c2 e_i2 (coordinates 1-based).
type dRoot struct{ i1, c1, i2, c2 int }

func dImage(w []int, r dRoot) dRoot {
	j1, j2 := w[r.i1-1], w[r.i2-1]
	return dRoot{absInt(j1), sgn(j1) * r.c1, absInt(j2), sgn(j2) * r.c2}
}

func dPositive(r dRoot) bool {
	if r.c1 > 0 && r.c2 > 0 {
		return true
	}
	if r.c1 < 0 && r.c2 < 0 {
		return false
	}
	hi, lo := r.i2, r.i1
	if r.c1 > 0 {
		hi, lo = r.i1, r.i2
	}
	return hi > lo
}

// enumerateBadDm lists all noncommuting terminal elements of D_m as signed
// permutations (Gern labels), sorted by (length, one-line word).
type badElt struct {
	length int
	w      []int
}

func enumerateBadDm(m int) []badElt {
	simple := make([]map[int]int, m+1)
	simple[1] = map[int]int{1: 1, 2: 1}
	for i := 2; i <= m; i++ {
		simple[i] = map[int]int{i: 1, i - 1: -1}
	}
	adj := make([][]int, m+1)
	addEdge := func(a, b int) {
		adj[a] = append(adj[a], b)
		adj[b] = append(adj[b], a)
	}
	addEdge(1, 3)
	addEdge(2, 3)
	for i := 3; i < m; i++ {
		addEdge(i, i+1)
	}
	toRoot := func(d map[int]int) dRoot {
		var idx []int
		for k, v := range d {
			if v != 0 {
				idx = append(idx, k)
			}
		}
		if len(idx) != 2 {
			panic("root of D_m must have exactly two nonzero coordinates")
		}
		sort.Ints(idx)
		return dRoot{idx[0], d[idx[0]], idx[1], d[idx[1]]}
	}
	simpleRoot := make([]dRoot, m+1)
	for s := 1; s <= m; s++ {
		simpleRoot[s] = toRoot(simple[s])
	}
	sumRoot := make([]map[int]dRoot, m+1)
	for s := 1; s <= m; s++ {
		sumRoot[s] = map[int]dRoot{}
		for _, t := range adj[s] {
			r := map[int]int{}
			for k, v := range simple[s] {
				r[k] += v
			}
			for k, v := range simple[t] {
				r[k] += v
			}
			sumRoot[s][t] = toRoot(r)
		}
	}
	var posRoots []dRoot
	for i := 1; i <= m; i++ {
		for j := i + 1; j <= m; j++ {
			posRoots = append(posRoots, dRoot{i, -1, j, 1}, dRoot{i, 1, j, 1})
		}
	}
	descents := func(w []int) []int {
		var out []int
		for s := 1; s <= m; s++ {
			if !dPositive(dImage(w, simpleRoot[s])) {
				out = append(out, s)
			}
		}
		return out
	}
	inv := func(w []int) []int {
		out := make([]int, m)
		for i, j := range w {
			if j > 0 {
				out[j-1] = i + 1
			} else {
				out[-j-1] = -(i + 1)
			}
		}
		return out
	}
	rightTerminal := func(w []int) bool {
		for _, s := range descents(w) {
			for _, t := range adj[s] {
				if !dPositive(dImage(w, sumRoot[s][t])) {
					return false
				}
			}
		}
		return true
	}
	length := func(w []int) int {
		c := 0
		for _, r := range posRoots {
			if !dPositive(dImage(w, r)) {
				c++
			}
		}
		return c
	}
	commutingProduct := func(w []int) bool {
		R := descents(w)
		inR := map[int]bool{}
		for _, s := range R {
			inR[s] = true
		}
		for _, s := range R {
			for _, t := range adj[s] {
				if inR[t] {
					return false
				}
			}
		}
		return length(w) == len(R) && keyInts(signedPerm(R, m)) == keyInts(w)
	}

	var bad []badElt
	count := 0
	perm := make([]int, m)
	for i := range perm {
		perm[i] = i + 1
	}
	w := make([]int, m)
	for {
		for mask := 0; mask < 1<<uint(m); mask++ {
			neg := 0
			for b := 0; b < m; b++ {
				if mask>>uint(b)&1 == 1 {
					neg++
				}
			}
			if neg%2 == 1 {
				continue
			}
			for b := 0; b < m; b++ {
				if mask>>uint(b)&1 == 1 {
					w[b] = -perm[b]
				} else {
					w[b] = perm[b]
				}
			}
			count++
			if rightTerminal(w) && rightTerminal(inv(w)) && !commutingProduct(w) {
				bad = append(bad, badElt{length(w), append([]int(nil), w...)})
			}
		}
		if !nextPermutation(perm) {
			break
		}
	}
	fact := 1
	for i := 2; i <= m; i++ {
		fact *= i
	}
	if count != (1<<uint(m-1))*fact {
		panic("enumerate_bad_dm: wrong number of elements of D_m")
	}
	sort.Slice(bad, func(a, b int) bool {
		if bad[a].length != bad[b].length {
			return bad[a].length < bad[b].length
		}
		return cmpVec(bad[a].w, bad[b].w) < 0
	})
	return bad
}

// nextPermutation advances p to the next lexicographic permutation.
func nextPermutation(p []int) bool {
	i := len(p) - 2
	for i >= 0 && p[i] >= p[i+1] {
		i--
	}
	if i < 0 {
		return false
	}
	j := len(p) - 1
	for p[j] <= p[i] {
		j--
	}
	p[i], p[j] = p[j], p[i]
	for a, b := i+1, len(p)-1; a < b; a, b = a+1, b-1 {
		p[a], p[b] = p[b], p[a]
	}
	return true
}

// gernPred is one element predicted by Gern Theorem 2.3.6.
type gernPred struct {
	name string
	k    int
	U    []int
	perm []int
}

// gernPredicted: bad elements of D_m are w_k u, k even, 4 <= k <= m, u a
// commuting product of generators neither in nor adjacent to supp w_k =
// {1,...,k}.  Returned sorted by name.
func gernPredicted(m int) []gernPred {
	var out []gernPred
	for k := 4; k <= m; k += 2 {
		var rest []int
		for x := k + 2; x <= m; x++ {
			rest = append(rest, x)
		}
		var rec func(start int, cur []int)
		rec = func(start int, cur []int) {
			ok := true
			for _, a := range cur {
				for _, b := range cur {
					if absInt(a-b) == 1 {
						ok = false
					}
				}
			}
			if ok {
				name := fmt.Sprintf("w%d", k)
				for _, u := range cur {
					name += fmt.Sprintf("*s%d", u)
				}
				U := append([]int(nil), cur...)
				word := append(gernWnWord(k), U...)
				out = append(out, gernPred{name, k, U, signedPerm(word, m)})
			}
			for i := start; i < len(rest); i++ {
				rec(i+1, append(cur, rest[i]))
			}
		}
		rec(0, nil)
	}
	sort.Slice(out, func(a, b int) bool { return out[a].name < out[b].name })
	return out
}
