// Command review-infinite preserves the independent checks made for the
// 10 October 2026 review. It verifies affine root data, cover braid witnesses,
// translation length, printed lowering words, and exact samples of the
// rank-uniform and E10 formulas. The samples do not prove parameter identities.
//
// This command uses only the standard library, no supplement code or data
// files. Small matrix calculations use bounded int64 arithmetic; large root
// samples and their norms use math/big, including all intermediate products.
//
// Run from the repository root:
//
//	go run -C tools ./cmd/review-infinite
package main

import (
	"flag"
	"fmt"
	"math/big"
	"os"
	"reflect"
	"sort"
	"strings"
)

type vector []int64
type matrix []vector // Columns in the simple-root basis.
type diagram [][]int
type bigVector []*big.Int

func check(ok bool, format string, args ...any) {
	if !ok {
		fmt.Fprintf(os.Stderr, "review-infinite: "+format+"\n", args...)
		os.Exit(1)
	}
}

func en(n int) diagram {
	a := make(diagram, n)
	add := func(s, t int) {
		a[s] = append(a[s], t)
		a[t] = append(a[t], s)
	}
	for s := 0; s < n-2; s++ {
		add(s, s+1)
	}
	add(2, n-1)
	return a
}

func adjacent(a diagram, s, t int) bool {
	for _, u := range a[s] {
		if u == t {
			return true
		}
	}
	return false
}

// Bound every small-arithmetic input/output. With rank <= 13, subsequent
// pairings and matrix products of bounded vectors are far inside int64.
func bounded(v vector) {
	check(len(v) <= 13, "small-arithmetic rank exceeds 13")
	for _, x := range v {
		check(x > -(1<<20) && x < 1<<20, "small-arithmetic coordinate exceeds bound: %v", v)
	}
}

func pairing(x, y vector, a diagram) int64 {
	bounded(x)
	bounded(y)
	var out int64
	for s := range x {
		c := 2 * y[s]
		for _, t := range a[s] {
			c -= y[t]
		}
		out += x[s] * c
	}
	return out
}

func reflectRoot(x vector, s int, a diagram) vector {
	bounded(x)
	y := append(vector(nil), x...)
	y[s] = -x[s]
	for _, t := range a[s] {
		y[s] += x[t]
	}
	bounded(y)
	return y
}

func positive(v vector) bool {
	bounded(v)
	pos, neg := false, false
	for _, x := range v {
		pos = pos || x > 0
		neg = neg || x < 0
	}
	check(pos != neg, "vector is not a signed root: %v", v)
	return pos
}

func identity(n int) matrix {
	m := make(matrix, n)
	for j := range m {
		m[j] = make(vector, n)
		m[j][j] = 1
	}
	return m
}

func apply(m matrix, x vector) vector {
	bounded(x)
	out := make(vector, len(x))
	for j, c := range x {
		bounded(m[j])
		for i, y := range m[j] {
			out[i] += c * y
		}
	}
	bounded(out)
	return out
}

func right(m matrix, s int, a diagram) matrix {
	out := make(matrix, len(m))
	bounded(m[s])
	for j, col := range m {
		bounded(col)
		out[j] = append(vector(nil), col...)
		for i, x := range col {
			if j == s {
				out[j][i] = -x
			} else if adjacent(a, s, j) {
				out[j][i] = x + m[s][i]
			}
		}
		bounded(out[j])
	}
	return out
}

func digits(s string) []int {
	out := make([]int, len(s))
	for j, c := range s {
		check(c >= '0' && c <= '9', "invalid digit %q", c)
		out[j] = int(c - '0')
	}
	return out
}

func wordMatrix(word []int, a diagram, requireReduced bool) matrix {
	m := identity(len(a))
	for j, s := range word {
		check(s >= 0 && s < len(a), "invalid generator %d", s)
		if requireReduced {
			check(positive(m[s]), "word is not reduced at position %d: %v", j, word)
		}
		m = right(m, s, a)
	}
	return m
}

func reduced(m matrix, a diagram) []int {
	var word []int
	unit := identity(len(a))
	for !reflect.DeepEqual(m, unit) {
		s := -1
		for j, col := range m {
			if !positive(col) {
				s = j
				break
			}
		}
		check(s >= 0, "nonidentity matrix has no descent")
		word = append(word, s)
		check(len(word) < 10000, "descent stripping exceeded length bound")
		m = right(m, s, a)
	}
	for i, j := 0, len(word)-1; i < j; i, j = i+1, j-1 {
		word[i], word[j] = word[j], word[i]
	}
	return word
}

// A convex sts chain in the dependency poset can be made consecutive by
// commutations. We verify the resulting order and its reduced word matrix.
func braidWord(word []int, a diagram) ([]int, int, bool) {
	n := len(word)
	pred := make([][]bool, n)
	for j, s := range word {
		pred[j] = make([]bool, n)
		for i := 0; i < j; i++ {
			if word[i] == s || adjacent(a, s, word[i]) {
				pred[j][i] = true
				for q := 0; q < i; q++ {
					pred[j][q] = pred[j][q] || pred[i][q]
				}
			}
		}
	}
	for j, s := range word {
		for i := 0; i < j; i++ {
			if !pred[j][i] || word[i] != s {
				continue
			}
			var between []int
			for k := 0; k < j; k++ {
				if pred[j][k] && pred[k][i] {
					between = append(between, k)
				}
			}
			if len(between) != 1 || !adjacent(a, s, word[between[0]]) {
				continue
			}
			k := between[0]
			var order []int
			for q := 0; q < j; q++ {
				if pred[j][q] && q != i && q != k {
					order = append(order, q)
				}
			}
			position := len(order)
			order = append(order, i, k, j)
			for q := 0; q < n; q++ {
				if !pred[j][q] && q != j {
					order = append(order, q)
				}
			}
			check(len(order) == n, "braid witness order has wrong length")
			seen := make([]bool, n)
			out := make([]int, n)
			for p, q := range order {
				check(!seen[q], "braid witness order repeats a position")
				for r, required := range pred[q] {
					check(!required || seen[r], "braid witness order violates dependency")
				}
				seen[q] = true
				out[p] = word[q]
			}
			check(reflect.DeepEqual(wordMatrix(out, a, true), wordMatrix(word, a, true)), "braid witness changes the element")
			return out, position, true
		}
	}
	return nil, 0, false
}

func finiteRoots(a diagram, omit int) ([]vector, map[string]bool) {
	found := make(map[string]bool)
	var queue []vector
	for j, e := range identity(len(a)) {
		if j != omit {
			queue = append(queue, e)
			found[fmt.Sprint(e)] = true
		}
	}
	for head := 0; head < len(queue); head++ {
		for s := range a {
			if s == omit {
				continue
			}
			y := reflectRoot(queue[head], s, a)
			key := fmt.Sprint(y)
			if !found[key] {
				found[key] = true
				queue = append(queue, y)
				check(len(queue) <= 240, "finite root closure exceeded 240 roots")
			}
		}
	}
	return queue, found
}

func negativeIndicator(v vector) int64 {
	if positive(v) {
		return 0
	}
	return 1
}

func dot(x, y vector) int64 {
	bounded(x)
	bounded(y)
	var out int64
	for j := range x {
		out += x[j] * y[j]
	}
	return out
}

func height(v vector) int64 {
	var out int64
	for _, x := range v {
		out += x
	}
	return out
}

func wordString(word []int) string {
	var out strings.Builder
	for _, s := range word {
		fmt.Fprint(&out, s)
	}
	return out.String()
}

func checkAffine() {
	a := en(9)
	delta := vector{2, 4, 6, 5, 4, 3, 2, 1, 3}
	gamma := vector{1, 1, 1, 1, 1, 0, 0, 0, 0}
	word := digits("312875645234123012856745231")
	c := wordMatrix(word, a, true)
	unit := identity(9)
	for j, col := range c {
		check(reflect.DeepEqual(apply(c, col), unit[j]), "affine C is not an involution")
	}
	check(reflect.DeepEqual(apply(c, delta), delta), "C does not fix delta")
	cg := apply(c, gamma)
	difference := make(vector, 9)
	for j := range difference {
		difference[j] = gamma[j] - cg[j]
	}
	d := make(vector, 9)
	var descents []int
	for j, e := range unit {
		d[j] = pairing(e, difference, a)
		if !positive(c[j]) {
			descents = append(descents, j)
		}
	}
	check(reflect.DeepEqual(d, vector{2, -1, 1, -1, 1, -1, 1, -1, -1}), "affine d differs: %v", d)
	check(reflect.DeepEqual(descents, []int{1, 3, 5, 7, 8}), "affine descents differ: %v", descents)
	for _, s := range descents {
		for _, t := range a[s] {
			image := make(vector, 9)
			for j := range image {
				image[j] = c[s][j] + c[t][j]
			}
			check(positive(image) && d[s]+d[t] >= 0, "affine edge sign fails: %d,%d", s, t)
		}
	}
	roots, rootSet := finiteRoots(a, 7)
	check(len(roots) == 240 && rootSet[fmt.Sprint(gamma)], "affine finite root system differs")
	type class struct{ intercept, slope int64 }
	classes := make(map[class]int64)
	for _, eta := range roots {
		img := apply(c, eta)
		h := img[7]
		eta1 := make(vector, 9)
		for j := range eta1 {
			eta1[j] = img[j] - h*delta[j]
		}
		check(rootSet[fmt.Sprint(eta1)], "affine root finite part missing: %v", eta1)
		intercept := negativeIndicator(eta1) - h - negativeIndicator(eta)
		slope := -dot(eta, d)
		classes[class{intercept, slope}]++
	}
	var keys []class
	var totalIntercept, totalSlope int64
	for key, multiplicity := range classes {
		keys = append(keys, key)
		check(key.intercept*key.slope >= 0, "affine inversion class changes sign: %v", key)
		if key.intercept > 0 {
			totalIntercept += key.intercept * multiplicity
		}
		if key.slope > 0 {
			totalSlope += key.slope * multiplicity
		}
	}
	check(totalIntercept == 27 && totalSlope == 92, "affine length totals differ: %d,%d", totalIntercept, totalSlope)
	sort.Slice(keys, func(i, j int) bool {
		if keys[i].intercept != keys[j].intercept {
			return keys[i].intercept < keys[j].intercept
		}
		return keys[i].slope < keys[j].slope
	})
	fmt.Println("Affine inversion classes (intercept, slope): multiplicity")
	for _, key := range keys {
		fmt.Printf("  (%d, %d): %d\n", key.intercept, key.slope, classes[key])
	}
	covers := make(map[string]bool)
	for i := range word {
		deleted := append([]int(nil), word[:i]...)
		deleted = append(deleted, word[i+1:]...)
		m := wordMatrix(deleted, a, false)
		rw := reduced(m, a)
		if len(rw) == 26 {
			_, _, found := braidWord(rw, a)
			check(found, "affine base cover has no braid witness: deletion %d", i)
			covers[fmt.Sprint(m)] = true
		}
	}
	check(len(covers) == 21, "affine distinct base-cover count: %d", len(covers))
	fmt.Println("All 21 base covers have explicit braid words")
	for _, s := range descents {
		left := make(matrix, 9)
		for j, col := range c {
			left[j] = reflectRoot(col, s, a)
		}
		rw := reduced(left, a)
		bw, p, found := braidWord(rw, a)
		check(found, "left descent s%d has no braid witness", s)
		fmt.Printf("s=%d left descent braid witness: %s[%s]%s\n", s, wordString(bw[:p]), wordString(bw[p:p+3]), wordString(bw[p+3:]))
	}
	r := identity(9)
	for j := range r {
		for i := range r[j] {
			r[j][i] += delta[i] * d[j]
		}
	}
	length := len(reduced(r, a))
	check(length == 92, "translation R length: %d", length)
	for j, col := range r {
		want := make(vector, 9)
		for i := range want {
			want[i] = c[j][i] + delta[i]*d[j]
		}
		check(reflect.DeepEqual(apply(c, col), want), "affine CR factorization differs at column %d", j)
	}
	var slope int64
	for _, eta := range roots {
		if b := -dot(eta, d); b > 0 {
			slope += b
		}
	}
	check(slope == 92 && dot(delta, d) == 0, "translation inversion slope or nilpotence differs")
	fmt.Println("Translation R length 92; affine inversion slope for R^k = 92")
}

func lowerCheck(v vector, word []int, a diagram, target int) {
	for j, s := range word {
		y := reflectRoot(v, s, a)
		check(positive(y) && height(y) < height(v), "lowering witness fails at position %d, generator %d", j, s)
		v = y
	}
	check(reflect.DeepEqual(v, identity(len(a))[target]), "lowering witness has wrong endpoint: %v", v)
}

func embedded(n int, initial vector, branch int64) vector {
	v := make(vector, n)
	copy(v, initial)
	v[n-1] = branch
	return v
}

func checkLowering() {
	sigma := digits("210432165432")
	for _, n := range []int{13, 10} {
		p, a := n-1, en(n)
		delta := embedded(n, vector{2, 4, 6, 5, 4, 3, 2, 1}, 3)
		gamma := embedded(n, vector{1, 2, 3, 2, 2, 1, 1, 0}, 1)
		tau := append([]int(nil), sigma...)
		tau = append(tau, p)
		tau = append(tau, digits("210321432543654765")...)
		tau = append(tau, p)
		tau = append(tau, digits("210321432")...)
		lowerCheck(gamma, sigma, a, p)
		sum := make(vector, n)
		for j := range sum {
			sum[j] = gamma[j] + delta[j]
		}
		lowerCheck(sum, tau, a, p)
	}
	lowerCheck(vector{3, 7, 10, 9, 7, 6, 4, 3, 1, 6}, digits("13547659210321432567892103214325436547659210321432"), en(10), 9)
	fmt.Println("All printed lowering witnesses correct")
}

// The helpers return new integers, never mutating their arguments.
func add(x, y *big.Int) *big.Int { return new(big.Int).Add(x, y) }
func sub(x, y *big.Int) *big.Int { return new(big.Int).Sub(x, y) }
func mul(x, y *big.Int) *big.Int { return new(big.Int).Mul(x, y) }
func neg(x *big.Int) *big.Int    { return new(big.Int).Neg(x) }

func toBig(v vector) bigVector {
	out := make(bigVector, len(v))
	for j, x := range v {
		out[j] = big.NewInt(x)
	}
	return out
}

func simplePairings(v bigVector, a diagram) bigVector {
	out := make(bigVector, len(v))
	for s, x := range v {
		out[s] = mul(big.NewInt(2), x)
		for _, t := range a[s] {
			out[s].Sub(out[s], v[t])
		}
	}
	return out
}

func bigPairing(x, y bigVector, a diagram) *big.Int {
	p := simplePairings(y, a)
	out := new(big.Int)
	for j := range x {
		out.Add(out, mul(x[j], p[j]))
	}
	return out
}

func sameBigVector(x, y bigVector) bool {
	if len(x) != len(y) {
		return false
	}
	for j := range x {
		if x[j].Cmp(y[j]) != 0 {
			return false
		}
	}
	return true
}

func checkSamples() {
	one := big.NewInt(1)
	for _, rankParameter := range []int64{3, 4, 5, 6, 10, 25, 100} {
		n := int(4*rankParameter + 1)
		a := en(n)
		r, coefficient := big.NewInt(rankParameter), big.NewInt(2*rankParameter-4)
		seedSmall := vector{rankParameter - 1, rankParameter, 2*rankParameter - 1}
		for j := 3; j < n-1; j++ {
			seedSmall = append(seedSmall, 2*rankParameter-int64(j/2))
		}
		seedSmall = append(seedSmall, rankParameter)
		seed := toBig(seedSmall)
		delta := toBig(embedded(n, vector{2, 4, 6, 5, 4, 3, 2, 1}, 3))
		gamma := toBig(embedded(n, vector{1, 2, 3, 2, 2, 1, 1, 0}, 1))
		check(bigPairing(seed, seed, a).Cmp(big.NewInt(2)) == 0, "uniform seed norm, r=%d", rankParameter)
		check(bigPairing(seed, delta, a).Cmp(neg(coefficient)) == 0, "uniform seed/delta pairing, r=%d", rankParameter)
		check(bigPairing(seed, gamma, a).Cmp(neg(r)) == 0, "uniform seed/gamma pairing, r=%d", rankParameter)
		for _, parameter := range []int64{0, 1, 2, 3, 25, 1000000} {
			k := big.NewInt(parameter)
			ak := mul(coefficient, k)
			t := add(mul(ak, k), mul(r, k))
			beta := make(bigVector, n)
			for j := range beta {
				beta[j] = add(sub(seed[j], mul(ak, gamma[j])), mul(t, delta[j]))
				check(beta[j].Sign() > 0, "uniform root not strictly positive: r=%d, k=%d, j=%d", rankParameter, parameter, j)
			}
			m := simplePairings(beta, a)
			expected := bigVector{big.NewInt(rankParameter - 2), big.NewInt(2 - rankParameter), neg(add(one, ak))}
			for j := 3; j < 8; j++ {
				x := add(one, ak)
				if j%2 == 0 {
					x = neg(x)
				}
				expected = append(expected, x)
			}
			expected = append(expected, neg(add(one, t)))
			for j := 9; j < n-1; j++ {
				x := int64(1)
				if j%2 == 0 {
					x = -1
				}
				expected = append(expected, big.NewInt(x))
			}
			expected = append(expected, add(one, ak))
			check(sameBigVector(m, expected), "uniform pairing formula: r=%d, k=%d", rankParameter, parameter)
			check(bigPairing(beta, beta, a).Cmp(big.NewInt(2)) == 0, "uniform root norm: r=%d, k=%d", rankParameter, parameter)
			for s := range a {
				for _, t := range a[s] {
					check(add(m[s], m[t]).Sign() <= 0, "uniform edge inequality: r=%d, k=%d, edge=%d,%d", rankParameter, parameter, s, t)
				}
			}
		}
	}
	fmt.Println("Uniform formula samples verified through r=100, k=1000000 (exact arithmetic; not a symbolic proof)")
	a := en(10)
	seed := toBig(vector{3, 7, 10, 9, 7, 6, 4, 3, 1, 6})
	delta := toBig(vector{2, 4, 6, 5, 4, 3, 2, 1, 0, 3})
	gamma := toBig(vector{1, 2, 3, 2, 2, 1, 1, 0, 0, 1})
	for _, parameter := range []int64{0, 1, 2, 25, 1000000} {
		k := big.NewInt(parameter)
		t := add(mul(k, k), mul(big.NewInt(3), k))
		beta := make(bigVector, 10)
		for j := range beta {
			beta[j] = add(sub(seed[j], mul(k, gamma[j])), mul(t, delta[j]))
		}
		u := add(one, k)
		expected := bigVector{big.NewInt(-1), one, neg(add(big.NewInt(2), k)), u, neg(u), u, neg(u), u, neg(add(one, t)), add(big.NewInt(2), k)}
		check(sameBigVector(simplePairings(beta, a), expected), "E10 pairing formula: k=%d", parameter)
		check(bigPairing(beta, beta, a).Cmp(big.NewInt(2)) == 0, "E10 root norm: k=%d", parameter)
	}
	fmt.Println("E10 pairing formula samples verified through k=1000000")
}

func main() {
	flag.Usage = func() {
		fmt.Fprintln(flag.CommandLine.Output(), "Usage: review-infinite\nIndependent affine, lowering-witness and exact family-sample checks; no input files.")
	}
	flag.Parse()
	check(flag.NArg() == 0, "unexpected arguments: %v", flag.Args())
	checkAffine()
	checkLowering()
	checkSamples()
	fmt.Println("All checks passed.")
}
