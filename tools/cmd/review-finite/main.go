// Command review-finite independently checks the eleven finite terminal table
// words and the exceptional orthogonal-reflection descriptors reviewed on
// 10 October 2026. It uses only the standard library, no supplement code and
// no input files. It does not replace the exhaustive terminal enumerations.
//
// Run from the repository root:
//
//	go run -C tools ./cmd/review-finite
package main

import (
	"flag"
	"fmt"
	"math/bits"
	"os"
	"reflect"
)

type vector []int64

// A matrix is stored by columns in the simple-root basis.
type matrix []vector
type edge struct{ s, t int }
type tableRow struct {
	rank        int
	word, lower string
}

var rows = []tableRow{
	{6, "1325213", "135"},
	{7, "1326213", "136"},
	{7, "13256213", "1356"},
	{7, "132543621324356", "1356"},
	{7, "1325436210321432543621324356", "1356"},
	{8, "1327213", "137"},
	{8, "13257213", "1357"},
	{8, "61327213", "1367"},
	{8, "132543721324357", "1357"},
	{8, "1325437210321432543721324357", "1357"},
	{8, "75342312701234562103214325437210321432543721324357", "1357"},
}

var orthogonal = map[int][]vector{
	7: {{0, 0, 1, 1, 1, 1, 1}, {0, 1, 1, 1, 1, 0, 1},
		{0, 1, 2, 1, 0, 0, 1}, {1, 2, 2, 2, 1, 1, 1}},
	8: {{1, 2, 3, 3, 2, 2, 1, 2}, {1, 3, 4, 3, 2, 1, 0, 2}},
}

func check(ok bool, format string, args ...any) {
	if !ok {
		fmt.Fprintf(os.Stderr, "review-finite: "+format+"\n", args...)
		os.Exit(1)
	}
}

func edges(n int) []edge {
	var out []edge
	for s := 0; s < n-2; s++ {
		out = append(out, edge{s, s + 1})
	}
	return append(out, edge{2, n - 1})
}

func identity(n int) matrix {
	m := make(matrix, n)
	for j := range m {
		m[j] = make(vector, n)
		m[j][j] = 1
	}
	return m
}

func element(n int, word string) matrix {
	m := identity(n)
	for _, letter := range word {
		s := int(letter - '0')
		check(s >= 0 && s < n, "invalid generator %q in rank %d", letter, n)
		old := append(vector(nil), m[s]...)
		for i := range old {
			m[s][i] = -old[i]
		}
		for _, e := range edges(n) {
			t := -1
			if e.s == s {
				t = e.t
			} else if e.t == s {
				t = e.s
			}
			if t >= 0 {
				for i := range old {
					m[t][i] += old[i]
				}
			}
		}
	}
	return m
}

func apply(m matrix, v vector) vector {
	out := make(vector, len(v))
	for j, c := range v {
		for i, x := range m[j] {
			out[i] += c * x
		}
	}
	return out
}

// Root coordinates must have one sign and cannot all vanish.
func positive(v vector) bool {
	pos, neg := false, false
	for _, x := range v {
		pos = pos || x > 0
		neg = neg || x < 0
	}
	check(pos != neg, "vector is not a signed root: %v", v)
	return pos
}

func roots(n int) ([]vector, map[string]bool) {
	found := make(map[string]bool)
	queue := []vector(identity(n))
	for _, v := range queue {
		found[fmt.Sprint(v)] = true
	}
	for head := 0; head < len(queue); head++ {
		v := queue[head]
		for s := 0; s < n; s++ {
			b := append(vector(nil), v...)
			b[s] = -v[s]
			for _, e := range edges(n) {
				if e.s == s {
					b[s] += v[e.t]
				} else if e.t == s {
					b[s] += v[e.s]
				}
			}
			key := fmt.Sprint(b)
			if !found[key] {
				positive(b)
				found[key] = true
				queue = append(queue, b)
			}
		}
	}
	check(len(queue) == map[int]int{6: 72, 7: 126, 8: 240}[n],
		"E%d root count: %d", n, len(queue))
	return queue, found
}

func pairing(n int, u, v vector) int64 {
	var out int64
	for i := range u {
		out += 2 * u[i] * v[i]
	}
	for _, e := range edges(n) {
		out -= u[e.s]*v[e.t] + u[e.t]*v[e.s]
	}
	return out
}

func letterMask(word string) uint {
	var out uint
	for _, s := range word {
		out |= 1 << uint(s-'0')
	}
	return out
}

func main() {
	flag.Usage = func() {
		fmt.Fprintln(flag.CommandLine.Output(), "Usage: review-finite\nIndependent finite table and orthogonal-reflection checks; no input files.")
	}
	flag.Parse()
	check(flag.NArg() == 0, "unexpected arguments: %v", flag.Args())
	rootLists := make(map[int][]vector)
	rootSets := make(map[int]map[string]bool)
	for _, n := range []int{6, 7, 8} {
		rootLists[n], rootSets[n] = roots(n)
	}
	for _, row := range rows {
		n, word := row.rank, row.word
		m := element(n, word)
		reverse := []byte(word)
		for i, j := 0, len(reverse)-1; i < j; i, j = i+1, j-1 {
			reverse[i], reverse[j] = reverse[j], reverse[i]
		}
		inverse := element(n, string(reverse))
		check(reflect.DeepEqual(m, inverse), "E%d word %s is not an involution", n, word)
		length := 0
		for _, v := range rootLists[n] {
			if positive(v) && !positive(apply(m, v)) {
				length++
			}
		}
		check(length == len(word), "E%d word %s: length %d, printed %d", n, word, length, len(word))
		want := letterMask(row.lower)
		for _, action := range []matrix{m, inverse} {
			var descents uint
			for s, v := range action {
				if !positive(v) {
					descents |= 1 << uint(s)
				}
			}
			check(descents == want, "E%d word %s: descent mask %b, want %b", n, word, descents, want)
			for _, e := range edges(n) {
				if descents&(1<<uint(e.s)|1<<uint(e.t)) == 0 {
					continue
				}
				image := make(vector, n)
				for i := range image {
					image[i] = action[e.s][i] + action[e.t][i]
				}
				check(positive(image), "E%d word %s: terminality fails at edge %v", n, word, e)
			}
		}
		support, maximum := letterMask(word), 0
		for subset := uint(0); subset < 1<<uint(n); subset++ {
			if subset & ^support != 0 {
				continue
			}
			independent := true
			for _, e := range edges(n) {
				mask := uint(1)<<uint(e.s) | uint(1)<<uint(e.t)
				independent = independent && subset&mask != mask
			}
			if independent && bits.OnesCount(subset) > maximum {
				maximum = bits.OnesCount(subset)
			}
		}
		check(maximum == bits.OnesCount(want), "E%d word %s: maximum independent size %d", n, word, maximum)
		var descents []int
		for s := 0; s < n; s++ {
			if want&(1<<uint(s)) != 0 {
				descents = append(descents, s)
			}
		}
		fmt.Printf("E%d, length %d, descents %v: reduced, involution, terminal, maximum independent set\n", n, length, descents)
	}
	for _, n := range []int{7, 8} {
		vectors := orthogonal[n]
		for i, u := range vectors {
			check(rootSets[n][fmt.Sprint(u)], "E%d descriptor is not a root: %v", n, u)
			for j, v := range vectors {
				var want int64
				if i == j {
					want = 2
				}
				check(pairing(n, u, v) == want, "E%d descriptor pairing (%d,%d)", n, i, j)
			}
		}
		columns := identity(n)
		for _, v := range columns {
			for _, beta := range vectors {
				c := pairing(n, v, beta)
				for i := range v {
					v[i] -= c * beta[i]
				}
			}
		}
		matched := false
		for _, row := range rows {
			if row.rank == n && len(row.word) == map[int]int{7: 28, 8: 50}[n] {
				check(reflect.DeepEqual(columns, element(n, row.word)), "E%d orthogonal-reflection product differs", n)
				matched = true
			}
		}
		check(matched, "missing E%d exceptional table word", n)
		fmt.Printf("E%d exceptional terminal: orthogonal-reflection descriptor verified\n", n)
	}
	fmt.Println("All checks passed.")
}
