// Command exceptional-cosets prints exact root certificates for exceptional
// cosets carrying Gern's D6 interval in the Weyl groups E7 and E8.
//
// It ports research/exceptional_cosets.py to Go. Matrices are stored by
// columns in the basis of simple roots, in the Python original's layout. The
// numbering of E_n is 1--3--4--5--...--n with 2 attached to 4. The program
// checks Coxeter length and descent statements only; it makes no Kazhdan-
// Lusztig computation. All arithmetic is exact integer arithmetic (entries are
// tiny, far below the int range). The output is deterministic and is
// byte-identical to the Python original's standard output (the original prints
// text lines, with Python tuple/list reprs, not JSON); every assert of the
// original is an explicit check here that exits non-zero with a message.
//
// Invocation (from the repository root):
//
//	go run -C tools ./cmd/exceptional-cosets
//
// Imports: standard library only; no other package of this repository.
package main

import (
	"fmt"
	"os"
	"strconv"
	"strings"
)

// matrix is stored by columns: m[j] is the image of simple root alpha_j.
type matrix [][]int

func check(ok bool, format string, args ...any) {
	if !ok {
		fmt.Fprintf(os.Stderr, "assertion failed: "+format+"\n", args...)
		os.Exit(1)
	}
}

func isPositive(root []int) bool {
	for _, c := range root {
		if c < 0 {
			return false
		}
	}
	return true
}

func isNegative(root []int) bool {
	for _, c := range root {
		if c > 0 {
			return false
		}
	}
	return true
}

// tuple formats a root like a Python tuple repr.
func tuple(root []int) string {
	parts := make([]string, len(root))
	for i, c := range root {
		parts[i] = strconv.Itoa(c)
	}
	return "(" + strings.Join(parts, ", ") + ")"
}

// list formats ints like a Python list repr.
func list(v []int) string {
	parts := make([]string, len(v))
	for i, c := range v {
		parts[i] = strconv.Itoa(c)
	}
	return "[" + strings.Join(parts, ", ") + "]"
}

type edge struct{ s, t int }

func certificate(n int, prescribedStrip []int) {
	edges := []edge{{0, 2}, {2, 3}, {1, 3}}
	for i := 3; i < n-1; i++ {
		edges = append(edges, edge{i, i + 1})
	}
	cartan := make([][]int, n)
	for i := range cartan {
		cartan[i] = make([]int, n)
		cartan[i][i] = 2
	}
	for _, e := range edges {
		cartan[e.s][e.t] = -1
		cartan[e.t][e.s] = -1
	}
	identity := make(matrix, n)
	for j := range identity {
		identity[j] = make([]int, n)
		identity[j][j] = 1
	}

	// rightMultiply returns matrix * s_s: column j becomes
	// column_j - cartan[s][j] * column_s.
	rightMultiply := func(m matrix, s int) matrix {
		image := m[s]
		out := make(matrix, n)
		for j, column := range m {
			out[j] = make([]int, n)
			for i := 0; i < n; i++ {
				out[j][i] = column[i] - cartan[s][j]*image[i]
			}
		}
		return out
	}

	longest := func(nodes []int) (matrix, []int) {
		m, word := identity, []int{}
		for {
			s := -1
			for _, j := range nodes { // nodes ascending, so the first hit is the minimum
				if isPositive(m[j]) {
					s = j
					break
				}
			}
			if s < 0 {
				return m, word
			}
			m = rightMultiply(m, s)
			word = append(word, s)
			check(len(word) <= 1000, "longest-word search did not terminate")
		}
	}

	all := make([]int, n)
	for i := range all {
		all[i] = i
	}
	w0, fullWord := longest(all)
	nodes := []int{1, 2, 3, 4, 5, 6}
	_, parabolicWord := longest(nodes)
	check(len(parabolicWord) == 30, "len(parabolic_word) = %d", len(parabolicWord))
	a := w0
	for _, s := range parabolicWord {
		check(isNegative(a[s]), "a[%d] not negative while stripping parabolic word", s)
		a = rightMultiply(a, s)
	}
	lengthA := len(fullWord) - 30
	for _, j := range nodes {
		check(isPositive(a[j]), "a[%d] not positive", j)
	}
	remainder := a
	for _, j := range prescribedStrip {
		check(isNegative(remainder[j-1]), "(%d, %d, %s)", n, j, tuple(remainder[j-1]))
		remainder = rightMultiply(remainder, j-1)
	}
	var adjacent []edge
	for _, e := range edges {
		if isNegative(remainder[e.s]) && isNegative(remainder[e.t]) {
			adjacent = append(adjacent, edge{e.s + 1, e.t + 1})
		}
	}
	check(len(adjacent) > 0, "no adjacent right descents remain")
	fmt.Printf("E%d: length(a)=%d; lengths(a*x6,a*w6)=(%d,%d); interval rank=11\n",
		n, lengthA, lengthA+4, lengthA+15)
	fmt.Println("  stripped right descents:", list(prescribedStrip))
	pairs := make([]string, len(adjacent))
	for i, e := range adjacent {
		pairs[i] = fmt.Sprintf("(%d, %d)", e.s, e.t)
	}
	fmt.Println("  adjacent right descents remaining: [" + strings.Join(pairs, ", ") + "]")
	for _, e := range adjacent {
		fmt.Printf("  image(alpha_%d) = %s\n", e.s, tuple(remainder[e.s-1]))
		fmt.Printf("  image(alpha_%d) = %s\n", e.t, tuple(remainder[e.t-1]))
	}
}

func main() {
	certificate(7, []int{1, 3, 4, 2, 5, 4, 3, 1, 6, 5, 4, 2, 3, 4, 5})
	certificate(8, []int{1, 3, 4, 2, 5, 4, 3, 1, 6})
}
