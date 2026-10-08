// Package d6 holds the primitives that the Python D6 recurrence-certificate
// verifier shares with the original KL evaluator: polynomial arithmetic from
// computations/coxeter.py (trim, add, mul) and the sparse signed-permutation
// Coxeter group of computations/sparse_kl.py (length, right multiplication,
// descents, Bruhat comparison, lower ideals, R-polynomials).
//
// Deliberately absent: the Kazhdan-Lusztig evaluator (kl, corrections,
// with_descent, verify_reciprocity, the cache of evaluated values and the
// certificate writer). The certificate checker in cmd/d6-certificate never
// calls an evaluator because this package does not contain one.
package d6

import "fmt"

// Poly is a polynomial in q with ascending powers; the zero polynomial is the
// empty slice. Arithmetic is exact int64 with overflow panics.
type Poly []int64

// Trim returns a copy of p without trailing zero coefficients (never nil).
func Trim(p []int64) Poly {
	n := len(p)
	for n > 0 && p[n-1] == 0 {
		n--
	}
	out := make(Poly, n)
	copy(out, p[:n])
	return out
}

func mulChecked(a, b int64) int64 {
	if a == 0 || b == 0 {
		return 0
	}
	c := a * b
	if c/b != a || (a == -1 && b == -1<<63) || (b == -1 && a == -1<<63) {
		panic(fmt.Sprintf("int64 overflow in %d*%d", a, b))
	}
	return c
}

func addChecked(a, b int64) int64 {
	c := a + b
	if (c > a) != (b > 0) {
		panic(fmt.Sprintf("int64 overflow in %d+%d", a, b))
	}
	return c
}

// Add returns p + scale*q^shift as a trimmed polynomial (coxeter.add).
// shift must be non-negative.
func Add(p, q []int64, shift int, scale int64) Poly {
	if shift < 0 {
		panic("negative shift")
	}
	n := len(p)
	if len(q)+shift > n {
		n = len(q) + shift
	}
	a := make([]int64, n)
	copy(a, p)
	for i, c := range q {
		a[i+shift] = addChecked(a[i+shift], mulChecked(scale, c))
	}
	return Trim(a)
}

// Mul returns the product p*q (coxeter.mul).
func Mul(p, q []int64) Poly {
	a := Poly{}
	for i, c := range p {
		a = Add(a, q, i, c)
	}
	return a
}

// Equal reports whether two coefficient lists are identical (no trimming).
func Equal(p, q []int64) bool {
	if len(p) != len(q) {
		return false
	}
	for i := range p {
		if p[i] != q[i] {
			return false
		}
	}
	return true
}
