// Command d6-certificate verifies the saved D6 Kazhdan-Lusztig recurrence
// certificate (results/d6-recurrence-certificate.json) without invoking any
// KL computation routine.
//
// For every supplied KL dependency it checks the standard recurrence from the
// polynomial records, recalculating Bruhat comparisons, the entire correction
// index set, degree/positivity constraints and all required references.
// Separately it checks q^d P(q^-1) = sum R_xz P_zw at the root. The JSON
// summary is printed to stdout.
//
// Ports computations/verify_certificate.py. Imports only
// github.com/tygern/kl/supplement/internal/d6 (the primitives that the original
// verifier shared with sparse_kl.py and coxeter.py: polynomial add/mul/trim,
// signed-permutation length, right action, descents, Bruhat order, lower
// ideals, R-polynomials). That package contains no KL evaluator, so the
// original's "evaluator never ran" assertion holds by construction; it imports
// none of the other engines.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"sort"

	"github.com/tygern/kl/supplement/internal/d6"
)

type record struct {
	X               *int       `json:"x"`
	W               *int       `json:"w"`
	Polynomial      *[]int64   `json:"polynomial"`
	RightDescent    *int       `json:"right_descent"`
	CorrectionTerms *[][]int64 `json:"correction_terms"`
}

type rootRecurrence struct {
	RightDescent    *int       `json:"right_descent"`
	CorrectionTerms *[][]int64 `json:"correction_terms"`
}

type certificate struct {
	Group struct {
		Type string `json:"type"`
		Rank int    `json:"rank"`
	} `json:"group"`
	Elements        [][]int64        `json:"elements"`
	Root            []int            `json:"root"`
	Records         []record         `json:"records"`
	RootRecurrences []rootRecurrence `json:"root_recurrences"`
	NumberRecords   *int             `json:"number_records"`
}

type result struct {
	Records                int     `json:"records"`
	ComparableRecords      int     `json:"comparable_records"`
	StrictRecurrenceChecks int     `json:"strict_recurrence_checks"`
	IncomparableZeroChecks int     `json:"incomparable_zero_checks"`
	RootPolynomial         []int64 `json:"root_polynomial"`
	RootRightDescentChecks int     `json:"root_right_descent_checks"`
	RootRReciprocity       bool    `json:"root_R_reciprocity"`
	EvaluatorCalled        bool    `json:"evaluator_called"`
}

type check struct{ msg string }

func assert(ok bool, format string, args ...any) {
	if !ok {
		panic(check{fmt.Sprintf(format, args...)})
	}
}

// term is a correction term (z, exponent, mu).
type term struct {
	z    d6.Perm
	e    int64
	mu   int64
	zIdx int
}

func sortTerms(t []term) {
	sort.SliceStable(t, func(i, j int) bool {
		if t[i].z != t[j].z {
			return t[i].z.Less(t[j].z)
		}
		if t[i].e != t[j].e {
			return t[i].e < t[j].e
		}
		return t[i].mu < t[j].mu
	})
}

func termsEqual(a, b []term) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i].z != b[i].z || a[i].e != b[i].e || a[i].mu != b[i].mu {
			return false
		}
	}
	return true
}

type pairKey struct{ x, w d6.Perm }

type row struct {
	poly    []int64
	descent *int
	terms   *[][]int64
}

func verify(filename string) result {
	raw, err := os.ReadFile(filename)
	assert(err == nil, "cannot read certificate: %v", err)
	var data certificate
	assert(json.Unmarshal(raw, &data) == nil, "certificate is not valid JSON of the expected shape")

	g, err := d6.NewSparseCoxeter(data.Group.Type, data.Group.Rank)
	assert(err == nil, "%v", err)
	elements := make([]d6.Perm, len(data.Elements))
	seen := map[d6.Perm]bool{}
	for i, w := range data.Elements {
		p, err := d6.NewPerm(w)
		assert(err == nil, "%v", err)
		elements[i] = p
		seen[p] = true
	}
	assert(len(seen) == len(elements), "duplicate elements")
	for i, w := range elements {
		assert(w.N == elements[0].N && int(w.N) >= g.MinLetters(), "element %d has the wrong number of letters", i)
		abs := make([]int, w.N)
		negatives := 0
		for k := 0; k < int(w.N); k++ {
			a := int(w.V[k])
			if a < 0 {
				negatives++
				a = -a
			}
			abs[k] = a
		}
		sort.Ints(abs)
		for k, a := range abs {
			assert(a == k+1, "element %d is not a signed permutation", i)
		}
		assert(negatives%2 == 0, "element %d has an odd number of negative entries", i)
	}

	index := func(i int) d6.Perm {
		assert(i >= 0 && i < len(elements), "element index %d out of range", i)
		return elements[i]
	}

	records := map[pairKey]row{}
	var order []pairKey
	for _, r := range data.Records {
		assert(r.X != nil && r.W != nil && r.Polynomial != nil, "record missing x, w or polynomial")
		k := pairKey{index(*r.X), index(*r.W)}
		if _, ok := records[k]; !ok {
			order = append(order, k)
		}
		records[k] = row{*r.Polynomial, r.RightDescent, r.CorrectionTerms}
	}
	assert(data.NumberRecords != nil && len(records) == len(data.Records) && len(data.Records) == *data.NumberRecords,
		"record count mismatch")

	p := func(x, w d6.Perm) []int64 {
		r, ok := records[pairKey{x, w}]
		assert(ok, "missing KL record for a required reference")
		return r.poly
	}

	listedTerms := func(raw *[][]int64) []term {
		assert(raw != nil, "missing correction_terms")
		out := make([]term, 0, len(*raw))
		for _, t := range *raw {
			assert(len(t) == 3, "malformed correction term")
			assert(t[0] >= 0 && int(t[0]) < len(elements), "correction term index out of range")
			out = append(out, term{z: elements[t[0]], e: t[1], mu: t[2]})
		}
		sortTerms(out)
		return out
	}

	recurrence := func(x, w d6.Perm, s int) (d6.Poly, []term) {
		assert(g.HasDescent(w, s), "right_descent %d is not a descent of the top element", s)
		v, xs := g.Right(w, s), g.Right(x, s)
		var c int
		if g.Length(xs) < g.Length(x) {
			c = 1
		}
		value := d6.Add(nil, p(xs, v), 1-c, 1)
		value = d6.Add(value, p(x, v), c, 1)
		terms := []term{}
		for z := range g.Lower(v) {
			delta := g.Length(v) - g.Length(z)
			if delta <= 0 || delta%2 == 0 || !g.HasDescent(z, s) {
				continue
			}
			polynomial := p(z, v)
			// delta is positive and odd here, so (delta-1)/2 is exact.
			idx := (delta - 1) / 2
			var mu int64
			if len(polynomial) > idx {
				mu = polynomial[idx]
			}
			if mu != 0 && g.Leq(x, z) {
				exponent := (delta + 1) / 2
				terms = append(terms, term{z: z, e: int64(exponent), mu: mu})
				value = d6.Add(value, p(x, z), exponent, -mu)
			}
		}
		sortTerms(terms)
		return value, terms
	}

	comparable, strict, zero := 0, 0, 0
	for _, k := range order {
		x, w := k.x, k.w
		r := records[k]
		polynomial := r.poly
		switch {
		case x == w:
			assert(d6.Equal(polynomial, []int64{1}), "P(x,x) != 1")
			comparable++
		case !g.Leq(x, w):
			assert(len(polynomial) == 0, "nonzero polynomial for incomparable pair")
			zero++
		default:
			comparable++
			strict++
			assert(len(polynomial) > 0 && polynomial[0] == 1, "constant term is not 1")
			for _, c := range polynomial {
				assert(c >= 0, "negative coefficient")
			}
			assert(2*(len(polynomial)-1) < g.Length(w)-g.Length(x), "degree bound violated")
			assert(r.descent != nil, "missing right_descent")
			value, terms := recurrence(x, w, *r.descent)
			listed := listedTerms(r.terms)
			assert(termsEqual(terms, listed), "correction terms differ from the recomputed index set")
			assert(d6.Equal(polynomial, value), "recurrence value differs from the record")
		}
	}

	assert(len(data.Root) == 2, "root must be a pair of element indices")
	x, w := index(data.Root[0]), index(data.Root[1])
	rootDescents := map[int]bool{}
	for _, rr := range data.RootRecurrences {
		assert(rr.RightDescent != nil, "missing right_descent in root recurrence")
		value, terms := recurrence(x, w, *rr.RightDescent)
		assert(d6.Equal(value, p(x, w)), "root recurrence value differs from the root polynomial")
		assert(termsEqual(terms, listedTerms(rr.CorrectionTerms)), "root correction terms differ")
		rootDescents[*rr.RightDescent] = true
	}
	wd := g.Descents(w)
	assert(len(rootDescents) == len(wd), "root recurrences do not cover exactly the right descents of the top")
	for _, s := range wd {
		assert(rootDescents[s], "root recurrences miss descent %d", s)
	}

	polynomial, delta := p(x, w), g.Length(w)-g.Length(x)
	lhsRaw := make([]int64, max(delta+1, 0))
	for i := 0; i <= delta; i++ {
		if j := delta - i; j >= 0 && j < len(polynomial) {
			lhsRaw[i] = polynomial[j]
		}
	}
	lhs := d6.Trim(lhsRaw)
	rhs := d6.Poly{}
	for z := range g.Lower(w) {
		if g.Leq(x, z) {
			rhs = d6.Add(rhs, d6.Mul(g.RPoly(x, z), p(z, w)), 0, 1)
		}
	}
	assert(d6.Equal(lhs, rhs), "R-polynomial reciprocity identity fails at the root")
	// The verifier's KL evaluator must never have run: internal/d6 contains
	// none (see its TestNoEvaluator), so evaluator_called is false by
	// construction.
	if polynomial == nil {
		polynomial = []int64{}
	}
	return result{
		Records: len(records), ComparableRecords: comparable,
		StrictRecurrenceChecks: strict, IncomparableZeroChecks: zero,
		RootPolynomial: polynomial, RootRightDescentChecks: len(wd),
		RootRReciprocity: true, EvaluatorCalled: false,
	}
}

func main() {
	cert := flag.String("cert", "results/d6-recurrence-certificate.json", "certificate to verify (relative to the working directory)")
	flag.Parse()
	if flag.NArg() != 0 {
		fmt.Fprintln(os.Stderr, "d6-certificate: unexpected positional arguments")
		os.Exit(2)
	}
	var res result
	func() {
		defer func() {
			if r := recover(); r != nil {
				if c, ok := r.(check); ok {
					fmt.Fprintf(os.Stderr, "d6-certificate: VERIFICATION FAILED: %s\n", c.msg)
				} else {
					fmt.Fprintf(os.Stderr, "d6-certificate: VERIFICATION FAILED: %v\n", r)
				}
				os.Exit(1)
			}
		}()
		res = verify(*cert)
	}()
	out, err := json.MarshalIndent(res, "", "  ")
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	fmt.Println(string(out))
}
