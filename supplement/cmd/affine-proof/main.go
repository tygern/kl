// Command affine-proof checks the proved affine E8 family of odd-length
// terminal conjugates (length 27+92k for every k>=0) and writes
// research/en_affine_referee/proved-affine-certificate.json.
//
// It ports supplement/affine_proof.py (which uses the unchanged arithmetic
// of research/en_affine_referee/verify_families.py). Reads the seed
// research/en_families/e9_full_support_eligible.json and the FC catalogue
// research/en_independent/e9-fc.json relative to the working directory.
//
// Imports: internal/affine only (shared with cmd/affine-fc-covers, as the
// originals shared verify_families.py). It imports none of the other engines.
package main

import (
	"flag"

	q "github.com/tygern/kl/supplement/internal/affine"
)

const outPath = "research/en_affine_referee/proved-affine-certificate.json"

func main() {
	flag.Parse()
	q.Main(run)
}

func run() {
	for _, e := range q.E {
		q.Assert(q.Pair(q.Delta, q.Vec(e)) == 0, "delta not orthogonal to the simple roots")
	}
	q.Assert(q.FiniteRoots[q.Gamma], "gamma is not a finite root")
	q.Assert(q.MM(q.Reflection(q.Gamma), q.Reflection(q.Add(q.Gamma, q.Delta, 1))) == q.Translate(q.Gamma, 1),
		"t_gamma != r_gamma r_(gamma+delta)")
	word := q.BaseWord27()
	a := q.WordMatrix(word, true)
	q.Assert(q.MM(a, a) == q.E, "base element is not an involution")
	var d q.Vec
	shift := q.Add(q.Gamma, q.MV(a, q.Gamma), -1)
	for i, e := range q.E {
		d[i] = q.Pair(q.Vec(e), shift)
	}
	q.Assert(d == q.Vec{2, -1, 1, -1, 1, -1, 1, -1, -1}, "column slopes %v", d)
	q.Assert(q.MM(q.MM(q.Translate(q.Gamma, 1), a), q.Translate(q.Gamma, -1)) == q.Shifted(a, d, 1),
		"conjugation formula")
	result := q.AuditFamily(word, d, [2]int64{27, 92}, []int{1, 3, 5, 7, 8}, "proved odd-length affine conjugates")
	q.Assert(result.CompleteFCMaskCandidates == 1, "mask candidates %d", result.CompleteFCMaskCandidates)
	q.Assert(len(result.EligibleFCBottomsForAllK) >= 1 && intsEq(result.EligibleFCBottomsForAllK[0].Word, []int{1, 3, 5, 7, 8}),
		"eligible bottom is not {1,3,5,7,8}")
	q.WriteJSON(outPath, result)
	q.WriteJSON("", map[string]any{
		"affine_length":      []int{27, 92},
		"finite_roots":       len(q.FiniteRoots),
		"terminal_for_all_k": result.TerminalForAllK,
	})
}

func intsEq(a, b []int) bool {
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
