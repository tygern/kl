// Command affine-proof checks the proved affine E8 family of odd-length
// terminal conjugates (length 27+92k for every k>=0) and writes
// research/en_affine_referee/proved-affine-certificate.json.
//
// The displayed base word and five reduced braid witnesses are checked
// directly. The translation length formula propagates their braids to every
// parameter, without reading a seed file or the FC catalogue. The separate
// affine-fc-covers command retains catalogue closure and all 21 base covers
// as supplementary checks.
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
	result := prove()
	q.WriteJSON(outPath, result)
	q.WriteJSON("", map[string]any{
		"affine_length":            []int{27, 92},
		"finite_roots":             len(q.FiniteRoots),
		"terminal_for_all_k":       result.TerminalForAllK,
		"base_descent_braid_count": len(result.LeftDescentBraidWitnesses),
		"no_FC_covers_for_all_k":   result.NoFCCoversForAllK,
	})
}

func word(s string) []int {
	w := make([]int, len(s))
	for i := range s {
		q.Assert(s[i] >= '0' && s[i] < '0'+q.N, "invalid word digit")
		w[i] = int(s[i] - '0')
	}
	return w
}

func prove() q.Family {
	for _, e := range q.E {
		q.Assert(q.Pair(q.Delta, q.Vec(e)) == 0, "delta not orthogonal to the simple roots")
	}
	q.Assert(q.FiniteRoots[q.Gamma], "gamma is not a finite root")
	q.Assert(q.MM(q.Reflection(q.Gamma), q.Reflection(q.Add(q.Gamma, q.Delta, 1))) == q.Translate(q.Gamma, 1),
		"t_gamma != r_gamma r_(gamma+delta)")
	baseWord := word("312875645234123012856745231")
	a := q.WordMatrix(baseWord, true)
	q.Assert(q.MM(a, a) == q.E, "base element is not an involution")
	var d q.Vec
	shift := q.Add(q.Gamma, q.MV(a, q.Gamma), -1)
	for i, e := range q.E {
		d[i] = q.Pair(q.Vec(e), shift)
	}
	q.Assert(d == q.Vec{2, -1, 1, -1, 1, -1, 1, -1, -1}, "column slopes %v", d)
	q.Assert(q.MM(q.MM(q.Translate(q.Gamma, 1), a), q.Translate(q.Gamma, -1)) == q.Shifted(a, d, 1),
		"conjugation formula")
	witnesses := []q.DescentBraidWitness{
		// Concatenated strings make the verified three-letter braids visible.
		{Generator: 1, Word: word("3" + "282" + "7564534123012856745231"), BraidStart: 1},
		{Generator: 3, Word: word("1" + "282" + "7564534123012856745231"), BraidStart: 1},
		{Generator: 5, Word: word("31" + "282" + "764534123012856745231"), BraidStart: 2},
		{Generator: 7, Word: word("31" + "282" + "564534123012856745231"), BraidStart: 2},
		{Generator: 8, Word: word("3" + "121" + "8756453423012856745231"), BraidStart: 1},
	}
	return q.AuditFamily(baseWord, d, [2]int64{27, 92}, []int{1, 3, 5, 7, 8}, witnesses, "proved odd-length affine conjugates")
}
