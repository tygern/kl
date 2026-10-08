package main

// The stable summary of the run, built from the same certificate fields as
// supplement/run_proofs.py and compared with expected-summary.json.

import (
	"fmt"
	"path/filepath"
)

func (r *runner) summary() *object {
	summary := newObject()
	finite := newObject()
	for _, n := range []int{6, 7, 8} {
		row := r.read(sprintf("research/en_e8/e%d-recursive.json", n))
		entry := pick(row, "group_order", "right_terminal_count", "commuting_terminals", "noncommuting_terminal_count", "fc_count")
		var lengths []int64
		for _, x := range items(field(row, "bad")) {
			lengths = append(lengths, integer(field(x, "length")))
		}
		entry.set("bad_lengths", sortedInt64(lengths))
		finite.set(sprintf("E%d", n), entry)
	}
	summary.set("finite", finite)
	d := r.read("research/en_independent/e8-d7-terminals.json")
	summary.set("E8_D7", pick(d, "cosets", "parabolic_right_terminals", "candidates_tested", "terminal_count"))
	d, err := readJSON(filepath.Join(r.logs, "D6-recurrence-certificate.stdout.log"))
	must(err)
	summary.set("D6", pick(d, "records", "root_polynomial", "root_right_descent_checks", "root_R_reciprocity", "evaluator_called"))
	d = r.read("research/en_uniform/referee-certificate.json")
	var ranks []any
	for _, row := range items(field(d, "rows")) {
		ranks = append(ranks, field(row, "r"))
	}
	summary.set("uniform", newObject().set("symbolic_assertions", field(field(d, "symbolic"), "all_assertions_pass")).set("checked_ranks", ranks))
	d = r.read("research/en_affine_referee/proved-affine-certificate.json")
	affine := newObject().set("length", []any{field(field(d, "length_formula"), "intercept"), field(field(d, "length_formula"), "slope")}).
		set("terminal_for_all_k", field(d, "terminal_for_all_k"))
	d = r.read("research/en_affine_referee/fc-cover-certificate.json")
	for _, k := range []string{"FC_count", "maximum_FC_length", "base_Bruhat_cover_count", "FC_base_covers"} {
		affine.set(k, field(d, k))
	}
	summary.set("affine", affine)
	d = r.read("research/en_affine_referee/e10-certificate.json")
	summary.set("E10", pick(d, "rank", "base_length", "T_minus_identity_cube_zero", "full_support_for_all_k", "maximum_independent_size", "eligible_mu"))
	summary.set("finite_matchings", r.read("research/ai-review-notes/finite-descents.json"))
	d = r.read("results/e6-mu-table.json")
	summary.set("E6_mu", pick(d, "group_order", "fully_commutative_count", "max_mu", "mu_histogram", "pairs_with_mu_at_least_2",
		"fully_commutative_lower_endpoint", "fully_commutative_upper_endpoint", "pairs_attaining_max_mu",
		"largest_mu_with_fully_commutative_upper_endpoint"))
	d = r.read("results/e9-odd-gap-certificate.json")
	oddGap := []any{}
	for _, e := range items(field(d, "elements")) {
		oddGap = append(oddGap, pick(e, "name", "word_string", "length", "L", "R", "involution", "terminal", "quotient_size", "lower_ideal_size", "eligible_fully_commutative_bottoms"))
	}
	summary.set("E9_odd_gap", oddGap)
	d = r.read("results/d8-gern-certificate.json")
	gern := newObject()
	ranksMap := mapping(field(d, "ranks"))
	for _, rank := range sortedKeys(ranksMap) {
		gern.set(rank, pick(ranksMap[rank], "w_length", "x_length", "lower_ideal_size", "interval_size", "P_x_w_ascending", "P_e_w_ascending", "mu_x_w"))
	}
	summary.set("D6_D8_Gern", gern)
	d = r.read("results/fc-maxima-certificate.json")
	maxima := newObject()
	ranksMap = mapping(field(d, "ranks"))
	for _, rank := range sortedKeys(ranksMap) {
		maxima.set(rank, pick(ranksMap[rank], "fully_commutative_count", "maximum_length"))
	}
	summary.set("FC_maxima", maxima)
	summary.set("terminal_structure", terminalStructureSummary(r.read("results/terminal-structure-certificate.json")))
	// Default rows, identical in both modes (checked above).
	d = canonical(filepath.Join(r.root, "payload/results/uniform-family-certificate.json"))
	summary.set("uniform_family", uniformFamilySummary(d))
	d = r.read("results/affine-d4-independent.json")
	summary.set("affine_D4", pick(d, "lengths", "lower_ideal_size", "polynomial", "mu"))
	return summary
}

// terminalStructureSummary extracts the stable parts of the structure
// certificate exactly as the original runner and packager do.
func terminalStructureSummary(d any) *object {
	out := newObject()
	rows := newObject()
	gern := mapping(field(d, "gern_plus_two"))
	for _, k := range sortedKeys(gern) {
		if v, ok := gern[k].(map[string]any); ok {
			rows.set(k, field(v, "rows_not_from_type_D"))
		}
	}
	out.set("rows_not_from_type_D", rows)
	typeD := newObject()
	enumeration := mapping(field(d, "type_D_enumeration"))
	for _, k := range sortedKeys(enumeration) {
		typeD.set(k, field(enumeration[k], "lengths"))
	}
	out.set("type_D_noncommuting_terminals", typeD)
	negated, layers := newObject(), newObject()
	structure := mapping(field(d, "structure"))
	for _, k := range sortedKeys(structure) {
		negated.set(k, field(structure[k], "negated_positive_roots"))
		layers.set(k, field(structure[k], "layers"))
	}
	out.set("negated_roots", negated)
	out.set("layers", layers)
	chain := []any{}
	for _, rel := range items(field(field(field(d, "right_weak_order_chain"), "E8"), "relations")) {
		// Python's `prefix and suffix` on two booleans.
		both := boolean(field(rel, "prefix")) && boolean(field(rel, "suffix"))
		chain = append(chain, []any{field(rel, "lower"), field(rel, "upper"), both})
	}
	out.set("chain_E8", chain)
	inner, w0 := newObject(), newObject()
	factors := mapping(field(d, "inner_longest_parabolic_factors"))
	for _, k := range sortedKeys(factors) {
		inner.set(k, field(factors[k], "max_inner_longest_parabolic_length"))
		if f, ok := mapping(factors[k])["w0_factorization"]; ok {
			w0.set(k, f)
		}
	}
	out.set("inner_longest_parabolic_lengths", inner)
	out.set("w0_factorizations", w0)
	return out
}

func uniformFamilySummary(d any) *object {
	lengths := []any{}
	for _, row := range items(field(d, "pairs")) {
		lengths = append(lengths, field(row, "length"))
	}
	covers := []any{}
	for _, c := range items(field(d, "bruhat_covers")) {
		covers = append(covers, []any{field(c, "r"), field(c, "k"), field(c, "distinct_covers"), field(c, "fully_commutative_covers")})
	}
	return newObject().set("checked_pairs", field(d, "checked_pairs")).set("lengths", lengths).
		set("length_formula_holds", field(d, "length_formula_holds_for_all_checked_pairs")).set("bruhat_covers", covers)
}

func sprintf(format string, args ...any) string { return fmt.Sprintf(format, args...) }
