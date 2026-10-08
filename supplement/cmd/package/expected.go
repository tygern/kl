package main

// expected computes the stable mathematical outcomes from the shipped
// certificates in the repository; it must agree with the summary built by
// cmd/proofs. It is the Go counterpart of expected() in
// supplement/build_package.py, including its assertions of the values stated
// in the manuscript, so that a changed certificate cannot pass unnoticed.

import (
	"fmt"
	"sort"
)

func expected(read func(string) any) *object {
	result := newObject()
	finite := newObject()
	for _, n := range []int{6, 7, 8} {
		d := read(fmt.Sprintf("research/en_e8/e%d-recursive.json", n))
		entry := pick(d, "group_order", "right_terminal_count", "commuting_terminals", "noncommuting_terminal_count", "fc_count")
		var lengths []int64
		for _, row := range items(field(d, "bad")) {
			lengths = append(lengths, integer(field(row, "length")))
		}
		entry.set("bad_lengths", sortedInt64(lengths))
		finite.set(fmt.Sprintf("E%d", n), entry)
	}
	result.set("finite", finite)
	d := read("research/en_independent/e8-d7-terminals.json")
	result.set("E8_D7", pick(d, "cosets", "parabolic_right_terminals", "candidates_tested", "terminal_count"))
	result.set("D6", newObject().set("records", 24245).set("root_polynomial", []int64{1, 6, 11, 6, 1, 1}).
		set("root_right_descent_checks", 4).set("root_R_reciprocity", true).set("evaluator_called", false))
	ranks := make([]int64, 0, 38)
	for r := int64(3); r <= 40; r++ {
		ranks = append(ranks, r)
	}
	result.set("uniform", newObject().set("symbolic_assertions", true).set("checked_ranks", ranks))
	affine := newObject().set("length", []int64{27, 92}).set("terminal_for_all_k", true)
	d = read("research/en_affine_referee/fc-cover-certificate.json")
	for _, k := range []string{"FC_count", "maximum_FC_length", "base_Bruhat_cover_count", "FC_base_covers"} {
		affine.set(k, field(d, k))
	}
	result.set("affine", affine)
	d = read("research/en_affine_referee/e10-certificate.json")
	result.set("E10", pick(d, "rank", "base_length", "T_minus_identity_cube_zero", "full_support_for_all_k", "maximum_independent_size", "eligible_mu"))
	result.set("finite_matchings", read("research/ai-review-notes/finite-descents.json"))
	d = read("results/e6-mu-table.json")
	e6mu := pick(d, "group_order", "fully_commutative_count", "max_mu", "mu_histogram", "pairs_with_mu_at_least_2",
		"fully_commutative_lower_endpoint", "fully_commutative_upper_endpoint", "pairs_attaining_max_mu",
		"largest_mu_with_fully_commutative_upper_endpoint")
	result.set("E6_mu", e6mu)
	d = read("results/e9-odd-gap-certificate.json")
	oddGap := []any{}
	for _, e := range items(field(d, "elements")) {
		oddGap = append(oddGap, pick(e, "name", "word_string", "length", "L", "R", "involution", "terminal", "quotient_size", "lower_ideal_size", "eligible_fully_commutative_bottoms"))
	}
	result.set("E9_odd_gap", oddGap)
	d = read("results/d8-gern-certificate.json")
	gern := newObject()
	ranksMap := mapping(field(d, "ranks"))
	for _, rank := range sortedKeys(ranksMap) {
		gern.set(rank, pick(ranksMap[rank], "w_length", "x_length", "lower_ideal_size", "interval_size", "P_x_w_ascending", "P_e_w_ascending", "mu_x_w"))
	}
	result.set("D6_D8_Gern", gern)
	d = read("results/fc-maxima-certificate.json")
	maxima := newObject()
	ranksMap = mapping(field(d, "ranks"))
	for _, rank := range sortedKeys(ranksMap) {
		maxima.set(rank, pick(ranksMap[rank], "fully_commutative_count", "maximum_length"))
	}
	result.set("FC_maxima", maxima)
	result.set("terminal_structure", terminalStructureSummary(read("results/terminal-structure-certificate.json")))
	uniform := read("results/uniform-family-certificate.json")
	result.set("uniform_family", uniformFamilySummary(uniform))
	d = read("results/affine-d4-independent.json")
	result.set("affine_D4", pick(d, "lengths", "lower_ideal_size", "polynomial", "mu"))

	// Values stated in the manuscript, asserted here so that a changed
	// certificate cannot pass unnoticed.
	assertEqual(e6mu.values["max_mu"], literal(`10`), "E6 max_mu")
	assertEqual(field(e6mu.values["mu_histogram"], "10"), literal(`8`), "E6 mu_histogram[10]")
	assertEqual(e6mu.values["fully_commutative_lower_endpoint"], literal(`{"nonzero_mu_pairs":6431,"max_mu":1}`), "E6 fully_commutative_lower_endpoint")
	quotients, ideals := []any{}, []any{}
	for _, e := range oddGap {
		quotients = append(quotients, e.(*object).values["quotient_size"])
		ideals = append(ideals, e.(*object).values["lower_ideal_size"])
	}
	assertEqual(quotients, literal(`[364156,216990,207866]`), "E9 quotient sizes")
	assertEqual(ideals, literal(`[5826496,6943680,6651712]`), "E9 lower ideal sizes")
	assertEqual(gern.values["D8"].(*object).values["P_x_w_ascending"], literal(`[1,12,59,154,233,221,147,70,20,2]`), "D8 P_x_w")
	assertEqual(gern.values["D6"].(*object).values["P_x_w_ascending"], literal(`[1,6,11,6,1,1]`), "D6 P_x_w")
	maxLengths := map[string]any{}
	for _, k := range []string{"E10", "E11", "E12", "E13"} {
		if v, ok := maxima.values[k]; ok {
			maxLengths[k] = v.(*object).values["maximum_length"]
		}
	}
	assertEqual(maxLengths, literal(`{"E10":55,"E11":66,"E12":78,"E13":92}`), "FC maximum lengths E10-E13")
	affineD4 := result.values["affine_D4"].(*object)
	assertEqual(affineD4.values["polynomial"], literal(`[1,3,2]`), "affine D4 polynomial")
	assertEqual(affineD4.values["mu"], literal(`2`), "affine D4 mu")
	for _, cert := range []string{"results/e6-mu-table.json", "results/e9-odd-gap-certificate.json", "results/d8-gern-certificate.json", "results/fc-maxima-certificate.json",
		"results/terminal-structure-certificate.json", "results/uniform-family-certificate.json"} {
		assertEqual(field(read(cert), "status"), "passed", cert+" status")
	}
	if len(items(field(uniform, "pairs"))) != 24 {
		panic(failure{"ship the default-mode uniform certificate"})
	}
	return result
}

func assertEqual(got, want any, what string) {
	if equal, diff := sameJSON(got, want); !equal {
		panic(failure{"assertion failed: " + what + ": " + diff})
	}
}

// terminalStructureSummary and uniformFamilySummary select the stable parts
// of two certificates exactly as the original packager and runner do.
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

func sortedInt64(values []int64) []int64 {
	out := append([]int64(nil), values...)
	sort.Slice(out, func(i, j int) bool { return out[i] < out[j] })
	return out
}
