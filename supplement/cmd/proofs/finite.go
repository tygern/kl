package main

import (
	"fmt"
	"path/filepath"
)

// finiteSummaryKeys is the exact part of the distributed expected summary
// recomputed by -finite. In particular, no saved infinite-family or full KL
// table result is reported as a successful computation in this mode.
var finiteSummaryKeys = []string{
	"finite", "E8_D7", "D6", "finite_matchings", "terminal_data_checks", "E8_chains", "D6_ambient",
}

func proofMode(full, finite bool) (string, error) {
	if full && finite {
		return "", fmt.Errorf("-finite and -full are mutually exclusive")
	}
	if finite {
		return "finite", nil
	}
	if full {
		return "full", nil
	}
	return "default", nil
}

// verifyFinite regenerates the finite proof inputs without building any
// affine or indefinite-family program. The historical E6/E7 full-group
// snapshots remain comparison inputs, as in default mode; the flat and
// recursive classifications and the independent E8/D7 search run afresh.
func (r *runner) verifyFinite(manuscriptSHA256AtBuild string) {
	// This mode has no uniform verifier to record the audited manuscript.
	// Check the isolated manuscript snapshot against INPUTS.json directly.
	if digest(filepath.Join(r.work, "results/exceptional-leading.tex")) != manuscriptSHA256AtBuild {
		panic(failure{"Finite verifier manuscript differs from the one recorded at build time."})
	}

	flat := r.build("terminals-flat")
	recursive := r.build("terminals-recursive")
	d7 := r.build("e8-d7")
	verifyExceptional := r.build("verify-exceptional")
	verifyOutputs := r.build("verify-outputs")
	finiteDescents := r.build("finite-descents")
	terminalDataCheck := r.build("terminal-data-check")
	d6Certificate := r.build("d6-certificate")
	d6Ambient := r.build("d6-ambient")

	for _, n := range []int{6, 7, 8} {
		flatFile := fmt.Sprintf("research/en_e8/e%d-validation.json", n)
		if n == 8 {
			flatFile = "research/en_e8/e8-terminals.json"
		}
		r.execute(fmt.Sprintf("E%d-flat", n), []string{flat, "-rank", fmt.Sprint(n)}, flatFile, "")
		r.execute(fmt.Sprintf("E%d-recursive", n), []string{recursive, "-rank", fmt.Sprint(n)}, fmt.Sprintf("research/en_e8/e%d-recursive.json", n), "")
	}
	r.execute("E8-D7", []string{d7}, "research/en_independent/e8-d7-terminals.json", "")
	r.execute("E8-recursive-chains", []string{recursive, "-rank", "8", "-chains-certificate", "research/en_e8/e8-recursive-chains.json"}, "", "")
	r.compareSnapshot("research/en_e8/e8-recursive-chains.json", nil)
	r.execute("finite-source-snapshot-audit", []string{verifyExceptional}, "", "")
	r.execute("finite-method-comparison", []string{verifyOutputs}, "", "")
	r.execute("finite-support-matchings", []string{finiteDescents}, "", "")
	r.compareSnapshot("research/ai-review-notes/finite-descents.json", nil)
	r.execute("terminal-data-check", []string{terminalDataCheck}, "research/exceptional_referee/checks.json", "")
	r.compareSnapshot("research/exceptional_referee/checks.json", nil)
	r.execute("D6-recurrence-certificate", []string{d6Certificate}, "", "")
	r.execute("D6-ambient-E7-E8", []string{d6Ambient, "-out", "results/d6-ambient-certificate.json"}, "", "")
	r.compareSnapshot("results/d6-ambient-certificate.json", nil)

	r.compareSummary(finiteProofSummary(r.read, r.readD6RecurrenceReport()), finiteSummaryKeys)
	r.verifyPrintedTable(verifyOutputs)
}

func (r *runner) readD6RecurrenceReport() any {
	d, err := readJSON(filepath.Join(r.logs, "D6-recurrence-certificate.stdout.log"))
	must(err)
	return d
}

// finiteCoreSummary is also used by default/full mode, so the finite
// outcomes cannot drift between the different modes' summary builders.
func finiteCoreSummary(read func(string) any, d6 any) *object {
	summary := newObject()
	finite := newObject()
	for _, n := range []int{6, 7, 8} {
		row := read(sprintf("research/en_e8/e%d-recursive.json", n))
		entry := pick(row, "group_order", "right_terminal_count", "commuting_terminals", "noncommuting_terminal_count", "fc_count")
		var lengths []int64
		for _, x := range items(field(row, "bad")) {
			lengths = append(lengths, integer(field(x, "length")))
		}
		entry.set("bad_lengths", sortedInt64(lengths))
		finite.set(sprintf("E%d", n), entry)
	}
	summary.set("finite", finite)
	d := read("research/en_independent/e8-d7-terminals.json")
	summary.set("E8_D7", pick(d, "cosets", "parabolic_right_terminals", "candidates_tested", "terminal_count"))
	summary.set("D6", pick(d6, "records", "root_polynomial", "root_right_descent_checks", "root_R_reciprocity", "evaluator_called"))
	return summary
}

func finiteProofSummary(read func(string) any, d6 any) *object {
	summary := finiteCoreSummary(read, d6)
	summary.set("finite_matchings", read("research/ai-review-notes/finite-descents.json"))
	summary.set("terminal_data_checks", read("research/exceptional_referee/checks.json"))
	summary.set("E8_chains", finiteChainsSummary(read("research/en_e8/e8-recursive-chains.json")))
	summary.set("D6_ambient", d6AmbientSummary(read("results/d6-ambient-certificate.json")))
	return summary
}
