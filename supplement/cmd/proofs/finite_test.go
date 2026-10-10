package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestProofMode(t *testing.T) {
	for _, tc := range []struct {
		full, finite bool
		want         string
	}{
		{false, false, "default"},
		{true, false, "full"},
		{false, true, "finite"},
		{true, true, ""},
	} {
		got, err := proofMode(tc.full, tc.finite)
		if got != tc.want || (err != nil) != (tc.full && tc.finite) {
			t.Fatalf("proofMode(%v, %v) = %q, %v", tc.full, tc.finite, got, err)
		}
	}
}

func finiteSummaryFixture() (map[string]any, any) {
	inputs := map[string]any{}
	for _, n := range []int{6, 7, 8} {
		inputs[sprintf("research/en_e8/e%d-recursive.json", n)] = literal(`{
			"group_order": 51840, "right_terminal_count": 72, "commuting_terminals": 22,
			"noncommuting_terminal_count": 1, "fc_count": 662, "bad": [{"length": 7}]}`)
	}
	inputs["research/en_independent/e8-d7-terminals.json"] = literal(`{
		"cosets": 2160, "parabolic_right_terminals": 280, "candidates_tested": 604800, "terminal_count": 64}`)
	inputs["research/ai-review-notes/finite-descents.json"] = literal(`[{"type":"E6", "word":"1325213"}]`)
	inputs["research/exceptional_referee/checks.json"] = literal(`{"status":"passed"}`)
	inputs["research/en_e8/e8-recursive-chains.json"] = literal(`{
		"group_order":696729600, "right_terminal_count":2160, "terminal_count":64,
		"commuting_terminals":58, "noncommuting_terminal_count":6,
		"right_terminal_sets_equal":true, "terminal_sets_equal":true,
		"chains":[{"chain":"e7", "added_generators":[7,2,3,1,0,4,5,6], "coset_counts":[2,3,4,8,10,27,56,240],
		"right_terminal_counts":[2,3,6,24,40,72,576,2160], "candidates_tested":[2,6,12,48,240,1080,4032,138240], "terminal_count":64}]}`)
	inputs["results/d6-ambient-certificate.json"] = literal(`{
		"polynomial":[1,6,11,6,1,1], "mu":1, "ideal_size":3184, "interval_size":1676,
		"interval_rank_vector":[1], "length_gap":11, "models":[{"model":"E7", "ideal_size":3184,
		"interval_size":1676, "P_xb":[1,6,11,6,1,1], "P_eb":[1,6,11,6,1,1], "mu":1, "b_length":15, "x_length":4}]}`)
	d6 := literal(`{"records":24245, "root_polynomial":[1,6,11,6,1,1],
		"root_right_descent_checks":4, "root_R_reciprocity":true, "evaluator_called":false}`)
	return inputs, d6
}

func TestFiniteSummaryReadsOnlyFiniteInputs(t *testing.T) {
	inputs, d6 := finiteSummaryFixture()
	readCounts := map[string]int{}
	summary := finiteProofSummary(func(path string) any {
		v, ok := inputs[path]
		if !ok {
			t.Fatalf("finite summary read an unrelated or unavailable input: %s", path)
		}
		readCounts[path]++
		return v
	}, d6)
	if equal, diff := sameJSON(summary.keys, finiteSummaryKeys); !equal {
		t.Fatalf("unexpected finite summary scope: %s", diff)
	}
	for path := range inputs {
		if readCounts[path] != 1 {
			t.Errorf("%s read %d times, want once", path, readCounts[path])
		}
	}
	core, err := normalize(finiteCoreSummary(func(path string) any { return inputs[path] }, d6))
	if err != nil {
		t.Fatal(err)
	}
	for _, key := range []string{"finite", "E8_D7", "D6"} {
		if equal, diff := sameJSON(summary.values[key], field(core, key)); !equal {
			t.Errorf("%s differs from the shared default/full summary: %s", key, diff)
		}
	}
}

func writeTestJSON(t *testing.T, path string, value any) {
	t.Helper()
	data, err := encodeJSON(value)
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, data, 0o644); err != nil {
		t.Fatal(err)
	}
}

func requireFailure(t *testing.T, contains string, action func()) {
	t.Helper()
	defer func() {
		rec := recover()
		f, ok := rec.(failure)
		if !ok || !strings.Contains(f.msg, contains) {
			t.Fatalf("got panic %#v, want failure containing %q", rec, contains)
		}
	}()
	action()
}

func TestFiniteExpectedSummaryScopeAndTampering(t *testing.T) {
	inputs, d6 := finiteSummaryFixture()
	summary := finiteProofSummary(func(path string) any { return inputs[path] }, d6)
	r := &runner{root: t.TempDir(), run: t.TempDir()}
	expected, err := normalize(summary)
	if err != nil {
		t.Fatal(err)
	}
	// Saved results outside the finite proof are deliberately different and
	// need not even have the schema required by the full runner.
	mapping(expected)["uniform"] = "not recomputed"
	mapping(expected)["affine"] = false
	path := filepath.Join(r.root, "expected-summary.json")
	writeTestJSON(t, path, expected)
	r.compareSummary(summary, finiteSummaryKeys)
	// Default/full comparison still requires the entire expected summary.
	requireFailure(t, "affine", func() { r.compareSummary(summary, nil) })
	saved, err := readJSON(filepath.Join(r.run, "summary.json"))
	if err != nil {
		t.Fatal(err)
	}
	if _, exists := mapping(saved)["affine"]; exists {
		t.Fatal("finite run reported an uncomputed affine outcome")
	}
	// A finite polynomial mismatch must fail; a subset comparison must not
	// silently accept incorrect or absent finite results.
	mapping(field(expected, "D6"))["root_polynomial"] = literal(`[1,6,11,6,1,2]`)
	writeTestJSON(t, path, expected)
	requireFailure(t, "D6.root_polynomial", func() { r.compareSummary(summary, finiteSummaryKeys) })
	delete(mapping(expected), "D6")
	writeTestJSON(t, path, expected)
	requireFailure(t, `missing key "D6"`, func() { r.compareSummary(summary, finiteSummaryKeys) })
}

func TestFiniteManuscriptHashCheckedBeforeBuilds(t *testing.T) {
	r := &runner{work: t.TempDir()}
	if err := os.Mkdir(filepath.Join(r.work, "results"), 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(r.work, "results/exceptional-leading.tex"), []byte("changed manuscript"), 0o644); err != nil {
		t.Fatal(err)
	}
	requireFailure(t, "manuscript differs", func() { r.verifyFinite("different recorded digest") })
	if len(r.steps) != 0 {
		t.Fatal("runner built programs before checking the manuscript hash")
	}
}
