package main

import (
	"testing"

	q "github.com/tygern/kl/supplement/internal/affine"
)

func TestAffineProofWithoutCatalogue(t *testing.T) {
	// The package test directory has neither the research seed nor an FC
	// catalogue: the proof is built entirely from its printed finite data.
	result := prove()
	if len(result.LeftDescentBraidWitnesses) != 5 || !result.NoFCCoversForAllK ||
		!result.DescentPrefixesForAllK || !result.RightTranslationPowersForAllK {
		t.Fatal("incomplete cover certificate")
	}
	if result.MaximumIndependentSize != 5 || result.RightTranslationLengthFormula.Intercept != 0 ||
		result.RightTranslationLengthFormula.Slope != 92 {
		t.Fatal("incorrect maximum-descent or translation certificate")
	}
	// Supplement the all-parameter identities with direct checks that actual
	// concatenated descent words remain reduced and have the required matrix.
	for _, k := range []int{0, 1, 3} {
		ck := q.Shifted(result.BaseRowMatrix, result.ColumnSlopes, int64(k))
		for _, witness := range result.LeftDescentBraidWitnesses {
			w := append([]int{}, witness.Word...)
			for j := 0; j < k; j++ {
				w = append(w, result.RightTranslationReducedWord...)
			}
			if len(w) != 26+92*k || q.WordMatrix(w, true) != q.Left(ck, witness.Generator) {
				t.Fatalf("descent prefix failed at k=%d, s=%d", k, witness.Generator)
			}
		}
	}
}

func TestRejectMalformedBraidWitnesses(t *testing.T) {
	result := prove()
	tests := []struct {
		name   string
		mutate func([]q.DescentBraidWitness) []q.DescentBraidWitness
	}{
		{"missing descent", func(w []q.DescentBraidWitness) []q.DescentBraidWitness { return w[:4] }},
		{"wrong generator", func(w []q.DescentBraidWitness) []q.DescentBraidWitness { w[0].Generator = 3; return w }},
		{"wrong braid position", func(w []q.DescentBraidWitness) []q.DescentBraidWitness { w[0].BraidStart = 0; return w }},
		{"out-of-range braid", func(w []q.DescentBraidWitness) []q.DescentBraidWitness { w[0].BraidStart = 25; return w }},
		{"changed product", func(w []q.DescentBraidWitness) []q.DescentBraidWitness { w[0].Word[0] = 0; return w }},
	}
	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			witnesses := append([]q.DescentBraidWitness{}, result.LeftDescentBraidWitnesses...)
			for i := range witnesses {
				witnesses[i].Word = append([]int{}, witnesses[i].Word...)
			}
			defer func() {
				r := recover()
				if _, ok := r.(q.AssertionError); !ok {
					t.Fatalf("expected a checked rejection, got %v", r)
				}
			}()
			q.VerifyDescentBraids(result.BaseRowMatrix, len(result.BaseWord), result.Descents, test.mutate(witnesses))
		})
	}
}
