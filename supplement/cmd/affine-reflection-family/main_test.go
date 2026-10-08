package main

import "testing"

func TestFamilyCertificate(t *testing.T) {
	out := build()
	if out.FiniteRootCount != 240 || out.LengthSlope != 58 || out.LengthOffsetForKPlus1 != -25 {
		t.Fatalf("symbolic counts wrong: %d %d %d", out.FiniteRootCount, out.LengthSlope, out.LengthOffsetForKPlus1)
	}
	if len(out.Checks) != 11 || out.Checks[0].Length != 33 || out.Checks[10].Length != 33+58*10 {
		t.Fatalf("lengths wrong")
	}
	if len(out.BaseReflectionWord) != 33 || len(out.TranslationReducedWord) != 58 {
		t.Fatalf("word lengths wrong: %d %d", len(out.BaseReflectionWord), len(out.TranslationReducedWord))
	}
	for _, r := range out.Checks {
		if !r.Terminal || len(r.Extensions) != 2 {
			t.Fatalf("k=%d not terminal or wrong extension count", r.K)
		}
	}
}

func TestGroupBasics(t *testing.T) {
	g := newCoxeter()
	if g.length(g.elt([]int{0, 1, 0})) != 3 || g.length(g.elt([]int{0, 0})) != 0 {
		t.Fatal("length")
	}
	if g.pairing(g.e[0], g.e[0]) != 2 || g.pairing(g.e[2], g.e[8]) != -1 {
		t.Fatal("Cartan form")
	}
}
