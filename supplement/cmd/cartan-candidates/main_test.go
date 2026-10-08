package main

import "testing"

func TestE10MaxOnly(t *testing.T) {
	res := search(10, 2, true, nil)
	if res.IndependenceNumber != 5 || res.IndependentSetsTested != 9 || res.PairingsTested != 1557 ||
		res.PositiveIntegerVectors != 1522 || res.NormTwoVectors != 8 || len(res.RealTerminalRoots) != 8 {
		t.Fatalf("unexpected E10 m2 max1 counts: %+v", res)
	}
	r := res.RealTerminalRoots[0]
	if r.Length != 137 || r.Height != 83 || r.RootReduction.SimpleRoot != 9 {
		t.Fatalf("unexpected first root: %+v", r)
	}
}

func TestE10Unrestricted(t *testing.T) {
	res := search(10, 1, false, nil)
	if res.IndependentSetsTested != 151 || res.PairingsTested != 3272 || res.PositiveIntegerVectors != 3263 ||
		res.NormTwoVectors != 1 || len(res.RealTerminalRoots) != 1 {
		t.Fatalf("unexpected E10 m1 max0 counts: %+v", res)
	}
}

func TestGroupBasics(t *testing.T) {
	g := en(6)
	w := g.elt([]int{0, 1, 2, 3, 2, 1, 0})
	if got := g.word(w); len(got) != 7 {
		t.Fatalf("word length %d", len(got))
	}
	if !equalElt(g.inv(g.inv(w)), w) {
		t.Fatal("inverse is not an involution")
	}
}
