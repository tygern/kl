package main

import "testing"

func TestShape(t *testing.T) {
	c := computeShape()
	if c.IntervalSize != 1676 || c.Rank != 11 || c.AtomCount != 12 {
		t.Fatalf("unexpected shape: %d %d %d", c.IntervalSize, c.Rank, c.AtomCount)
	}
	if c.NonLatticeWitness.AllCommonUpperBoundsChecked != 720 {
		t.Fatalf("common bounds %d", c.NonLatticeWitness.AllCommonUpperBoundsChecked)
	}
	if len(c.NonLatticeWitness.MinimalCommonUpperBounds) != 2 {
		t.Fatal("expected two minimal common upper bounds")
	}
}

func TestSubword(t *testing.T) {
	c := computeSubword()
	if c.IntervalSize != 1676 || c.Atoms != 12 || c.AllCommonBound != 720 {
		t.Fatalf("unexpected subword result: %+v", c)
	}
	if len(c.RelativeRanks) != 2 || c.RelativeRanks[0] != 2 || c.RelativeRanks[1] != 3 {
		t.Fatalf("relative ranks %v", c.RelativeRanks)
	}
}

func TestLessLex(t *testing.T) {
	if !lessLex(elem{-4, 1, 1, 1, 1, 1}, elem{-3, 0, 0, 0, 0, 0}) || lessLex(elem{1, 2, 3, 4, 5, 6}, elem{1, 2, 3, 4, 5, 6}) {
		t.Fatal("lexicographic order wrong")
	}
}
