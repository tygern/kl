package main

import "testing"

func TestD6ParabolicPair(t *testing.T) {
	d6 := newE("E7", 7).restrict("D6", []int{1, 2, 3, 4, 5, 6})
	rec, _ := record(d6, 7, "132543621324356", "1356")
	if rec.IdealSize != 3184 || rec.IntervalSize != 1676 || rec.Mu != 1 || len(rec.PXB) != 6 {
		t.Fatalf("unexpected record: %+v", rec)
	}
}

func TestGernSignedPermutations(t *testing.T) {
	w := gernWordToPerm(gernW6Word)
	x := gernWordToPerm(gernX6Word)
	if !sameInts(intsOf(w), gernW6Expected) || !sameInts(intsOf(x), gernX6Expected) {
		t.Fatalf("signed permutations %v %v", w, x)
	}
	if gernLen(w) != 15 || gernLen(x) != 4 {
		t.Fatalf("lengths %d %d", gernLen(w), gernLen(x))
	}
}
