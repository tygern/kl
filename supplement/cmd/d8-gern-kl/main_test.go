package main

import "testing"

func TestD6BothModels(t *testing.T) {
	a := runModel("sp", 6, 2)
	b := runModel("geo", 6, 3)
	if a.ideal != 3184 || a.interval != 1676 || a.lw != 15 || a.lx != 4 || a.mu != 1 {
		t.Fatalf("unexpected D6 sp data: %+v", a)
	}
	want := []int64{1, 6, 11, 6, 1, 1}
	if !eqVec(a.Px[:deg(&a.Px)+1], want) {
		t.Fatalf("wrong P: %v", a.Px)
	}
	if a.Px != b.Px || a.Pe != b.Pe || a.ideal != b.ideal || !eqVec(a.rankInterval, b.rankInterval) {
		t.Fatal("models disagree")
	}
	if a.mismatches != 0 || b.mismatches != 0 || a.polys != 32 {
		t.Fatalf("mismatches %d %d polys %d", a.mismatches, b.mismatches, a.polys)
	}
}
