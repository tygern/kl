package main

import "testing"

func TestSymbolic(t *testing.T) {
	if s := symbolic(); !s.AllPass {
		t.Fatal("symbolic assertions did not pass")
	}
}

func TestCheckRank3(t *testing.T) {
	row, w := checkRank(3)
	if row.Rank != 13 || row.Alpha != 7 || len(row.SampleChecks) != 5 || len(w) != 2 {
		t.Fatalf("unexpected rank-3 row: %+v", row)
	}
	if row.SampleChecks[1].Height != 162 || row.SampleChecks[4].Coordinate0 != 4000004000002 {
		t.Fatalf("unexpected samples: %+v", row.SampleChecks)
	}
}

func TestPolyRing(t *testing.T) {
	r := variable(0)
	if !r.mul(r).sub(r.mul(r)).eq(poly{}) || !cst(2).mul(frac(1, 2)).eq(cst(1)) {
		t.Fatal("polynomial ring arithmetic broken")
	}
}
