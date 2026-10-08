package main

import "testing"

func TestVerifySmall(t *testing.T) {
	row := verify(3, []int64{0, 1, 2, 10})
	if row.Rank != 13 || len(row.I) != 7 || row.A != 2 {
		t.Fatalf("unexpected row for r=3: rank %d, |I| %d, a %d", row.Rank, len(row.I), row.A)
	}
	if len(row.Checks) != 4 || !row.Checks[3].Terminal {
		t.Fatalf("unexpected checks: %+v", row.Checks)
	}
	if len(seed(3)) != 13 {
		t.Fatalf("seed length %d", len(seed(3)))
	}
}
