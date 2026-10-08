package main

import "testing"

func TestSmallPairAndCovers(t *testing.T) {
	failures = failures[:0]
	row := checkPair(3, 0)
	if row.Length != 75 || row.Rank != 13 || len(row.Ir) != 7 {
		t.Fatalf("unexpected row: length %d rank %d |I|=%d", row.Length, row.Rank, len(row.Ir))
	}
	cov := checkCovers(3, 0)
	if cov.DistinctCovers != 74 || cov.FullyCommutativeCov != 0 {
		t.Fatalf("unexpected covers: %+v", cov)
	}
	if len(failures) != 0 {
		t.Fatalf("expectations failed: %v", failures)
	}
}
