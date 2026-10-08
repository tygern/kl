package main

import "testing"

func TestGernAndWeyl(t *testing.T) {
	if got := len(gernWnWord(6)); got != (3*36+12)/8 {
		t.Fatalf("length of w_6 word = %d", got)
	}
	rs := en(6)
	if rs.length(rs.longestElement([]int{0, 1, 2, 3, 4, 5})) != 36 {
		t.Fatal("E6 longest element should have length 36")
	}
	if len(rs.positiveRoots()) != 36 {
		t.Fatal("E6 has 36 positive roots")
	}
	if got := len(enumerateBadDm(5)); got != 1 {
		t.Fatalf("D5 noncommuting terminals = %d, want 1", got)
	}
	if rankOf([][]int{{1, 2}, {2, 4}}) != 1 {
		t.Fatal("rank")
	}
}
