package main

import "testing"

func TestRootSystem(t *testing.T) {
	initRoots()
	if len(roots) != 240 {
		t.Fatalf("got %d roots, want 240", len(roots))
	}
	// s_0 applied to the identity element negates root 0 and moves node 1.
	e := identity()
	s0 := right(e, 0)
	if desc(s0) != 1 || desc(right(s0, 0)) != 0 {
		t.Fatal("right multiplication by a simple reflection is not an involution with the right descent")
	}
	if !lessWord([]int{1, 2}, []int{1, 3}) || lessWord([]int{1, 3}, []int{1, 2}) || !lessWord([]int{1}, []int{1, 0}) {
		t.Fatal("lessWord is not lexicographic")
	}
}
