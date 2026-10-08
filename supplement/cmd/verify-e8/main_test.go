package main

import "testing"

func TestRootsAndBasics(t *testing.T) {
	if len(roots()) != 240 {
		t.Fatal("E8 has 240 roots")
	}
	indep := 0
	for b := uint(0); b < 1<<n; b++ {
		if independent(b) {
			indep++
		}
	}
	if indep != 58 {
		t.Fatalf("independent subsets %d", indep)
	}
	for s := 0; s < n; s++ {
		if right(right(identity, s), s) != identity {
			t.Fatal("simple reflection is not an involution")
		}
		if desc(right(identity, s)) != 1<<uint(s) {
			t.Fatal("descent of a simple reflection")
		}
	}
	p := poincare([]int{1, 7, 11, 13, 17, 19, 23, 29})
	var sum int64
	for _, c := range p {
		sum += c
	}
	if sum != 696729600 {
		t.Fatalf("order of W(E8) = %d", sum)
	}
}
