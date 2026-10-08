package main

import "testing"

// The quotient engine must agree with the naive recursion on random elements;
// the pair counts are those recorded in the certificate (same mt19937 stream).
func TestValidateA4D4(t *testing.T) {
	if got := validate("A4", 1, 200, 10); got != 1270 {
		t.Fatalf("A4 pairs = %d, want 1270", got)
	}
	if got := validate("D4", 2, 200, 12); got != 1911 {
		t.Fatalf("D4 pairs = %d, want 1911", got)
	}
}

func TestMT19937(t *testing.T) {
	// std::mt19937 with default seed 5489: the 10000th output is 4123659995.
	m := newMT(5489)
	var v uint32
	for i := 0; i < 10000; i++ {
		v = m.next()
	}
	if v != 4123659995 {
		t.Fatalf("mt19937 10000th output = %d", v)
	}
}
