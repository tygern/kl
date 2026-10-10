package main

import "testing"

func TestWordRoundTrip(t *testing.T) {
	// Check the rset/rget round trip.
	var x uint64
	x = rset(x, 3, 77)
	if rget(x, 3) != 77 || rget(x, 2) != 0 {
		t.Fatal("rset/rget mismatch")
	}
}
