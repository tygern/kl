package main

import "testing"

func TestRows(t *testing.T) {
	for _, r := range rows {
		res := check(r)
		if res.IndependenceNumber != len(r.lower) {
			t.Fatalf("bad independence number for %s", r.word)
		}
	}
}
