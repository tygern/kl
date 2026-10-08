package main

import "testing"

func TestE6Certificate(t *testing.T) {
	c := e6Compute()
	if c.ElementCount != 51840 || c.FcCount != 662 || c.CommutingTerminals != 22 || len(c.NonFcWeakBad) != 1 {
		t.Fatalf("unexpected E6 certificate: %+v", c)
	}
}
