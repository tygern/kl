package main

import "testing"

func TestE6(t *testing.T) {
	e := &engine{T: 2}
	o := e.enumerate(6)
	if o.count != 662 || o.maxlen != 16 || len(o.longest) != 2 {
		t.Fatalf("E6: got %d/%d/%d", o.count, o.maxlen, len(o.longest))
	}
}
