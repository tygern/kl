package main

import (
	"bufio"
	"io"
	"testing"
)

// The reference certificate e6-fc.json records 662 FC elements for rank 6.
func TestE6Count(t *testing.T) {
	e := newEngine(6)
	e.run(bufio.NewWriter(io.Discard))
	if len(e.els) != 662 {
		t.Fatalf("E6 FC count = %d, want 662", len(e.els))
	}
}
