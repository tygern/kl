package main

import (
	"io"
	"os"
	"strings"
	"testing"
)

func TestRightMultiplyInvolution(t *testing.T) {
	// Sign predicates and formatting helpers.
	if !isPositive([]int{0, 1}) || !isNegative([]int{0, -1}) || !isNegative([]int{0, 0}) {
		t.Fatal("sign predicates")
	}
	if tuple([]int{-1, 0, 2}) != "(-1, 0, 2)" || list([]int{1, 3}) != "[1, 3]" {
		t.Fatal("formatting")
	}
}

func TestCertificateOutput(t *testing.T) {
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatal(err)
	}
	old := os.Stdout
	os.Stdout = w
	certificate(7, []int{1, 3, 4, 2, 5, 4, 3, 1, 6, 5, 4, 2, 3, 4, 5})
	certificate(8, []int{1, 3, 4, 2, 5, 4, 3, 1, 6})
	os.Stdout = old
	w.Close()
	out, _ := io.ReadAll(r)
	lines := strings.Split(string(out), "\n")
	want := map[int]string{
		0: "E7: length(a)=33; lengths(a*x6,a*w6)=(37,48); interval rank=11",
		2: "  adjacent right descents remaining: [(6, 7)]",
		5: "E8: length(a)=90; lengths(a*x6,a*w6)=(94,105); interval rank=11",
	}
	for i, l := range want {
		if lines[i] != l {
			t.Fatalf("line %d: got %q want %q", i, lines[i], l)
		}
	}
}
