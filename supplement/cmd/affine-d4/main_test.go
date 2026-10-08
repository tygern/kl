package main

import "testing"

func TestEvaluate(t *testing.T) {
	x := evaluate(qx)
	if x.P != [4]int{-1, -2, -3, -4} || x.T != [4]int{0, 0, 1, 1} {
		t.Fatalf("bottom element wrong: %v", x)
	}
	w := evaluate(qw)
	if w.P != [4]int{1, 3, 2, 4} || w.T != [4]int{0, -1, 1, 0} {
		t.Fatalf("top element wrong: %v", w)
	}
}

func TestPolyArithmetic(t *testing.T) {
	// (1+q)(1+q) = 1+2q+q^2
	if !polyEqual(multiply(poly{1, 1}, poly{1, 1}), poly{1, 2, 1}) {
		t.Fatal("multiply")
	}
	if floorDiv(-1, 2) != -1 || floorDiv(8, 2) != 4 {
		t.Fatal("floorDiv")
	}
}
