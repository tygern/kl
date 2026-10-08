package indefinite

import "testing"

func TestE10Basics(t *testing.T) {
	c := Configure(10)
	if len(c.Edges) != 9 || len(c.Adj[2]) != 3 {
		t.Fatal("bad diagram")
	}
	beta := Vec{3, 7, 10, 9, 7, 6, 4, 3, 1, 6}
	if c.Pair(beta, beta) != 2 {
		t.Fatal("beta0 is not a real root")
	}
	w := c.RealRootWitness(beta)
	if w.Height != 56 {
		t.Fatalf("height %d", w.Height)
	}
	if !MatEq(c.MM(c.Reflection(beta), c.Reflection(beta)), c.E) {
		t.Fatal("reflection not an involution")
	}
}
