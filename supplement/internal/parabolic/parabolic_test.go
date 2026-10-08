package parabolic

import "testing"

func TestE6(t *testing.T) {
	g := New(6)
	par, order := g.ParabolicOneSided()
	if order != 1920 || len(par) != 40 {
		t.Fatalf("parabolic order %d, terminals %d", order, len(par))
	}
	if len(g.Cosets()) != 27 {
		t.Fatal("coset count")
	}
	fc := g.FullyCommutative()
	if len(fc) != 662 || g.FCStarChecks != 3620 {
		t.Fatalf("fc %d, star checks %d", len(fc), g.FCStarChecks)
	}
}
