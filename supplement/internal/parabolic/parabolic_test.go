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

func TestParabolicOrder(t *testing.T) {
	g := New(8)
	for _, c := range []struct {
		mask uint32
		want int64
	}{
		{0, 1},
		{1 << 7, 2},
		{1<<7 | 1<<2, 6},                         // A2
		{1<<0 | 1<<1 | 1<<7, 12},                 // A2 x A1
		{1<<7 | 1<<2 | 1<<3 | 1<<1, 192},         // D4
		{1<<7 | 1<<2 | 1<<3 | 1<<1 | 1<<0, 1920}, // D5
		{255 &^ (1 << 6), 2903040},               // E7
		{255 &^ 1, 322560},                       // D7
		{255 &^ (1 << 7), 40320},                 // A7
		{255, 696729600},                         // E8
	} {
		if got := g.ParabolicOrder(c.mask); got != c.want {
			t.Fatalf("order of mask %b: got %d, want %d", c.mask, got, c.want)
		}
	}
	if got := New(6).ParabolicOrder(63); got != 51840 {
		t.Fatalf("E6 order %d", got)
	}
}
