package affine

import "testing"

func TestFiniteRootsAndTranslation(t *testing.T) {
	if len(FiniteRoots) != 240 || len(FiniteRootList) != 240 {
		t.Fatalf("finite roots: %d", len(FiniteRoots))
	}
	for _, e := range E {
		if Pair(Delta, Vec(e)) != 0 {
			t.Fatal("delta not null")
		}
	}
	if MM(Reflection(Gamma), Reflection(Add(Gamma, Delta, 1))) != Translate(Gamma, 1) {
		t.Fatal("translation identity")
	}
}

func TestReducedWordRoundTrip(t *testing.T) {
	w := []int{0, 1, 2, 3, 2, 8}
	a := WordMatrix(w, true)
	r := ReducedWord(a)
	if WordMatrix(r, true) != a || len(r) != len(w) {
		t.Fatalf("round trip %v", r)
	}
}
