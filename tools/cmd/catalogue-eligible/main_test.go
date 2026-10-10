package main

import (
	"encoding/json"
	"io"
	"os"
	"path/filepath"
	"testing"
)

func TestGroupBasics(t *testing.T) {
	g := newEn(6)
	if got := g.length(g.elt([]int{0, 1, 0})); got != 3 {
		t.Fatalf("length of s0 s1 s0 = %d", got)
	}
	if g.fc(g.elt([]int{0, 1, 0})) {
		t.Fatal("s0 s1 s0 must not be fully commutative")
	}
	if !g.fc(g.elt([]int{0, 1, 2})) {
		t.Fatal("s0 s1 s2 must be fully commutative")
	}
	if !g.leq(g.elt([]int{0}), g.elt([]int{0, 1, 0})) || g.leq(g.elt([]int{2}), g.elt([]int{0, 1, 0})) {
		t.Fatal("Bruhat order wrong")
	}
	w := g.elt([]int{0, 1, 2, 3, 2, 4})
	if !g.equal(g.elt(g.word(w)), w) || !g.equal(g.inv(g.inv(w)), w) {
		t.Fatal("word / inverse round trip")
	}
	if floorDiv2(-3) != -2 || floorDiv2(3) != 1 {
		t.Fatal("floor division")
	}
}

func TestRunSynthetic(t *testing.T) {
	g := newEn(4) // type A4 path 0-1-2-3
	// The program does not re-derive full support, so a small terminal word
	// (s0 s2, commuting) flagged full_support exercises the whole computation.
	word := []int{0, 2}
	w := g.elt(word)
	if !g.terminal(w) || !g.fc(w) {
		t.Fatal("s0 s2 must be a terminal fully commutative element")
	}
	dir := t.TempDir()
	rows := map[string]any{"terminals": []any{
		map[string]any{"rank": 4, "word": word, "length": len(word), "support": []int{0, 1, 2, 3}, "R": g.desc(w), "full_support": true, "trial": 7},
		map[string]any{"rank": 4, "word": []int{0}, "length": 1, "support": []int{0}, "R": []int{0}, "full_support": false, "trial": 1},
	}}
	cat := map[string]any{"elements": []any{
		map[string]any{"word": word, "length": len(word), "Rmask": g.descMask(w), "Lmask": g.descMask(g.inv(w))},
	}}
	write := func(name string, v any) string {
		p := filepath.Join(dir, name)
		b, _ := json.Marshal(v)
		if err := os.WriteFile(p, b, 0o644); err != nil {
			t.Fatal(err)
		}
		return p
	}
	out := filepath.Join(dir, "out.json")
	run(4, write("c.json", rows), write("f.json", cat), out, io.Discard)
	var got []outRow
	b, _ := os.ReadFile(out)
	if err := json.Unmarshal(b, &got); err != nil || len(got) != 1 {
		t.Fatalf("output: %v %s", err, b)
	}
	if got[0].MaskCandidates != 1 || len(got[0].Eligible) != 1 || got[0].Eligible[0].Rank != 0 {
		t.Fatalf("unexpected row %+v", got[0])
	}
}
