package main

import (
	"fmt"
	"testing"
)

func TestUnfolding(t *testing.T) {
	r := runUnfolding()
	if r.IntervalSize != 448 || fmt.Sprint(r.RankVector) != "[1 5 14 32 63 98 109 80 36 9 1]" {
		t.Fatalf("unexpected interval data: %d %v", r.IntervalSize, r.RankVector)
	}
	if fmt.Sprint(r.DetectedRank3Edges) != "[[0 2] [1 2] [2 3] [2 4]]" {
		t.Fatalf("edges %v", r.DetectedRank3Edges)
	}
}

func TestGernD4D6(t *testing.T) {
	rows := runGernCoatoms(4, 6)
	if len(rows) != 2 {
		t.Fatalf("rows %d", len(rows))
	}
	if rows[0].Length != 7 || rows[0].Rank != 4 || rows[0].Coatoms != 6 || rows[0].AllLowerCovers != 6 {
		t.Fatalf("D4 row %+v", rows[0])
	}
	if rows[1].Length != 15 || rows[1].Rank != 11 || rows[1].Coatoms != 12 || rows[1].AllLowerCovers != 12 {
		t.Fatalf("D6 row %+v", rows[1])
	}
}

func TestAffineRelations(t *testing.T) {
	g := newAffineD(4, 3)
	if g.word([]int{0, 0}) != g.e || g.word([]int{2, 3, 2, 3, 2, 3}) != g.e {
		t.Fatal("relations")
	}
	if len(g.els) < 2 {
		t.Fatal("ball")
	}
}
