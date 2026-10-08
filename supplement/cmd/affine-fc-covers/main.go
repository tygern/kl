// Command affine-fc-covers independently certifies the complete E9 FC
// catalogue (by FC-ascent closure, not by trusting its status field) and
// that none of the 21 Bruhat covers of the length-27 base element is fully
// commutative. It writes research/en_affine_referee/fc-cover-certificate.json
// and prints a summary (without the cover words) on stdout.
//
// It ports research/en_affine_referee/verify_fc_catalogue.py. Reads
// research/en_independent/e9-fc.json and the seed
// research/en_families/e9_full_support_eligible.json relative to the
// working directory.
//
// Imports: internal/affine only (shared with cmd/affine-proof, as the
// originals shared verify_families.py). It imports none of the other engines.
package main

import (
	"flag"
	"sort"

	q "github.com/tygern/kl/supplement/internal/affine"
)

const outPath = "research/en_affine_referee/fc-cover-certificate.json"

func main() {
	flag.Parse()
	q.Main(run)
}

func desc(a q.Mat) int {
	m := 0
	for s := 0; s < q.N; s++ {
		if !q.Positive(q.Col(a, s)) {
			m |= 1 << s
		}
	}
	return m
}

func independent(mask int) bool {
	for _, e := range q.Edges {
		if mask&(1<<e[0]) != 0 && mask&(1<<e[1]) != 0 {
			return false
		}
	}
	return true
}

func wordKey(w []int) string {
	b := make([]byte, len(w))
	for i, s := range w {
		b[i] = byte(s)
	}
	return string(b)
}

type entry struct{ a, ai q.Mat }

type certificate struct {
	Status                           string  `json:"status"`
	FCCatalogueIndependentlyComplete bool    `json:"FC_catalogue_independently_complete"`
	FCCount                          int     `json:"FC_count"`
	MaximumFCLength                  int     `json:"maximum_FC_length"`
	FCAscentTests                    int64   `json:"FC_ascent_tests"`
	MissingNoncommuting              int64   `json:"missing_with_noncommuting_descents"`
	MissingNonFCPredecessor          int64   `json:"missing_with_nonFC_predecessor"`
	BaseBruhatCoverCount             int     `json:"base_Bruhat_cover_count"`
	FCBaseCovers                     int     `json:"FC_base_covers"`
	BaseCoverReducedWords            [][]int `json:"base_cover_reduced_words"`
	Conclusion                       string  `json:"conclusion"`
}

type summary struct {
	Status                           string `json:"status"`
	FCCatalogueIndependentlyComplete bool   `json:"FC_catalogue_independently_complete"`
	FCCount                          int    `json:"FC_count"`
	MaximumFCLength                  int    `json:"maximum_FC_length"`
	FCAscentTests                    int64  `json:"FC_ascent_tests"`
	MissingNoncommuting              int64  `json:"missing_with_noncommuting_descents"`
	MissingNonFCPredecessor          int64  `json:"missing_with_nonFC_predecessor"`
	BaseBruhatCoverCount             int    `json:"base_Bruhat_cover_count"`
	FCBaseCovers                     int    `json:"FC_base_covers"`
	Conclusion                       string `json:"conclusion"`
}

func run() {
	data := q.LoadCatalogue()
	states := map[q.Mat]int{}
	order := []q.Mat{}
	byWord := map[string]entry{}
	histogram := map[int]int64{}
	previous := -1
	for _, row := range data.Elements {
		word := row.Word
		q.Assert(len(word) == row.Length && row.Length >= previous, "length order violated at %v", word)
		previous = row.Length
		var a, ai q.Mat
		if len(word) == 0 {
			a, ai = q.E, q.E
		} else {
			parent, ok := byWord[wordKey(word[:len(word)-1])]
			q.Assert(ok, "missing parent of %v", word)
			last := word[len(word)-1]
			q.Assert(q.Positive(q.Col(parent.a, last)), "nonreduced %v", word)
			a = q.Right(parent.a, last)
			ai = q.Left(parent.ai, last)
		}
		_, seen := states[a]
		q.Assert(!seen, "duplicate element %v", word)
		d := desc(a)
		q.Assert(d == row.Rmask && desc(ai) == row.Lmask, "descent masks differ at %v", word)
		q.Assert(independent(d), "dependent descents at %v", word)
		for s := 0; s < q.N; s++ {
			if d&(1<<s) != 0 {
				l, ok := states[q.Right(a, s)]
				q.Assert(ok && l == len(word)-1, "predecessor of %v not FC", word)
			}
		}
		states[a] = len(word)
		order = append(order, a)
		byWord[wordKey(word)] = entry{a, ai}
		histogram[len(word)]++
	}
	maxLen := 0
	for l := range histogram {
		if l > maxLen {
			maxLen = l
		}
	}
	q.Assert(len(states) == 44199 && maxLen == 44, "catalogue size %d max %d", len(states), maxLen)
	var ascents, missingNoncommuting, missingBad int64
	for _, a := range order {
		length := states[a]
		d := desc(a)
		for s := 0; s < q.N; s++ {
			if d&(1<<s) != 0 {
				continue
			}
			ascents++
			b := q.Right(a, s)
			if lb, ok := states[b]; ok {
				q.Assert(lb == length+1, "ascent length mismatch")
				continue
			}
			bd := desc(b)
			if !independent(bd) {
				missingNoncommuting++
				continue
			}
			found := false
			for t := 0; t < q.N; t++ {
				if bd&(1<<t) != 0 {
					if _, ok := states[q.Right(b, t)]; !ok {
						found = true
						break
					}
				}
			}
			q.Assert(found, "missing FC ascent")
			missingBad++
		}
	}
	// A shortest omitted FC element would have all shorter FC predecessors
	// present and would violate the last assertion. Therefore closure is complete.
	q.Assert(len(data.LengthDistribution) == 45, "length distribution size")
	for i := 0; i < 45; i++ {
		q.Assert(histogram[i] == data.LengthDistribution[i], "length distribution differs at %d", i)
	}
	q.Assert(ascents == data.AscentsTested, "ascents %d != %d", ascents, data.AscentsTested)
	word := q.BaseWord27()
	covers := map[q.Mat][]int{}
	for i := range word {
		sub := append(append([]int{}, word[:i]...), word[i+1:]...)
		a := q.WordMatrix(sub, false)
		rw := q.ReducedWord(a)
		if len(rw) == 26 {
			covers[a] = rw
		}
	}
	q.Assert(len(covers) == 21, "cover count %d", len(covers))
	words := make([][]int, 0, len(covers))
	for a, rw := range covers {
		_, inFC := states[a]
		q.Assert(!inFC, "a base cover is fully commutative")
		words = append(words, rw)
	}
	sort.Slice(words, func(i, j int) bool {
		for k := 0; k < len(words[i]) && k < len(words[j]); k++ {
			if words[i][k] != words[j][k] {
				return words[i][k] < words[j][k]
			}
		}
		return len(words[i]) < len(words[j])
	})
	const conclusion = "For every k>=0 and every fully commutative x, ordinary mu(x,b_k)=0 for the length 27+92k family"
	out := certificate{Status: "All assertions passed", FCCatalogueIndependentlyComplete: true,
		FCCount: len(states), MaximumFCLength: 44, FCAscentTests: ascents,
		MissingNoncommuting: missingNoncommuting, MissingNonFCPredecessor: missingBad,
		BaseBruhatCoverCount: len(covers), FCBaseCovers: 0, BaseCoverReducedWords: words, Conclusion: conclusion}
	q.WriteJSON(outPath, out)
	q.WriteJSON("", summary{out.Status, true, out.FCCount, 44, ascents, missingNoncommuting, missingBad,
		len(covers), 0, conclusion})
}
