package main

// Part 2 (port of research/verify_i6_subword.py): an independent verification
// by direct subwords. No lifting property, no order comparison, no code shared
// with shape.go beyond the element type.

type subwordCert struct {
	IntervalSize   int     `json:"interval_size"`
	Rank           int     `json:"rank"`
	Atoms          int     `json:"atoms"`
	AllCommonBound int     `json:"all_common_bounds"`
	Minima         [][]int `json:"minima"`
	RelativeRanks  []int   `json:"relative_ranks"`
	Method         string  `json:"method"`
	Status         string  `json:"status"`
}

var (
	swIdentity = elem{1, 2, 3, 4, 5, 6}
	swBottom   = elem{-1, -2, 4, 3, 6, 5}
	swTop      = elem{-1, -6, 3, -4, 5, -2}
	swAtomA    = elem{-4, -2, 1, 3, 6, 5}
	swAtomB    = elem{-1, -3, 4, 2, 6, 5}
	swExpected = []elem{{-4, -3, 1, 2, 6, 5}, {-2, -4, 3, 1, 6, 5}}
)

type subwordEngine struct {
	lengths map[elem]int
	words   map[elem][]int
	lowers  map[elem]map[elem]struct{}
}

// length is the type D length by signed inversion counts.
func (e *subwordEngine) length(v elem) int {
	if l, ok := e.lengths[v]; ok {
		return l
	}
	total := 0
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if v[i] > v[j] {
				total++
			}
		}
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if -v[i] > v[j] {
				total++
			}
		}
	}
	e.lengths[v] = total
	return total
}

// swRight is Gern's simple generators, with zero-based labels.
func swRight(v elem, s int) elem {
	u := v
	if s == 0 {
		u[0], u[1] = -u[1], -u[0]
	} else {
		u[s-1], u[s] = u[s], u[s-1]
	}
	return u
}

// word obtains one reduced word by removing a right descent.
func (e *subwordEngine) word(v elem) []int {
	if v == swIdentity {
		return nil
	}
	if w, ok := e.words[v]; ok {
		return w
	}
	for s := 0; s < n; s++ {
		z := swRight(v, s)
		if e.length(z) < e.length(v) {
			w := append(append([]int(nil), e.word(z)...), s)
			e.words[v] = w
			return w
		}
	}
	fail("no right descent for %v", v)
	return nil
}

// lower generates the lower ideal by all subwords, without order comparisons.
func (e *subwordEngine) lower(v elem) map[elem]struct{} {
	if l, ok := e.lowers[v]; ok {
		return l
	}
	elements := map[elem]struct{}{swIdentity: {}}
	for _, s := range e.word(v) {
		added := make([]elem, 0, len(elements))
		for z := range elements {
			added = append(added, swRight(z, s))
		}
		for _, z := range added {
			elements[z] = struct{}{}
		}
	}
	e.lowers[v] = elements
	return elements
}

func inLower(e *subwordEngine, v, of elem) bool {
	_, ok := e.lower(of)[v]
	return ok
}

func computeSubword() subwordCert {
	e := &subwordEngine{
		lengths: map[elem]int{},
		words:   map[elem][]int{},
		lowers:  map[elem]map[elem]struct{}{},
	}
	var interval []elem
	for z := range e.lower(swTop) {
		if inLower(e, swBottom, z) {
			interval = append(interval, z)
		}
	}
	var atoms []elem
	for _, z := range interval {
		if e.length(z) == e.length(swBottom)+1 {
			atoms = append(atoms, z)
		}
	}
	var common []elem
	for _, z := range interval {
		if inLower(e, swAtomA, z) && inLower(e, swAtomB, z) {
			common = append(common, z)
		}
	}
	var minimal []elem
	for _, z := range common {
		ok := true
		for _, t := range common {
			if t != z && inLower(e, t, z) {
				ok = false
				break
			}
		}
		if ok {
			minimal = append(minimal, z)
		}
	}
	sortLex(minimal)
	rank := e.length(swTop) - e.length(swBottom)
	if len(interval) != 1676 {
		fail("interval size %d, expected 1676", len(interval))
	}
	if rank != 11 {
		fail("rank %d, expected 11", rank)
	}
	if len(atoms) != 12 {
		fail("%d atoms, expected 12", len(atoms))
	}
	hasA, hasB := false, false
	for _, a := range atoms {
		hasA = hasA || a == swAtomA
		hasB = hasB || a == swAtomB
	}
	if !hasA || !hasB {
		fail("the two named atoms are not both atoms")
	}
	if len(common) != 720 {
		fail("%d common bounds, expected 720", len(common))
	}
	// set(minimal) == EXPECTED_MINIMA: minimal is duplicate-free, so compare
	// the sorted lists.
	sortedExpected := append([]elem(nil), swExpected...)
	sortLex(sortedExpected)
	if len(minimal) != len(sortedExpected) {
		fail("%d minimal common bounds, expected %d", len(minimal), len(sortedExpected))
	}
	for i := range minimal {
		if minimal[i] != sortedExpected[i] {
			fail("minimal common bounds %v differ from the expected set", minimal)
		}
	}
	ranks := make([]int, len(minimal))
	for i, z := range minimal {
		ranks[i] = e.length(z) - e.length(swBottom)
	}
	return subwordCert{
		IntervalSize:   len(interval),
		Rank:           rank,
		Atoms:          len(atoms),
		AllCommonBound: len(common),
		Minima:         toLists(minimal),
		RelativeRanks:  ranks,
		Method:         "Independent direct subword membership; no lifting-property comparisons or KL computation",
		Status:         "All assertions passed",
	}
}
