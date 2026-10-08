package main

import "sort"

// Part 1 (port of research/check_i6_shape.py with SparseCoxeter of type D):
// Bruhat order by the lifting property.

var (
	shapeX = elem{-1, -2, 4, 3, 6, 5}
	shapeW = elem{-1, -6, 3, -4, 5, -2}
)

type shapeCert struct {
	Bottom                 []int          `json:"bottom"`
	Top                    []int          `json:"top"`
	IntervalSize           int            `json:"interval_size"`
	Rank                   int            `json:"rank"`
	Atoms                  [][]int        `json:"atoms"`
	AtomCount              int            `json:"atom_count"`
	PrincipalLowerObstruct string         `json:"principal_lower_obstruction"`
	NonLatticeWitness      witnessPayload `json:"non_lattice_witness"`
	Checks                 string         `json:"checks"`
}

type witnessPayload struct {
	Atoms                       [][]int `json:"atoms"`
	MinimalCommonUpperBounds    [][]int `json:"minimal_common_upper_bounds"`
	RelativeRanks               []int   `json:"relative_ranks"`
	AllCommonUpperBoundsChecked int     `json:"all_common_upper_bounds_checked"`
}

// lifting holds the memoized structures of SparseCoxeter('D', 6).
type lifting struct {
	lengths map[elem]int
	lowers  map[elem]map[elem]struct{}
	leqMemo map[[2]elem]bool
}

func newLifting() *lifting {
	return &lifting{
		lengths: map[elem]int{},
		lowers:  map[elem]map[elem]struct{}{},
		leqMemo: map[[2]elem]bool{},
	}
}

// length is the type D length: inversions plus negative-sum pairs.
func (g *lifting) length(w elem) int {
	if l, ok := g.lengths[w]; ok {
		return l
	}
	l := 0
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if w[i] > w[j] {
				l++
			}
			if -w[i] > w[j] {
				l++
			}
		}
	}
	g.lengths[w] = l
	return l
}

// right is right multiplication by the simple generator s (0 is the D-type
// special generator; s >= 1 swaps positions s-1 and s).
func (g *lifting) right(w elem, s int) elem {
	a := w
	if s > 0 {
		a[s-1], a[s] = a[s], a[s-1]
	} else {
		a[0], a[1] = -a[1], -a[0]
	}
	return a
}

// firstDescent is the smallest s with l(ws) < l(w); w must not be the identity.
func (g *lifting) firstDescent(w elem) int {
	lw := g.length(w)
	for s := 0; s < n; s++ {
		if g.length(g.right(w, s)) < lw {
			return s
		}
	}
	fail("element %v has no right descent but nonzero length", w)
	return -1
}

func (g *lifting) leq(x, w elem) bool {
	if x == w {
		return true
	}
	if g.length(x) >= g.length(w) {
		return false
	}
	key := [2]elem{x, w}
	if r, ok := g.leqMemo[key]; ok {
		return r
	}
	s := g.firstDescent(w)
	xs, ws := g.right(x, s), g.right(w, s)
	next := x
	if g.length(xs) < g.length(x) {
		next = xs
	}
	r := g.leq(next, ws)
	g.leqMemo[key] = r
	return r
}

func (g *lifting) lower(w elem) map[elem]struct{} {
	if l, ok := g.lowers[w]; ok {
		return l
	}
	var out map[elem]struct{}
	if g.length(w) == 0 {
		out = map[elem]struct{}{w: {}}
	} else {
		s := g.firstDescent(w)
		base := g.lower(g.right(w, s))
		out = make(map[elem]struct{}, 2*len(base))
		for x := range base {
			out[x] = struct{}{}
			out[g.right(x, s)] = struct{}{}
		}
	}
	g.lowers[w] = out
	return out
}

func sortLex(es []elem) {
	sort.Slice(es, func(i, j int) bool { return lessLex(es[i], es[j]) })
}

func computeShape() shapeCert {
	g := newLifting()
	x, w := shapeX, shapeW
	var interval []elem
	for z := range g.lower(w) {
		if g.leq(x, z) {
			interval = append(interval, z)
		}
	}
	sortLex(interval)
	rank := g.length(w) - g.length(x)
	var atoms []elem
	for _, z := range interval {
		if g.length(z) == g.length(x)+1 {
			atoms = append(atoms, z)
		}
	}
	if !(rank == 11 && len(atoms) == 12) {
		fail("expected rank 11 and 12 atoms, got rank %d and %d atoms", rank, len(atoms))
	}
	// An identity-bottom interval of rank d has <=d atoms: the atoms are
	// distinct simple generators occurring in a reduced word of length d.
	if !(len(atoms) > rank) {
		fail("atom count %d does not exceed rank %d", len(atoms), rank)
	}

	var witness *witnessPayload
	for i := 0; i < len(atoms) && witness == nil; i++ {
		for j := i + 1; j < len(atoms); j++ {
			a, b := atoms[i], atoms[j]
			var common []elem
			for _, z := range interval {
				if g.leq(a, z) && g.leq(b, z) {
					common = append(common, z)
				}
			}
			sort.Slice(common, func(p, q int) bool {
				lp, lq := g.length(common[p]), g.length(common[q])
				if lp != lq {
					return lp < lq
				}
				return lessLex(common[p], common[q])
			})
			var minimal []elem
			for _, z := range common {
				dominated := false
				for _, t := range minimal {
					if g.leq(t, z) {
						dominated = true
						break
					}
				}
				if !dominated {
					minimal = append(minimal, z)
				}
			}
			if len(minimal) > 1 {
				// Independently verify every minimality claim over all common bounds.
				for _, z := range minimal {
					for _, t := range common {
						if t != z && g.leq(t, z) {
							fail("minimality claim fails: %v < %v", t, z)
						}
					}
				}
				for p := 0; p < len(minimal); p++ {
					for q := p + 1; q < len(minimal); q++ {
						if g.leq(minimal[p], minimal[q]) || g.leq(minimal[q], minimal[p]) {
							fail("minimal bounds %v and %v are comparable", minimal[p], minimal[q])
						}
					}
				}
				ranks := make([]int, len(minimal))
				for k, z := range minimal {
					ranks[k] = g.length(z) - g.length(x)
				}
				witness = &witnessPayload{
					Atoms:                       toLists([]elem{a, b}),
					MinimalCommonUpperBounds:    toLists(minimal),
					RelativeRanks:               ranks,
					AllCommonUpperBoundsChecked: len(common),
				}
				break
			}
		}
	}
	if witness == nil {
		fail("no pair of atoms has two minimal common upper bounds")
	}
	return shapeCert{
		Bottom:                 toLists([]elem{x})[0],
		Top:                    toLists([]elem{w})[0],
		IntervalSize:           len(interval),
		Rank:                   rank,
		Atoms:                  toLists(atoms),
		AtomCount:              len(atoms),
		PrincipalLowerObstruct: "atom_count > interval_rank",
		NonLatticeWitness:      *witness,
		Checks:                 "All interval vertices generated by subwords; membership and comparisons use exact lifting property.",
	}
}
