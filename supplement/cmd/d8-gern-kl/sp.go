package main

// Model sp: signed permutations in Gern's conventions.

type perm [8]int8

type spGroup struct{ n int }

func abs8(v int8) int {
	if v < 0 {
		return -int(v)
	}
	return int(v)
}

func (g spGroup) key(w perm) uint32 {
	var k uint32
	for i := 0; i < g.n; i++ {
		v := uint32(abs8(w[i]) - 1)
		if w[i] < 0 {
			v |= 8
		}
		k |= v << (4 * uint(i))
	}
	return k
}

func (g spGroup) length(w perm) int {
	l := 0
	for i := 0; i < g.n; i++ {
		for j := i + 1; j < g.n; j++ {
			if w[i] > w[j] {
				l++
			}
			if int(w[i])+int(w[j]) < 0 {
				l++
			}
		}
	}
	return l
}

func (g spGroup) rmul(w perm, s int) perm {
	if s == 0 {
		a, b := w[0], w[1]
		w[0] = -b
		w[1] = -a
	} else {
		w[s-1], w[s] = w[s], w[s-1]
	}
	return w
}

func (g spGroup) lmul(w perm, s int) perm {
	for i := 0; i < g.n; i++ {
		v := int(w[i])
		a := v
		sg := 1
		if v < 0 {
			a = -v
			sg = -1
		}
		if s == 0 {
			if a == 1 {
				w[i] = int8(-2 * sg)
			} else if a == 2 {
				w[i] = int8(-1 * sg)
			}
		} else {
			if a == s {
				w[i] = int8((s + 1) * sg)
			} else if a == s+1 {
				w[i] = int8(s * sg)
			}
		}
	}
	return w
}

func (g spGroup) refl(w perm, i, j int, neg bool) perm {
	w[i], w[j] = w[j], w[i]
	if neg {
		w[i] = -w[i]
		w[j] = -w[j]
	}
	return w
}

func (g spGroup) identity() perm {
	var e perm
	for i := 0; i < g.n; i++ {
		e[i] = int8(i + 1)
	}
	return e
}

// equal compares only the first n positions (the rest is unused).
func (g spGroup) equal(a, b perm) bool {
	for i := 0; i < g.n; i++ {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func spBuild(rank int, top perm, xword []int) *Model {
	g := spGroup{rank}
	n := rank
	var elems []perm
	idx := map[uint32]int32{}
	add := func(w perm) int32 {
		k := g.key(w)
		if id, ok := idx[k]; ok {
			return id
		}
		id := int32(len(elems))
		elems = append(elems, w)
		idx[k] = id
		return id
	}
	lookup := func(w perm) (int32, bool) {
		id, ok := idx[g.key(w)]
		return id, ok
	}
	must := func(w perm) int32 {
		id, ok := lookup(w)
		if !ok {
			fatal("internal error: element missing from ideal")
		}
		return id
	}
	add(top)
	for p := 0; p < len(elems); p++ {
		w := elems[p]
		lw := g.length(w)
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				for ng := 0; ng < 2; ng++ {
					u := g.refl(w, i, j, ng == 1)
					if g.length(u) < lw {
						add(u)
					}
				}
			}
		}
	}
	M := newModel(n, len(elems))
	for id := 0; id < M.N; id++ {
		w := elems[id]
		lw := g.length(w)
		M.length[id] = uint8(lw)
		for s := 0; s < 8; s++ {
			M.rmul[id][s] = -1
			M.lmul[id][s] = -1
		}
		for s := 0; s < n; s++ {
			u := g.rmul(w, s)
			if v, ok := lookup(u); ok {
				M.rmul[id][s] = v
			}
			if g.length(u) < lw {
				M.rdes[id] |= 1 << uint(s)
			}
			v := g.lmul(w, s)
			if t, ok := lookup(v); ok {
				M.lmul[id][s] = t
			}
			if g.length(v) < lw {
				M.ldes[id] |= 1 << uint(s)
			}
		}
		var rd uint8
		for s := 0; s < n; s++ {
			var d bool
			if s == 0 {
				d = int(w[0])+int(w[1]) < 0
			} else {
				d = w[s-1] > w[s]
			}
			if d {
				rd |= 1 << uint(s)
			}
		}
		if rd != M.rdes[id] { // Gern Prop. 2.2.4
			fatal("descent formula mismatch")
		}
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				for ng := 0; ng < 2; ng++ {
					u := g.refl(w, i, j, ng == 1)
					if g.length(u) == lw-1 {
						M.coatoms[id] = append(M.coatoms[id], must(u))
					}
				}
			}
		}
	}
	e := g.identity()
	M.idE = int(must(e))
	M.idTop = int(must(top))
	x := e
	for _, s := range xword {
		x = g.rmul(x, s-1)
	}
	if id, ok := lookup(x); ok {
		M.idX = int(id)
	} else {
		M.idX = -1
	}
	return M
}
