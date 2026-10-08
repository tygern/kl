package main

// Model geo: the geometric representation on the root lattice.

// mat[j*8+i] = coefficient of alpha_i in w(alpha_j).
type mat [64]int8
type vec [8]int

type geoGroup struct {
	n        int
	A        [8][8]int
	posroots []vec
}

func narrow(v int) int8 {
	if v < -128 || v > 127 {
		fatal("matrix entry out of int8 range")
	}
	return int8(v)
}

func (g *geoGroup) wapply(w *mat, v vec) vec {
	var r vec
	for j := 0; j < g.n; j++ {
		if v[j] != 0 {
			for i := 0; i < g.n; i++ {
				r[i] += v[j] * int(w[j*8+i])
			}
		}
	}
	return r
}

func (g *geoGroup) pairing(v, b vec) int {
	s := 0
	for i := 0; i < g.n; i++ {
		if v[i] != 0 {
			for j := 0; j < g.n; j++ {
				s += v[i] * g.A[i][j] * b[j]
			}
		}
	}
	return s
}

func (g *geoGroup) reflectv(v, b vec) vec {
	c := g.pairing(v, b)
	r := v
	for i := 0; i < g.n; i++ {
		r[i] -= c * b[i]
	}
	return r
}

func (g *geoGroup) isneg(v vec) bool {
	nz := false
	for i := 0; i < g.n; i++ {
		if v[i] > 0 {
			return false
		}
		if v[i] < 0 {
			nz = true
		}
	}
	return nz
}

func (g *geoGroup) length(w *mat) int {
	l := 0
	for _, b := range g.posroots {
		if g.isneg(g.wapply(w, b)) {
			l++
		}
	}
	return l
}

func (g *geoGroup) identity() mat {
	var m mat
	for i := 0; i < g.n; i++ {
		m[i*8+i] = 1
	}
	return m
}

func (g *geoGroup) rrefl(w *mat, b vec) mat {
	wb := g.wapply(w, b)
	r := *w
	for j := 0; j < g.n; j++ {
		var ej vec
		ej[j] = 1
		c := g.pairing(ej, b)
		if c != 0 {
			for i := 0; i < g.n; i++ {
				r[j*8+i] = narrow(int(r[j*8+i]) - c*wb[i])
			}
		}
	}
	return r
}

func (g *geoGroup) rmul(w *mat, s int) mat {
	var b vec
	b[s] = 1
	return g.rrefl(w, b)
}

func (g *geoGroup) lmul(w *mat, s int) mat {
	var b vec
	b[s] = 1
	r := *w
	for j := 0; j < g.n; j++ {
		var col vec
		for i := 0; i < g.n; i++ {
			col[i] = int(w[j*8+i])
		}
		c2 := g.reflectv(col, b)
		for i := 0; i < g.n; i++ {
			r[j*8+i] = narrow(c2[i])
		}
	}
	return r
}

func geoBuild(rank int, topword, xword []int) *Model {
	g := &geoGroup{n: rank}
	n := rank
	for i := 0; i < 8; i++ {
		for j := 0; j < 8; j++ {
			if i == j {
				g.A[i][j] = 2
			}
		}
	}
	edge := func(a, b int) { g.A[a-1][b-1] = -1; g.A[b-1][a-1] = -1 }
	edge(1, 3)
	edge(2, 3)
	for i := 3; i < n; i++ {
		edge(i, i+1)
	}
	{
		var roots []vec
		seen := map[vec]bool{}
		for i := 0; i < n; i++ {
			var v vec
			v[i] = 1
			roots = append(roots, v)
			seen[v] = true
		}
		for p := 0; p < len(roots); p++ {
			for s := 0; s < n; s++ {
				var b vec
				b[s] = 1
				r := g.reflectv(roots[p], b)
				pos, nz := true, false
				for i := 0; i < n; i++ {
					if r[i] < 0 {
						pos = false
					}
					if r[i] != 0 {
						nz = true
					}
				}
				if pos && nz && !seen[r] {
					seen[r] = true
					roots = append(roots, r)
				}
			}
		}
		g.posroots = roots
		if len(g.posroots) != n*(n-1) {
			fatal("wrong number of positive roots")
		}
	}
	top := g.identity()
	for _, s := range topword {
		top = g.rmul(&top, s-1)
	}
	var elems []mat
	idx := map[mat]int32{}
	add := func(w mat) int32 {
		if id, ok := idx[w]; ok {
			return id
		}
		id := int32(len(elems))
		elems = append(elems, w)
		idx[w] = id
		return id
	}
	must := func(w mat) int32 {
		id, ok := idx[w]
		if !ok {
			fatal("internal error: element missing from ideal")
		}
		return id
	}
	add(top)
	for p := 0; p < len(elems); p++ {
		w := elems[p]
		lw := g.length(&w)
		for _, b := range g.posroots {
			u := g.rrefl(&w, b)
			if g.length(&u) < lw {
				add(u)
			}
		}
	}
	M := newModel(n, len(elems))
	for id := 0; id < M.N; id++ {
		w := elems[id]
		lw := g.length(&w)
		M.length[id] = uint8(lw)
		for s := 0; s < 8; s++ {
			M.rmul[id][s] = -1
			M.lmul[id][s] = -1
		}
		for s := 0; s < n; s++ {
			u := g.rmul(&w, s)
			if t, ok := idx[u]; ok {
				M.rmul[id][s] = t
			}
			if g.length(&u) < lw {
				M.rdes[id] |= 1 << uint(s)
			}
			v := g.lmul(&w, s)
			if t, ok := idx[v]; ok {
				M.lmul[id][s] = t
			}
			if g.length(&v) < lw {
				M.ldes[id] |= 1 << uint(s)
			}
		}
		var rd uint8
		for s := 0; s < n; s++ {
			var col vec
			for i := 0; i < n; i++ {
				col[i] = int(w[s*8+i])
			}
			if g.isneg(col) {
				rd |= 1 << uint(s)
			}
		}
		if rd != M.rdes[id] {
			fatal("geo descent mismatch")
		}
		for _, b := range g.posroots {
			u := g.rrefl(&w, b)
			if g.length(&u) == lw-1 {
				M.coatoms[id] = append(M.coatoms[id], must(u))
			}
		}
	}
	e := g.identity()
	M.idE = int(must(e))
	M.idTop = int(must(top))
	x := e
	for _, s := range xword {
		x = g.rmul(&x, s-1)
	}
	if id, ok := idx[x]; ok {
		M.idX = int(id)
	} else {
		M.idX = -1
	}
	return M
}
