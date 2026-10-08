package main

import (
	"fmt"
	"math/bits"
	"sort"
	"sync"
	"sync/atomic"
)

const maxDeg = 14 // degrees stay below (l(w)-1)/2 = 12 for D8

// guardBound returns the overflow guard for the stored coefficients and mu
// values of a model with n elements: the largest power of two B with
// n*B*B + 2*B < 2^63, that is B = 2^floor((62 - bitlen(n))/2). Every stored
// polynomial and every mu is checked against B before it is used, and a
// recurrence step (compute) adds two stored polynomials and at most n
// products mu*coefficient of stored values, so no intermediate sum can
// overflow int64 before its result is checked by intern. For D8
// (n = 5,160,960) B = 2^19 = 524,288; the largest coefficient that occurs
// in the certificate is 233.
func guardBound(n int) int64 {
	return int64(1) << uint((62-bits.Len(uint(n)))/2)
}

type poly [maxDeg]int64

type entry struct {
	x int32
	p int32
}

type pending struct {
	x int32
	p poly
}

type muPair struct {
	z  int32
	mu int64
}

// KL computes Kazhdan-Lusztig polynomials on extremal pairs by layers of
// length.  Within a layer the workers only read tables of lower layers; their
// results go to per-element slots, and the polynomials are interned
// sequentially in layer order after the join, so neither scheduling nor the
// thread count influences any stored value (this replaces the C++ append-only
// chunked store that avoided a data race).
type KL struct {
	M       *Model
	table   [][]entry
	mulist  [][]muPair
	polys   []poly
	polyidx map[poly]int32
	guard   int64 // overflow guard for coefficients and mu values, see guardBound
}

func newKL(m *Model) *KL {
	return &KL{M: m, table: make([][]entry, m.N), mulist: make([][]muPair, m.N), polyidx: map[poly]int32{}, guard: guardBound(m.N)}
}

func (k *KL) intern(p poly) int32 {
	if id, ok := k.polyidx[p]; ok {
		return id
	}
	for _, c := range p {
		if c > k.guard || c < -k.guard {
			fatal("coefficient overflow guard")
		}
	}
	id := int32(len(k.polys))
	k.polys = append(k.polys, p)
	k.polyidx[p] = id
	return id
}

func deg(p *poly) int {
	for i := maxDeg - 1; i >= 0; i-- {
		if p[i] != 0 {
			return i
		}
	}
	return -1
}

func (k *KL) extremal(x, y int) bool {
	M := k.M
	return M.rdes[y]&^M.rdes[x] == 0 && M.ldes[y]&^M.ldes[x] == 0
}

func (k *KL) extremalize(u, z int) int {
	M := k.M
	for {
		if r := M.rdes[z] &^ M.rdes[u]; r != 0 {
			u = int(M.rmul[u][bits.TrailingZeros8(r)])
			if u < 0 {
				return -1
			}
			continue
		}
		if l := M.ldes[z] &^ M.ldes[u]; l != 0 {
			u = int(M.lmul[u][bits.TrailingZeros8(l)])
			if u < 0 {
				return -1
			}
			continue
		}
		return u
	}
}

func (k *KL) lookupExt(xe, z int) *poly {
	t := k.table[z]
	i := sort.Search(len(t), func(i int) bool { return int(t[i].x) >= xe })
	if i == len(t) || int(t[i].x) != xe {
		return nil
	}
	return &k.polys[t[i].p]
}

func (k *KL) lookup(u, z int) *poly {
	if u < 0 {
		return nil
	}
	xe := k.extremalize(u, z)
	if xe < 0 {
		return nil
	}
	return k.lookupExt(xe, z)
}

func (k *KL) compute(x, y, s, v int) poly {
	M := k.M
	var P poly
	c := int((M.rdes[x] >> uint(s)) & 1)
	xs := int(M.rmul[x][s])
	if a := k.lookup(xs, v); a != nil {
		for i := 0; i+(1-c) < maxDeg; i++ {
			P[i+1-c] += a[i]
		}
	}
	if b := k.lookup(x, v); b != nil {
		for i := 0; i+c < maxDeg; i++ {
			P[i+c] += b[i]
		}
	}
	for _, zm := range k.mulist[v] {
		z := int(zm.z)
		if (M.rdes[z]>>uint(s))&1 == 0 {
			continue
		}
		pz := k.lookup(x, z)
		if pz == nil {
			continue
		}
		sh := (int(M.length[y]) - int(M.length[z])) / 2
		mu := zm.mu
		for i := 0; i+sh < maxDeg; i++ {
			P[i+sh] -= mu * pz[i]
		}
	}
	return P
}

func (k *KL) finish(y int, res []entry, npairs *int64) {
	M := k.M
	k.table[y] = res
	*npairs += int64(len(res))
	var mu []muPair
	for _, z := range M.coatoms[y] {
		mu = append(mu, muPair{z, 1})
	}
	for _, e := range res {
		cod := int(M.length[y]) - int(M.length[e.x])
		if cod <= 1 || cod%2 == 0 {
			continue
		}
		p := &k.polys[e.p]
		d := (cod - 1) / 2
		if p[d] != 0 {
			if p[d] > k.guard || p[d] < -k.guard {
				fatal("mu overflow guard")
			}
			mu = append(mu, muPair{e.x, p[d]})
		}
	}
	sort.Slice(mu, func(a, b int) bool {
		if mu[a].z != mu[b].z {
			return mu[a].z < mu[b].z
		}
		return mu[a].mu < mu[b].mu
	})
	out := mu[:0]
	for i, m := range mu {
		if i == 0 || m != mu[i-1] {
			out = append(out, m)
		}
	}
	mu = out
	for i := 1; i < len(mu); i++ {
		if mu[i].z == mu[i-1].z {
			fatal("conflicting mu values")
		}
	}
	k.mulist[y] = mu
}

func (k *KL) run(nthreads int) {
	M := k.M
	order := make([]int, M.N)
	for i := range order {
		order[i] = i
	}
	sort.SliceStable(order, func(a, b int) bool { return M.length[order[a]] < M.length[order[b]] })
	maxlen := int(M.length[order[len(order)-1]])
	var npairs int64
	pos := 0
	for L := 0; L <= maxlen; L++ {
		start := pos
		for pos < len(order) && int(M.length[order[pos]]) == L {
			pos++
		}
		layer := order[start:pos]
		results := make([][]pending, len(layer))
		var next int64 = -1
		worker := func() {
			var buf []int32
			stamps := make([]uint32, M.N)
			var cur uint32
			for {
				i := int(atomic.AddInt64(&next, 1))
				if i >= len(layer) {
					break
				}
				y := layer[i]
				var res []pending
				if M.length[y] == 0 {
					var one poly
					one[0] = 1
					results[i] = []pending{{int32(y), one}}
					continue
				}
				s := bits.TrailingZeros8(M.rdes[y])
				v := int(M.rmul[y][s])
				cur++
				buf = buf[:0]
				for _, e := range k.table[v] {
					if stamps[e.x] != cur {
						stamps[e.x] = cur
						buf = append(buf, e.x)
					}
				}
				for p := 0; p < len(buf); p++ {
					x := int(buf[p])
					r := M.rdes[v] & M.rdes[x]
					for r != 0 {
						g := bits.TrailingZeros8(r)
						r &= r - 1
						u := M.rmul[x][g]
						if u >= 0 && stamps[u] != cur {
							stamps[u] = cur
							buf = append(buf, u)
						}
					}
					l := M.ldes[v] & M.ldes[x]
					for l != 0 {
						g := bits.TrailingZeros8(l)
						l &= l - 1
						u := M.lmul[x][g]
						if u >= 0 && stamps[u] != cur {
							stamps[u] = cur
							buf = append(buf, u)
						}
					}
				}
				nv := len(buf)
				for p := 0; p < nv; p++ {
					x := int(buf[p])
					xs := M.rmul[x][s]
					if xs >= 0 && stamps[xs] != cur {
						stamps[xs] = cur
						buf = append(buf, xs)
					}
				}
				for _, x32 := range buf {
					x := int(x32)
					if k.extremal(x, y) {
						P := k.compute(x, y, s, v)
						cod := int(M.length[y]) - int(M.length[x])
						if d := deg(&P); d > (cod-1)/2 || d >= maxDeg-1 {
							fatal("degree bound violated")
						}
						res = append(res, pending{x32, P})
					}
				}
				sort.Slice(res, func(a, b int) bool { return res[a].x < res[b].x })
				results[i] = res
			}
		}
		var wg sync.WaitGroup
		for t := 0; t < nthreads; t++ {
			wg.Add(1)
			go func() { defer wg.Done(); worker() }()
		}
		wg.Wait()
		for i, y := range layer {
			res := make([]entry, len(results[i]))
			for j, pe := range results[i] {
				res[j] = entry{pe.x, k.intern(pe.p)}
			}
			results[i] = nil
			k.finish(y, res, &npairs)
		}
	}
}

func polyJSON(p *poly) string {
	d := deg(p)
	if d < 0 {
		d = 0
	}
	s := "["
	for i := 0; i <= d; i++ {
		if i > 0 {
			s += ","
		}
		s += fmt.Sprint(p[i])
	}
	return s + "]"
}
