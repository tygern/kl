package main

// E7 integer-matrix terminal certificate: port of research/verify_e7_matrices.cpp.

import (
	"bufio"
	"math/bits"
	"os"
	"strconv"
)

const e7N = 7

type e7Mat [e7N * e7N]int8

var e7Edges = [6][2]int{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {2, 6}}

func e7Multiply(a *e7Mat, s int, adj *[e7N][]int) e7Mat {
	b := *a
	for k := 0; k < e7N; k++ {
		b[s*e7N+k] = -a[s*e7N+k]
	}
	for _, t := range adj[s] {
		for k := 0; k < e7N; k++ {
			v := int(a[t*e7N+k]) + int(a[s*e7N+k])
			check(v >= -128 && v <= 127, "int8 matrix entry overflow")
			b[t*e7N+k] = int8(v)
		}
	}
	return b
}

func e7Descents(a *e7Mat) uint {
	var mask uint
	for s := 0; s < e7N; s++ {
		pos, neg := false, false
		for k := 0; k < e7N; k++ {
			v := a[s*e7N+k]
			pos = pos || v > 0
			neg = neg || v < 0
		}
		check(pos != neg, "column %d is not of definite sign", s)
		if neg {
			mask |= 1 << uint(s)
		}
	}
	return mask
}

// e7Hash is FNV-1a over the entries, as in the original.
func e7Hash(a *e7Mat) uint64 {
	h := uint64(1469598103934665603)
	for _, c := range a {
		h = (h ^ uint64(uint8(c))) * 1099511628211
	}
	return h
}

func runE7() {
	const expected = 2903040
	var adj [e7N][]int
	for _, e := range e7Edges {
		adj[e[0]] = append(adj[e[0]], e[1])
		adj[e[1]] = append(adj[e[1]], e[0])
	}
	var identity e7Mat
	for s := 0; s < e7N; s++ {
		identity[s*e7N+s] = 1
	}
	el := make([]e7Mat, 0, expected)
	el = append(el, identity)
	// Open-addressing table of element ids (+1); 0 means empty.
	const tableSize = 1 << 23
	table := make([]uint32, tableSize)
	insert := func(a *e7Mat, id uint32) {
		i := e7Hash(a) & (tableSize - 1)
		for table[i] != 0 {
			i = (i + 1) & (tableSize - 1)
		}
		table[i] = id + 1
	}
	find := func(a *e7Mat) (uint32, bool) {
		i := e7Hash(a) & (tableSize - 1)
		for table[i] != 0 {
			if el[table[i]-1] == *a {
				return table[i] - 1, true
			}
			i = (i + 1) & (tableSize - 1)
		}
		return 0, false
	}
	insert(&identity, 0)
	action := make([][e7N]uint32, 0, expected)
	parent := make([]uint32, 1, expected)
	last := make([]uint8, 1, expected)
	length := make([]uint8, 1, expected)
	rd := make([]uint8, 0, expected)
	fc := make([]bool, 0, expected)
	for w := uint32(0); int(w) < len(el); w++ {
		var row [e7N]uint32
		for s := 0; s < e7N; s++ {
			b := e7Multiply(&el[w], s, &adj)
			if id, ok := find(&b); ok {
				row[s] = id
			} else {
				v := uint32(len(el))
				el = append(el, b)
				insert(&b, v)
				parent = append(parent, w)
				last = append(last, uint8(s))
				check(int(length[w])+1 < 256, "length overflow")
				length = append(length, length[w]+1)
				row[s] = v
			}
		}
		action = append(action, row)
		d := e7Descents(&el[w])
		rd = append(rd, uint8(d))
		good := true
		for s := 0; s < e7N; s++ {
			if d>>uint(s)&1 == 1 {
				check(int(length[row[s]])+1 == int(length[w]), "descent does not lower length")
				good = good && fc[row[s]]
			}
		}
		for _, e := range e7Edges {
			if d>>uint(e[0])&1 == 1 && d>>uint(e[1])&1 == 1 {
				good = false
			}
		}
		fc = append(fc, good)
	}
	check(len(el) == expected, "E7 element count %d != %d", len(el), expected)
	el = nil
	table = nil

	inverse := make([]uint32, expected)
	for w := uint32(0); w < expected; w++ {
		z := uint32(0)
		for v := w; v != 0; v = parent[v] {
			z = action[z][last[v]]
		}
		inverse[w] = z
	}
	left := func(w uint32, s int) uint32 { return inverse[action[inverse[w]][s]] }
	weak := func(w uint32) bool {
		for s := 0; s < e7N; s++ {
			if rd[w]>>uint(s)&1 == 1 {
				for _, t := range adj[s] {
					if rd[action[w][s]]>>uint(t)&1 == 1 {
						return false
					}
				}
			}
		}
		return true
	}
	leq := func(x, w uint32) bool {
		for x != w {
			if length[x] >= length[w] {
				return false
			}
			s := bits.TrailingZeros(uint(rd[w]))
			if rd[x]>>uint(s)&1 == 1 {
				x = action[x][s]
			}
			w = action[w][s]
		}
		return true
	}
	wordOf := func(w uint32) []int {
		var word []int
		for ; w != 0; w = parent[w] {
			word = append(word, int(last[w]))
		}
		for i, j := 0, len(word)-1; i < j; i, j = i+1, j-1 {
			word[i], word[j] = word[j], word[i]
		}
		return word
	}

	var fcs, bads []uint32
	var starChecks uint64
	var levels [64]uint64
	commuting := 0
	for w := uint32(0); w < expected; w++ {
		levels[length[w]]++
		if fc[w] {
			fcs = append(fcs, w)
			for side := 0; side < 2; side++ {
				for _, e := range e7Edges {
					s, t := e[0], e[1]
					pair := uint(1)<<uint(s) | uint(1)<<uint(t)
					ds := func(v uint32) uint {
						if side == 1 {
							return uint(rd[inverse[v]])
						}
						return uint(rd[v])
					}
					if bits.OnesCount(ds(w)&pair) != 1 {
						continue
					}
					valid := 0
					for _, u := range []int{s, t} {
						var v uint32
						if side == 1 {
							v = left(w, u)
						} else {
							v = action[w][u]
						}
						if bits.OnesCount(ds(v)&pair) == 1 {
							valid++
							check(fc[v], "star move leaves FC set")
						}
					}
					check(valid == 1, "star at element %d has %d moves", w, valid)
					starChecks++
				}
			}
		}
		if weak(w) && weak(inverse[w]) {
			if !fc[w] {
				bads = append(bads, w)
			} else {
				var support uint
				for v := w; v != 0; v = parent[v] {
					support |= 1 << uint(last[v])
				}
				check(bits.OnesCount(support) == int(length[w]), "commuting terminal repeats a letter")
				for _, e := range e7Edges {
					check(!(support>>uint(e[0])&1 == 1 && support>>uint(e[1])&1 == 1), "commuting terminal contains an edge")
				}
				commuting++
			}
		}
	}
	check(len(fcs) == 2670, "E7 fc count %d != 2670", len(fcs))
	check(len(bads) == 4, "E7 bad count %d != 4", len(bads))

	out := bufio.NewWriter(os.Stdout)
	defer out.Flush()
	printWord := func(w uint32) {
		out.WriteByte('[')
		for k, s := range wordOf(w) {
			if k > 0 {
				out.WriteByte(',')
			}
			out.WriteString(strconv.Itoa(s))
		}
		out.WriteByte(']')
	}
	out.WriteString(`{"representation":"independent integer matrices","order":` + strconv.Itoa(expected) +
		`,"fc_count":` + strconv.Itoa(len(fcs)) + `,"fc_star_checks":` + strconv.FormatUint(starChecks, 10) +
		`,"commuting_terminals":` + strconv.Itoa(commuting) + `,"length_distribution":[`)
	for k := 0; k < 64; k++ {
		if k > 0 {
			out.WriteByte(',')
		}
		out.WriteString(strconv.FormatUint(levels[k], 10))
	}
	out.WriteString(`],"bad":[`)
	for k, w := range bads {
		if k > 0 {
			out.WriteByte(',')
		}
		out.WriteString(`{"word":`)
		printWord(w)
		out.WriteString(`,"length":` + strconv.Itoa(int(length[w])) + `,"Rmask":` + strconv.Itoa(int(rd[w])) +
			`,"Lmask":` + strconv.Itoa(int(rd[inverse[w]])) + `,"bottoms":[`)
		first := true
		for _, x := range fcs {
			if rd[x]&rd[w] == rd[w] && rd[inverse[x]]&rd[inverse[w]] == rd[inverse[w]] && leq(x, w) {
				if !first {
					out.WriteByte(',')
				}
				first = false
				out.WriteString(`{"word":`)
				printWord(x)
				out.WriteString(`,"length":` + strconv.Itoa(int(length[x])) + `,"rank":` + strconv.Itoa(int(length[w])-int(length[x])) + `}`)
			}
		}
		out.WriteString(`]}`)
	}
	out.WriteString("]}\n")
}
