package main

// Finite D_n signed permutations in one-line notation, with the independent
// length formula, right action, descents and lifting-property Bruhat
// comparison of SparseCoxeter (kind 'D') in sparse_kl.py. Words multiply on
// the right; s_0 maps (a,b) in positions 0,1 to (-b,-a) and s_i, i>0, swaps
// positions i-1,i.

// signedPerm is a signed permutation of 1..n (n <= maxSigned).
type signedPerm []int

func (w signedPerm) clone() signedPerm { return append(signedPerm(nil), w...) }

func equalPerm(a, b signedPerm) bool {
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// lengthD is the type D length: inversions plus negative-sum pairs.
func lengthD(w signedPerm) int {
	n := len(w)
	total := 0
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if w[i] > w[j] {
				total++
			}
			if -w[i] > w[j] {
				total++
			}
		}
	}
	return total
}

func rightD(w signedPerm, s int) signedPerm {
	a := w.clone()
	if s > 0 {
		a[s-1], a[s] = a[s], a[s-1]
	} else {
		a[0], a[1] = -a[1], -a[0]
	}
	return a
}

// firstDescent returns the smallest right descent s of w (w must not be the
// identity; the original takes descents(w)[0]).
func firstDescent(w signedPerm) int {
	lw := lengthD(w)
	for s := 0; s < len(w); s++ {
		if lengthD(rightD(w, s)) < lw {
			return s
		}
	}
	fail("element has no right descent")
	return -1
}

// leqD reports x <= w in Bruhat order by the lifting property with the first
// right descent of w (iterative form of SparseCoxeter.leq).
func leqD(x, w signedPerm) bool {
	for {
		if equalPerm(x, w) {
			return true
		}
		lx, lw := lengthD(x), lengthD(w)
		if lx >= lw {
			return false
		}
		s := firstDescent(w)
		xs, ws := rightD(x, s), rightD(w, s)
		if lengthD(xs) < lx {
			x = xs
		}
		w = ws
	}
}

type coatomRow struct {
	N                int      `json:"n"`
	Length           int      `json:"length"`
	Rank             int      `json:"rank"`
	Coatoms          int      `json:"coatoms"`
	AllLowerCovers   int      `json:"all_lower_covers"`
	CoverReflections [][3]int `json:"cover_reflections"`
}

// gernPair returns the Gern bad pair (x, w) in D_n for even n.
func gernPair(n int) (x, w signedPerm) {
	w = make(signedPerm, n)
	for i := 1; i <= n; i++ {
		switch {
		case i == 1:
			if (n/2)%2 == 0 {
				w[i-1] = 1
			} else {
				w[i-1] = -1
			}
		case i%2 == 1:
			w[i-1] = i
		default:
			w[i-1] = -(n + 2 - i)
		}
	}
	x = signedPerm{-1, -2}
	for i := 3; i <= n; i++ {
		if i%2 == 1 {
			x = append(x, i+1)
		} else {
			x = append(x, i-1)
		}
	}
	return x, w
}

func gernRow(n int) coatomRow {
	x, w := gernPair(n)
	lw := lengthD(w)
	coatoms := [][3]int{}
	allCovers := 0
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			for _, sign := range []int{1, -1} {
				z := w.clone()
				z[i], z[j] = sign*w[j], sign*w[i]
				if lengthD(z) == lw-1 {
					allCovers++
					if leqD(x, z) {
						coatoms = append(coatoms, [3]int{i + 1, j + 1, sign})
					}
				}
			}
		}
	}
	return coatomRow{N: n, Length: lw, Rank: lw - lengthD(x), Coatoms: len(coatoms),
		AllLowerCovers: allCovers, CoverReflections: coatoms}
}

// runGernCoatoms computes the rows for n = from, from+2, ..., to.
func runGernCoatoms(from, to int) []coatomRow {
	rows := []coatomRow{}
	for n := from; n <= to; n += 2 {
		rows = append(rows, gernRow(n))
	}
	return rows
}
