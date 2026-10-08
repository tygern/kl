package main

import (
	"math/big"
)

// poly is an exact polynomial in Q[r,k,q]. Keys are exponent triples; zero
// coefficients are never stored, so map equality is polynomial equality.
type poly map[[3]int]*big.Rat

func normalize(d poly) poly {
	out := poly{}
	for a, b := range d {
		if b.Sign() != 0 {
			out[a] = new(big.Rat).Set(b)
		}
	}
	return out
}

// cst returns the constant polynomial n.
func cst(n int64) poly {
	return frac(n, 1)
}

// frac returns the constant polynomial num/den.
func frac(num, den int64) poly {
	v := big.NewRat(num, den)
	if v.Sign() == 0 {
		return poly{}
	}
	return poly{{0, 0, 0}: v}
}

// variable returns the i-th generator (0=r, 1=k, 2=q).
func variable(i int) poly {
	var e [3]int
	e[i] = 1
	return poly{e: big.NewRat(1, 1)}
}

func (p poly) add(x poly) poly {
	d := poly{}
	for a, b := range p {
		d[a] = new(big.Rat).Set(b)
	}
	for a, b := range x {
		if cur, ok := d[a]; ok {
			d[a] = new(big.Rat).Add(cur, b)
		} else {
			d[a] = new(big.Rat).Set(b)
		}
	}
	return normalize(d)
}

func (p poly) neg() poly {
	d := poly{}
	for a, b := range p {
		d[a] = new(big.Rat).Neg(b)
	}
	return normalize(d)
}

func (p poly) sub(x poly) poly { return p.add(x.neg()) }

func (p poly) mul(x poly) poly {
	d := poly{}
	for a, b := range p {
		for c, e := range x {
			f := [3]int{a[0] + c[0], a[1] + c[1], a[2] + c[2]}
			t := new(big.Rat).Mul(b, e)
			if cur, ok := d[f]; ok {
				d[f] = new(big.Rat).Add(cur, t)
			} else {
				d[f] = t
			}
		}
	}
	return normalize(d)
}

func (p poly) eq(x poly) bool {
	if len(p) != len(x) {
		return false
	}
	for a, b := range p {
		c, ok := x[a]
		if !ok || b.Cmp(c) != 0 {
			return false
		}
	}
	return true
}
