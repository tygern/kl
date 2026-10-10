package main

// A faithful emulation of the slot layout of CPython's set/frozenset hash
// table (Objects/setobject.c: linear probing of 9 slots, perturbation
// i = 5*i + 1 + perturb, resizing rules, and the exact construction paths of
// frozenset([x]) and frozenset | set). The Python original iterates frozensets
// of group elements when first computing Kazhdan-Lusztig values. That
// iteration order is the record order of the committed
// certificate affine_d4_kl_certificate.json; reproducing it byte for byte
// requires reproducing the iteration order of those sets.

const (
	linearProbes = 9
	perturbShift = 5
	minSetSize   = 8
)

// pySet stores element identifiers (interned, so identity is equality) with
// their CPython hash values. keys[i] < 0 marks an unused slot (the original
// never deletes, so there are no dummy entries).
type pySet struct {
	mask   int
	fill   int
	used   int
	keys   []int32
	hashes []uint64
}

func newPySet() *pySet {
	s := &pySet{mask: minSetSize - 1}
	s.keys = make([]int32, minSetSize)
	s.hashes = make([]uint64, minSetSize)
	for i := range s.keys {
		s.keys[i] = -1
	}
	return s
}

func (s *pySet) insertClean(key int32, hash uint64) {
	mask := uint64(s.mask)
	i := hash & mask
	perturb := hash
	for {
		probes := uint64(0)
		if i+linearProbes <= mask {
			probes = linearProbes
		}
		// CPython advances an entry pointer through up to probes+1 slots but
		// leaves i at the start of the group.
		for j := i; j <= i+probes; j++ {
			if s.keys[j] < 0 {
				s.keys[j] = key
				s.hashes[j] = hash
				return
			}
		}
		perturb >>= perturbShift
		i = (i*5 + 1 + perturb) & mask
	}
}

func (s *pySet) resize(minused int) {
	newsize := minSetSize
	for newsize <= minused {
		newsize <<= 1
	}
	oldKeys, oldHashes := s.keys, s.hashes
	s.mask = newsize - 1
	s.keys = make([]int32, newsize)
	s.hashes = make([]uint64, newsize)
	for i := range s.keys {
		s.keys[i] = -1
	}
	s.fill = s.used
	for i, k := range oldKeys {
		if k >= 0 {
			s.insertClean(k, oldHashes[i])
		}
	}
}

func (s *pySet) add(key int32, hash uint64) {
	mask := uint64(s.mask)
	i := hash & mask
	perturb := hash
	for {
		probes := uint64(0)
		if i+linearProbes <= mask {
			probes = linearProbes
		}
		for j := i; j <= i+probes; j++ {
			if s.keys[j] < 0 {
				s.fill++
				s.used++
				s.keys[j] = key
				s.hashes[j] = hash
				if uint64(s.fill)*5 < mask*3 {
					return
				}
				if s.used > 50000 {
					s.resize(s.used * 2)
				} else {
					s.resize(s.used * 4)
				}
				return
			}
			if s.hashes[j] == hash && s.keys[j] == key {
				return
			}
		}
		perturb >>= perturbShift
		i = (i*5 + 1 + perturb) & mask
	}
}

// merge is set_merge(so, other).
func (s *pySet) merge(o *pySet) {
	if o == s || o.used == 0 {
		return
	}
	if (s.fill+o.used)*5 >= s.mask*3 {
		s.resize((s.used + o.used) * 2)
	}
	if s.fill == 0 && s.mask == o.mask && o.fill == o.used {
		for i := 0; i <= o.mask; i++ {
			if o.keys[i] >= 0 {
				s.keys[i] = o.keys[i]
				s.hashes[i] = o.hashes[i]
			}
		}
		s.fill = o.fill
		s.used = o.used
		return
	}
	if s.fill == 0 {
		s.fill = o.used
		s.used = o.used
		for i := 0; i <= o.mask; i++ {
			if o.keys[i] >= 0 {
				s.insertClean(o.keys[i], o.hashes[i])
			}
		}
		return
	}
	for i := 0; i <= o.mask; i++ {
		if o.keys[i] >= 0 {
			s.add(o.keys[i], o.hashes[i])
		}
	}
}

// order lists the members in table (iteration) order.
func (s *pySet) order() []int32 {
	out := make([]int32, 0, s.used)
	for _, k := range s.keys {
		if k >= 0 {
			out = append(out, k)
		}
	}
	return out
}
