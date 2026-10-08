// Command fc-maxima enumerates all fully commutative (FC) elements of the
// generalized Coxeter group E_n (chain 0-1-...-(n-2), node n-1 attached to
// node 2), level by level, and records the count, the length distribution, the
// maximum length and the reduced words of the longest elements.  Elements are
// stored as height vectors h(w)_i = height of w(alpha_i); right multiplication
// by s_j acts as h_i += h_j for i adjacent to j, then h_j = -h_j, and s_j is a
// right descent iff h_j < 0.  Full commutativity is decided by the recurrence
// "w is FC iff R(w) is pairwise commuting and w t is FC for every t in R(w)".
//
// Ports research/review_checks/fc_maxima.cpp.  Imports only the Go standard
// library and no other engine of this module (self-contained by design).
//
// Usage: fc-maxima -out results/fc-maxima-certificate.json -ranks 6,7,8,9,10,11,12,13
package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
)

const maxN = 16

// Key is the height vector of an element (unused entries stay zero).
type Key [maxN]int16

type engine struct {
	n   int
	T   int
	adj [][]int
	A   [20][20]bool
}

func shardOf(k *Key) uint64 {
	x := uint64(1469598103934665603)
	for i := 0; i < maxN; i++ {
		x ^= uint64(uint16(k[i]))
		x *= 1099511628211
		x ^= x >> 29
	}
	return (x * 0x9E3779B97F4A7C15) >> 40
}

func (e *engine) rmul(w *Key, j int) {
	for _, i := range e.adj[j] {
		v := int(w[i]) + int(w[j])
		if v > 32000 || v < -32000 {
			fmt.Fprintln(os.Stderr, "height overflow")
			os.Exit(3)
		}
		w[i] = int16(v)
	}
	w[j] = -w[j]
}

type outcome struct {
	count   int64
	maxlen  int
	levels  []int64
	longest []string
}

func (e *engine) enumerate(rank int) outcome {
	n := rank
	e.n = n
	e.adj = make([][]int, n)
	for i := 0; i < n-2; i++ {
		e.adj[i] = append(e.adj[i], i+1)
		e.adj[i+1] = append(e.adj[i+1], i)
	}
	e.adj[n-1] = append(e.adj[n-1], 2)
	e.adj[2] = append(e.adj[2], n-1)
	e.A = [20][20]bool{}
	for i := 0; i < n; i++ {
		for _, j := range e.adj[i] {
			e.A[i][j] = true
		}
	}
	T := e.T
	var start Key
	for i := 0; i < n; i++ {
		start[i] = 1
	}
	prev := make([]map[Key]struct{}, T)
	nxt := make([]map[Key]struct{}, T)
	for t := 0; t < T; t++ {
		prev[t] = map[Key]struct{}{}
	}
	prev[shardOf(&start)%uint64(T)][start] = struct{}{}
	out := outcome{count: 1, maxlen: 0, levels: []int64{1}}
	for {
		buckets := make([][][]Key, T)
		for t := range buckets {
			buckets[t] = make([][]Key, T)
		}
		var wg sync.WaitGroup
		for t := 0; t < T; t++ {
			wg.Add(1)
			go func(t int) {
				defer wg.Done()
				for u := range prev[t] {
					for s := 0; s < n; s++ {
						if u[s] < 0 {
							continue // s in R(u)
						}
						v := u
						e.rmul(&v, s)
						var R [20]int
						nr := 0
						for j := 0; j < n; j++ {
							if v[j] < 0 {
								R[nr] = j
								nr++
							}
						}
						ok := true
						for a := 0; a < nr && ok; a++ {
							for b := a + 1; b < nr; b++ {
								if e.A[R[a]][R[b]] {
									ok = false
									break
								}
							}
						}
						if !ok {
							continue
						}
						for a := 0; a < nr; a++ {
							if R[a] == s {
								continue // v s = u is FC by construction
							}
							x := v
							e.rmul(&x, R[a])
							if _, in := prev[shardOf(&x)%uint64(T)][x]; !in {
								ok = false
								break
							}
						}
						if !ok {
							continue
						}
						buckets[t][shardOf(&v)%uint64(T)] = append(buckets[t][shardOf(&v)%uint64(T)], v)
					}
				}
			}(t)
		}
		wg.Wait()
		for t := 0; t < T; t++ {
			wg.Add(1)
			go func(t int) {
				defer wg.Done()
				tot := 0
				for s := 0; s < T; s++ {
					tot += len(buckets[s][t])
				}
				m := make(map[Key]struct{}, tot)
				for s := 0; s < T; s++ {
					for _, k := range buckets[s][t] {
						m[k] = struct{}{}
					}
					buckets[s][t] = nil
				}
				nxt[t] = m
			}(t)
		}
		wg.Wait()
		var cnt int64
		for t := 0; t < T; t++ {
			cnt += int64(len(nxt[t]))
		}
		if cnt == 0 {
			break
		}
		out.maxlen++
		out.count += cnt
		out.levels = append(out.levels, cnt)
		prev, nxt = nxt, prev
	}
	// Reduced words of the longest elements (strip the smallest right descent repeatedly).
	for t := 0; t < T; t++ {
		for k := range prev[t] {
			w := k
			var word []int
			for {
				d := -1
				for j := 0; j < n; j++ {
					if w[j] < 0 {
						d = j
						break
					}
				}
				if d < 0 {
					break
				}
				e.rmul(&w, d)
				word = append(word, d)
			}
			parts := make([]string, len(word))
			for i := range word {
				parts[len(word)-1-i] = strconv.Itoa(word[i])
			}
			out.longest = append(out.longest, strings.Join(parts, " "))
		}
	}
	sort.Strings(out.longest)
	return out
}

func main() {
	outPath := flag.String("out", "", "output JSON path")
	ranksArg := flag.String("ranks", "", "comma-separated ranks (between 6 and 15)")
	threads := flag.Int("threads", 0, "number of shards/goroutines (0 = min(8, CPUs))")
	flag.Parse()
	if *outPath == "" || *ranksArg == "" {
		fmt.Fprintln(os.Stderr, "usage: fc-maxima -out FILE -ranks 6,7,...")
		os.Exit(2)
	}
	var ranks []int
	for _, s := range strings.Split(*ranksArg, ",") {
		r, err := strconv.Atoi(strings.TrimSpace(s))
		if err != nil {
			fmt.Fprintf(os.Stderr, "bad rank %q\n", s)
			os.Exit(2)
		}
		if r < 6 || r > 15 {
			fmt.Fprintln(os.Stderr, "rank out of range")
			os.Exit(2)
		}
		ranks = append(ranks, r)
	}
	T := *threads
	if T <= 0 {
		T = runtime.NumCPU()
		if T > 8 {
			T = 8
		}
		if T < 1 {
			T = 1
		}
	}
	type kv struct {
		count  int64
		maxlen int
	}
	known := map[int]kv{6: {662, 16}, 7: {2670, 27}, 8: {10846, 29}, 9: {44199, 44},
		10: {180438, 55}, 11: {737762, 66}, 12: {3021000, 78}, 13: {12387990, 92}}
	failures := 0
	f, err := os.Create(*outPath)
	if err != nil {
		fmt.Fprintln(os.Stderr, "cannot open output")
		os.Exit(2)
	}
	w := bufio.NewWriter(f)
	fmt.Fprint(w, "{\n  \"labelling\": \"E_n: chain 0-1-...-(n-2), node n-1 attached to node 2\",\n")
	fmt.Fprint(w, "  \"method\": \"level-wise closure under right ascents with the recurrence: w FC iff R(w) commuting and wt FC for every t in R(w); elements stored as height vectors h(w)_i = ht(w(alpha_i))\",\n")
	fmt.Fprint(w, "  \"ranks\": {\n")
	for a, rank := range ranks {
		eng := &engine{T: T}
		o := eng.enumerate(rank)
		exp, hasExp := known[rank]
		ok := true
		if hasExp {
			ok = o.count == exp.count && o.maxlen == exp.maxlen
			if !ok {
				fmt.Fprintf(os.Stderr, "EXPECTATION FAILED for E%d: %d / %d\n", rank, o.count, o.maxlen)
				failures++
			}
		}
		fmt.Fprintf(w, "    \"E%d\": {\"fully_commutative_count\": %d, \"maximum_length\": %d, \"longest_elements\": %d,\n", rank, o.count, o.maxlen, len(o.longest))
		fmt.Fprint(w, "      \"length_distribution\": [")
		for i, l := range o.levels {
			if i > 0 {
				fmt.Fprint(w, ",")
			}
			fmt.Fprintf(w, "%d", l)
		}
		fmt.Fprint(w, "],\n      \"longest_reduced_words\": [")
		for i, s := range o.longest {
			if i > 0 {
				fmt.Fprint(w, ", ")
			}
			fmt.Fprintf(w, "\"%s\"", s)
		}
		me := "null"
		if hasExp {
			me = "false"
			if ok {
				me = "true"
			}
		}
		comma := ""
		if a+1 < len(ranks) {
			comma = ","
		}
		fmt.Fprintf(w, "],\n      \"matches_expected\": %s}%s\n", me, comma)
		note := ""
		if hasExp {
			note = " (UNEXPECTED)"
			if ok {
				note = " (as expected)"
			}
		}
		fmt.Printf("E%d: %d fully commutative elements, maximum length %d%s\n", rank, o.count, o.maxlen, note)
	}
	status := "passed"
	if failures > 0 {
		status = "FAILED"
	}
	fmt.Fprintf(w, "  },\n  \"status\": \"%s\"\n}\n", status)
	if err := w.Flush(); err != nil {
		fmt.Fprintln(os.Stderr, "write failed:", err)
		os.Exit(2)
	}
	f.Close()
	if failures > 0 {
		os.Exit(1)
	}
}
