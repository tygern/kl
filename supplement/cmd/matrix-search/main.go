// Command matrix-search is the independent integer-matrix terminal search for
// the simply laced finite Weyl groups E6 and E7.
//
// It ports two original programs:
//
//	research/verify_e6_independent.py   (-rank 6)
//	research/verify_e7_matrices.cpp     (-rank 7)
//
// Group elements are integral matrices (stored by columns) acted on by the
// simple reflections in the faithful geometric representation; no root-index
// engine is used. This command imports only the Go standard library: it
// imports no internal package and none of the other engines, preserving the
// claim that the integer-matrix engine and the root-index engine
// (cmd/enumerate-bad) share no code.
//
// Usage (run in the work directory):
//
//	matrix-search -rank 6   writes results/e6-independent-certificate.json and prints it
//	matrix-search -rank 7   prints the E7 certificate JSON on stdout
//
// Any failed assertion of the original terminates the program with a
// non-zero exit status and a message on standard error.
package main

import (
	"flag"
	"fmt"
	"os"
)

// check is the Go counterpart of the originals' assert statements.
func check(cond bool, format string, args ...any) {
	if !cond {
		fmt.Fprintf(os.Stderr, "matrix-search: assertion failed: "+format+"\n", args...)
		os.Exit(1)
	}
}

func main() {
	rank := flag.Int("rank", 0, "rank of the Weyl group: 6 (E6) or 7 (E7)")
	flag.Parse()
	if flag.NArg() != 0 {
		fmt.Fprintln(os.Stderr, "matrix-search: positional arguments are not used")
		os.Exit(2)
	}
	switch *rank {
	case 6:
		runE6()
	case 7:
		runE7()
	default:
		fmt.Fprintln(os.Stderr, "matrix-search: -rank must be 6 or 7")
		os.Exit(2)
	}
}
