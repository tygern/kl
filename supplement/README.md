# Proof supplement: v0.2.0

Author: Tyson Gern (Initial Capacity, tyson@initialcapacity.io, ORCID 0009-0003-8288-8786).

This supplement reproduces the computational parts of the manuscript `results/exceptional-leading.tex` (Leading Kazhdan--Lusztig coefficients with fully commutative lower endpoints in type E). It is prepared for release v0.2.0 of [github.com/tygern/kl](https://github.com/tygern/kl); release v0.1.0 remains the historical first version. All code in this archive is licensed under the MIT License (file `LICENSE`, copyright 2026 Tyson Gern); the manuscript text is licensed under CC BY 4.0.

This supplement, like the manuscript it accompanies, was produced with substantial assistance from large language model systems: OpenAI's GPT-6.1 Sol and GPT-6 Astra wrote the enumeration and verification programs and the internal review notes, Anthropic's Claude Fable 5.1 carried out an independent adversarial review whose certificate programs are included under `research/review_checks/`, and the author directed the work, checked the results and takes full responsibility for the content. The file `research/en_uniform/referee.txt` is an AI-generated internal review note, not peer review.

## Run

Unzip the archive and enter its `exceptional-leading-proof` directory. Use Python 3.10 or later and an installed C++17 compiler:

```sh
python3 run_proofs.py
```

The runner searches for `clang++`, `c++`, or `g++`; a particular installed executable can be selected with `--compiler /path/to/compiler`. No network, installation, third-party Python packages, or `uv` are required. Do not use Python's `-O` option: assertions are part of verification. Compiler assertions are also retained. The default run takes under two minutes on a 2026 laptop (about 45 s of wall-clock time when the 10-core machine is otherwise idle, roughly twice that when it is busy) and needs about 500 MB of memory.

```sh
python3 run_proofs.py --full
```

The full mode additionally reruns the complete E6 and E7 root-index and integer-matrix searches (all 2,903,040 E7 elements), the second (geometric) model of the D6/D8 Kazhdan--Lusztig computation, and the larger uniform-family checks (ranks up to E49 and more Bruhat cover sets); it compares every regenerated snapshot with the shipped one, ignoring only wall-clock timing fields. It takes under three minutes on the same machine (about 100 s when idle, longer under load) and uses slightly more memory than the default. Both modes retain their isolated working trees, binaries, per-step stdout/stderr, stable summary, and status under `runs/`. The initial `payload/` stays unchanged; repeated runs are independent. An interrupted or failed run is retained for diagnosis. The runner first checks the SHA256 inventory, then requires the stable proof outcomes to equal `expected-summary.json`. Timing, compiler paths, and generated hashes of mutable output files are not treated as stable mathematical output.

## What the default recomputes

* E6/E7/E8 complete terminal classifications by flat and recursive parabolic pruning, and the separate E8/D7 classification. It checks completeness through the exhaustive source algorithms and compares all terminal matrices; a `complete: true` field alone is not accepted as an exhaustiveness argument.
* E6/E7/E8/E9 FC catalogues, including independent E8 closure and E9 closure verification. The E9 catalogue and its maximum length 44 are used in the affine cover proof. A second, independent FC enumeration by height vectors (`research/review_checks/fc_maxima.cpp`) recomputes the counts and maximum lengths of E6 to E13 (E10 to E13: 180,438/55, 737,762/66, 3,021,000/78, 12,387,990/92).
* Descent sets and support matchings for all eleven printed finite table rows, and equality of their matrices with the regenerated terminal list.
* The structure of the exceptional terminals (`terminal_structure.py`): the rows of Table 1 of lengths 7, 8 and 15 are Gern's type-D elements under the label map, Gern's Theorem 2.3.6 is re-derived for D5, D6 and D7 by exhaustive enumeration, and the layered palindromic words, orthogonal-reflection decompositions, the right weak order chain w_4 < w_6 < w_7 < w_8 and the length-additive factorizations through w_0(J), l(w_0(J)) = 6, of w_7 and w_8 are verified.
* The uniform verifier's exact polynomial identities and real-root witnesses, plus its separate finite matrix diagnostics; and the finite ingredients of the cover lemma for b_{r,k} (`uniform_family_checks.py`): for r = 3..8 and k = 0..3 the elements b s have the stated descent sets and are not fully commutative, s b s = r_{s beta}, l(b_{r,k}) = 8r^2 + 3 + 116k, and b_{3,0}, b_{3,1}, b_{4,0} have no fully commutative Bruhat cover.
* Only the proved affine conjugate family: exact translation and terminal identities, all finite-root affine inversion strings, FC catalogue closure, and the complete base cover exclusion. Its certificate `research/en_affine_referee/proved-affine-certificate.json` (the matrix of c_0 stored by rows, the column slopes, the inversion data and the eligible lower endpoint) is shipped and compared with the regenerated file.
* The E10 all-parameter reflection proof, including real-root lowering witnesses and the reduced base word.
* The complete E6 Kazhdan--Lusztig table with every leading coefficient (`e6_mu_table.cpp`): maximum mu = 10 over all 51,840 elements, attained by 8 pairs; 1,386 pairs with mu >= 2 (histogram 2:400, 3:556, 4:310, 5:108, 6:4, 10:8), none with a fully commutative lower endpoint; 6,431 pairs with a fully commutative lower endpoint and nonzero mu, all equal to 1; and the example with fully commutative upper endpoint and mu = 5.
* The affine E8 odd-gap polynomials (`e9_quotient_kl.cpp`): Kazhdan--Lusztig polynomials in the parabolic quotient by W_{R(w)} for the length-33 reflection r_beta, beta = (1,2,3,3,2,2,1,1,2), and its two length-34 extensions, with their eligible fully commutative lower endpoints (ideals of 5,826,496, 6,943,680 and 6,651,712 elements; 364,156, 216,990 and 207,866 cosets). The engine is first validated against a naive recursion on random elements of A4, D4 and E6.
* The Kazhdan--Lusztig polynomials P_{x_6,w_6} = 1 + 6q + 11q^2 + 6q^3 + q^4 + q^5 and P_{x_8,w_8} = P_{e,w_8} = 1 + 12q + 59q^2 + 154q^3 + 233q^4 + 221q^5 + 147q^6 + 70q^7 + 20q^8 + 2q^9 of Gern's elements (`d8_gern_kl.cpp`), in the signed-permutation model by default and also in the geometric model with `--full`.
* The affine D4 witness P_{x,xcx} = 1 + 3q + 2q^2, mu = 2, by direct subwords and R-polynomial reciprocity (`verify_affine_d4_r.py`).

The default verifies (rather than rediscovers) the supplied D6 recurrence dependency certificate. It recalculates the recurrence correction sets, Bruhat comparisons, all target descents, and the target R-polynomial reciprocity identity without invoking its KL evaluator. The target polynomial is `1+6q+11q^2+6q^3+q^4+q^5`; the same value is recomputed independently by `d8_gern_kl.cpp`.

The default also compares the newly computed finite classifications with the supplied E6/E7 full-group snapshots. It does not rerun those four searches unless `--full` is selected. The full mode regenerates both root-index snapshots (written by `enumerate_bad.cpp`) and both integer-matrix snapshots and requires each to equal the shipped file.

## Source independence and dependencies

Original computation sources are copied unchanged under `payload/`, with their original directory layout. The small `affine_proof.py` entry point calls only the proved part of the existing affine referee arithmetic; the original module's exploratory `main()` is never executed. No discovery searches or capped interval probes are included or run.

* `recursive_terminals.cpp` includes `parabolic_terminals.cpp`: those two algorithms share arithmetic, so they are not independent implementations.
* The D7 enumeration is a separate implementation and decomposition. The arbitrary-precision `verify_e8.py` imports neither enumeration engine. Its completeness evidence still includes rerunning the complete enumeration engines, not just reading their small output files.
* `verify_fc_catalogue.py` imports the row-matrix arithmetic in `verify_families.py`. The E9 closure verifier is independent of the C++ FC enumerator, but shares affine arithmetic with the affine verifier.
* `verify_indefinite_e10.py` imports `verify_indefinite.py` as an arithmetic library. The latter's E13 discovery-snapshot-reading `main()` is not run; its E13 input is not needed or shipped.
* The D6 certificate checker imports `SparseCoxeter` from `sparse_kl.py` and polynomial arithmetic from `coxeter.py`. It does not call its KL evaluator; its length, Bruhat, and R-polynomial primitives are shared with the generator, so it is not an independent implementation. The independent recomputation of the D6 value is `d8_gern_kl.cpp`.
* The uniform verifier imports neither the construction nor the group enumerators. `referee.txt` (an AI-generated internal note) and the manuscript explain the symbolic proof; arithmetic samples supplement it.
* The finite support-matching checker imports no enumeration engine. Its transcribed table words are additionally compared with the regenerated terminal matrices by the runner.
* The programs under `research/review_checks/` import nothing from the rest of the archive; each uses its own group arithmetic (integer matrices in the simple-root basis, signed permutations, or height vectors) and asserts the values stated in the manuscript before the runner compares its output with the shipped certificate.

`INPUTS.json` records every collected file, its original repository path, its purpose, and the source hash, including the hash of the manuscript snapshot, which the runner compares with the hash recorded by the uniform verifier. Two discovery-result JSON files are reduced to the exact seed rows used by the proofs; their original source hashes and selection rules are recorded. The large FC catalogues are generated, not shipped. The D6 dependency certificate is shipped because checking that certificate is the intended verification task. `MANIFEST.json` hashes every archive input other than itself. These checksums provide integrity and provenance within this package, not authentication or a substitute for reviewing the proofs.

The terminal engines and verifiers produce their detailed certificates in the isolated working tree. `expected-summary.json` records the stable counts, polynomials, all finite table support matchings, and family outcomes. No conclusion about all coefficients of arbitrary KL polynomials, any unresolved exploratory family, or historical priority is asserted.

## Build from the original working tree

The local repository's `supplement/build_package.py` collects an explicit inventory, records the current manuscript snapshot, and writes a deterministic ZIP and staging directory under `supplement/dist/`. Its output can be rebuilt with:

```sh
python3 supplement/build_package.py
```

A changed source or manuscript produces a changed archive and manifest. The ZIP uses sorted paths, fixed timestamps, and fixed permissions. Releases are hosted at `https://github.com/tygern/kl/releases`; GitHub hosting supplies a fixed release tag and downloadable assets, it is not a DOI deposit.
