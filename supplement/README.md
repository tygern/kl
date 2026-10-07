# Proof supplement: v0.1.0

This supplement reproduces the computational parts of `results/exceptional-leading.tex`. It accompanies [GitHub release v0.1.0](https://github.com/tygern/kl/releases/tag/v0.1.0). No distribution license, copyright ownership, or author identity is assigned by this packaging work.

## Run

Unzip the archive and enter its `exceptional-leading-proof` directory. Use Python 3.10 or later and an installed C++17 compiler:

```sh
python3 run_proofs.py
```

The runner searches for `clang++`, `c++`, or `g++`; a particular installed executable can be selected with `--compiler /absolute/path/to/compiler`. No network, installation, third-party Python packages, or `uv` are required. Do not use Python's `-O` option: assertions are part of verification. Compiler assertions are also retained.

For optional historical full-group searches:

```sh
python3 run_proofs.py --full
```

The full mode adds the complete E6 and E7 root-index and integer-matrix searches, including all 2,903,040 E7 elements. It uses substantially more memory than the default. Both modes retain their isolated working trees, binaries, per-step stdout/stderr, stable summary, and status under `runs/`. The initial `payload/` stays unchanged; repeated runs are independent. An interrupted or failed run is retained for diagnosis. The runner first checks the SHA256 inventory, then requires the stable proof outcomes to equal `expected-summary.json`. Timing, compiler paths, and generated hashes of mutable output files are not treated as stable mathematical output.

## What the default recomputes

* E6/E7/E8 complete terminal classifications by flat and recursive parabolic pruning, and the separate E8/D7 classification. It checks completeness through the exhaustive source algorithms and compares all terminal matrices; a `complete: true` field alone is not accepted as an exhaustiveness argument.
* E6/E7/E8/E9 FC catalogues, including independent E8 closure and E9 closure verification. Finite catalogues are redundant checks for the revised finite theorem. The E9 catalogue and its maximum length are used in the stated affine cover proof.
* Descent sets and support matchings for all eleven printed finite table rows, and equality of their matrices with the regenerated terminal list.
* The uniform verifier's exact polynomial identities and real-root witnesses, plus its separate finite matrix diagnostics. Its all-rank argument is the symbolic argument and the manuscript proof, not extrapolation from sampled ranks.
* Only the proved affine conjugate family: exact translation and terminal identities, all finite-root affine inversion strings, FC catalogue closure, and the complete base cover exclusion.
* The E10 all-parameter reflection proof, including real-root lowering witnesses and the reduced base word.

The default verifies (rather than rediscovers) the supplied D6 recurrence dependency certificate. It recalculates the recurrence correction sets, Bruhat comparisons, all target descents, and the target R-polynomial reciprocity identity without invoking its KL evaluator. The target polynomial is `1+6q+11q^2+6q^3+q^4+q^5`.

The default also compares the newly computed finite classifications with supplied historical E6/E7 full-group snapshots. It does not rerun those four historical searches unless `--full` is selected. The full mode regenerates both the E6 and E7 root-index snapshots and both independent matrix snapshots before the same comparison. This corrects the ambiguity of calling snapshot comparison a fresh full-group search.

## Source independence and dependencies

Original computation sources are copied unchanged under `payload/`, with their original directory layout. The small new `affine_proof.py` entry point calls only the proved part of the existing affine referee arithmetic; the original module's exploratory `main()` is never executed. No discovery searches or capped interval probes are included or run.

* `recursive_terminals.cpp` includes `parabolic_terminals.cpp`: those two algorithms share arithmetic, so they are not independent implementations.
* The D7 enumeration is a separate implementation and decomposition. The arbitrary-precision `verify_e8.py` imports neither enumeration engine. Its completeness evidence still includes rerunning the complete enumeration engines, not just reading their small output files.
* `verify_fc_catalogue.py` imports the row-matrix arithmetic in `verify_families.py`. The E9 closure verifier is independent of the C++ FC enumerator, but shares affine arithmetic with the affine verifier.
* `verify_indefinite_e10.py` imports `verify_indefinite.py` as an arithmetic library. The latter's E13 discovery-snapshot-reading `main()` is not run; its E13 input is not needed or shipped.
* The D6 certificate checker imports `SparseCoxeter` from `sparse_kl.py` and polynomial arithmetic from `coxeter.py`. It does not call its KL evaluator; its length, Bruhat, and R-polynomial primitives are shared code, so this is not a claim of wholly independent implementations of every component.
* The uniform verifier imports neither the construction nor the group enumerators. `referee.txt` and the manuscript explain the symbolic proof; arithmetic samples supplement it.
* The finite support-matching checker imports no enumeration engine. Its transcribed table words are additionally compared with the regenerated terminal matrices by the runner.

`INPUTS.json` records every collected file, its original repository path, its purpose, and the source hash. Two discovery-result JSON files are reduced to the exact seed rows used by the proofs; their original source hashes and selection rules are recorded. The large FC catalogues are generated, not shipped. The D6 dependency certificate is shipped because checking that certificate is the intended verification task. `MANIFEST.json` hashes every archive input other than itself. These checksums provide integrity and provenance within this package, not authentication or a substitute for reviewing the proofs.

The terminal engines and verifiers produce their detailed certificates in the isolated working tree. `expected-summary.json` records the stable counts, target polynomial, all finite table support matchings, and family outcomes. No conclusion about all coefficients of arbitrary KL polynomials, any unresolved exploratory family, or historical priority is asserted.

## Build from the original working tree

The local repository's `supplement/build_package.py` collects an explicit inventory, records the current manuscript snapshot, and writes a deterministic ZIP and staging directory under `supplement/dist/`. Its output can be rebuilt with:

```sh
python3 supplement/build_package.py
```

A changed source or manuscript produces a changed archive and manifest. The ZIP uses sorted paths, fixed timestamps, and fixed permissions. The versioned release is hosted at `https://github.com/tygern/kl/releases/tag/v0.1.0`. GitHub hosting supplies a fixed release tag and downloadable assets; it is not a DOI deposit.
