# Proof supplement v0.6.1

Author: Tyson Gern (Initial Capacity, tyson@initialcapacity.io, ORCID 0009-0003-8288-8786).

The proof supplement reproduces the computational parts of `results/exceptional-leading.tex`. Release v0.6.1 includes the revised prose, two Dynkin diagrams, a matching PDF, the proof-study companion, and the citation and style audits. The affine braid-witness proof, exceptional-terminal descriptions, and finite-only verification mode introduced in v0.6.0 are retained. The style audit is at `payload/results/style-audit-2026-10-10.json`; its file hashes record the checkout at the end of the prose audit, before release preparation. The two companion documents are at `payload/results/proof-study-companion.md` and `payload/results/citation-audit.md` in the archive (`results/` in the repository). Release v0.5.1 preserves the preceding citation revision. The stated theorems and previously recorded mathematical values are unchanged. The affine verifier checks five braid witnesses and their occurrence in reduced words for all parameters, without a catalogue input.

The author directed development of the programs and notes with assistance from OpenAI and Anthropic models. AI-generated review notes are not peer review. The computational dependencies and shared routines are listed below; the repository's `PROVENANCE.md` records the development history. Code is licensed under MIT and the manuscript under CC BY 4.0.

## Run

Unzip the archive and enter its `exceptional-leading-proof` directory. The only requirement is Go 1.22 or later (`go version`); the module uses the standard library alone, so no network access, module download, third-party package or compiler other than Go is needed.

```sh
go run ./cmd/proofs
```

The runner checks the SHA256 inventory in `MANIFEST.json`, copies `payload/` to an isolated working tree, then builds the selected programs with `go build` into the run's `bin/` directory. It executes each program inside the working tree and compares regenerated snapshots with the shipped files, ignoring only the timing keys `seconds` and `elapsed_seconds`. It then compares the computed mathematical summary with `expected-summary.json`. The default run takes about 30 s on a 2026 laptop (10 cores, otherwise idle; longer under load) and uses less than 2 GB of memory. If the archive is extracted below a directory containing a `go.work` file, run `GOWORK=off go run ./cmd/proofs` (the runner's `go build` steps already ignore any enclosing workspace).

```sh
go run ./cmd/proofs -full
```

Full mode also runs the complete E6 and E7 root-index and integer-matrix searches (all 2,903,040 E7 elements), the second (geometric) model of the D6/D8 Kazhdan--Lusztig computation, and the larger uniform-family checks (ranks up to E49 and more Bruhat cover sets). The runner compares every regenerated snapshot with the shipped file. Full mode takes about 90-100 s on the same machine, most of it in the D6/D8 Kazhdan--Lusztig computation, and uses less than 2 GB of memory.

Each run retains its isolated working tree (`work/`), binaries (`bin/`), per-step stdout and stderr (`logs/`), step record with timings (`steps.json`), mathematical summary (`summary.json`) and final status (`status.json`) under `runs/<mode>-<random>/`. The initial `payload/` stays unchanged; repeated runs are independent. Failed and interrupted runs retain their working files and logs for diagnosis. Summary comparisons exclude timing, the Go version and generated-file paths. Recorded paths are relative to the archive root.

For the finite theorem alone:

```sh
go run ./cmd/proofs -finite
```

In finite mode, the runner regenerates the flat and recursive E6/E7/E8 terminal lists, the independent E8/D7 list and four E8 parabolic chains. It compares their outputs with the shipped E6/E7 full-group snapshots, checks printed words, descents and support matchings, verifies the saved D6 recurrence, and recomputes that polynomial in D6, E7 and E8. It checks the archive manifest and manuscript hash, then compares the finite part of the expected summary. FC catalogues, infinite-family checks and supplementary KL tables are omitted. Use `-full` for fresh full-group E6/E7 cross-checks; `-finite` and `-full` cannot be combined. All three modes retain separate logs and records.

`go vet ./...` and `go test ./...` run the module's static checks and its small unit tests (a few seconds); they are not part of the proof run.

## What the default recomputes

* E6/E7/E8 complete terminal classifications by flat and recursive parabolic pruning (`cmd/terminals-flat`, `cmd/terminals-recursive`), and the separate E8/D7 classification (`cmd/e8-d7`). The runner executes the exhaustive algorithms and compares all terminal matrices (`cmd/verify-outputs`). Exhaustiveness follows from the algorithms; the `complete: true` field records their completion.
* E6/E7/E8/E9 FC catalogues (`cmd/fc-catalogue`), including independent E8 closure (`cmd/verify-e8`) and E9 closure verification (`cmd/affine-fc-covers`). The E9 catalogue and its maximum length 44 are supplementary checks. The affine cover proof uses five braid witnesses. The independent program `cmd/fc-maxima` enumerates FC elements by height vectors and recomputes the counts and maximum lengths of E6 to E13 (E10 to E13: 180,438/55, 737,762/66, 3,021,000/78, 12,387,990/92).
* Descent sets and support matchings for all eleven printed finite table rows (`cmd/finite-descents`), and equality of their matrices with the regenerated terminal list (the runner compares the table words with the matrices printed by `verify-outputs -matrices N`).
* The structure of the exceptional terminals (`cmd/terminal-structure`). The program checks that the rows of Table 1 of lengths 7, 8 and 15 are Gern's type-D elements under the label map and re-derives Gern's Theorem 2.3.6 for D5, D6 and D7 by exhaustive enumeration. It verifies the layered palindromic words, orthogonal-reflection decompositions, the right weak order chain w_4 < w_6 < w_7^E < w_8^E and the length-additive factorizations through w_0(J), l(w_0(J)) = 6, of the exceptional w_7^E and w_8^E.
* The uniform verifier's exact polynomial identities and real-root witnesses, plus its separate finite matrix diagnostics (`cmd/uniform-verify`); and the finite ingredients of the cover lemma for b_{r,k} (`cmd/uniform-covers`): for r = 3..8 and k = 0..3 the elements b s have the stated descent sets and are not fully commutative, s b s = r_{s beta}, l(b_{r,k}) = 8r^2 + 3 + 116k, and b_{3,0}, b_{3,1}, b_{4,0} have no fully commutative Bruhat cover.
* The proved affine conjugate family (`cmd/affine-proof`): exact translation and terminal identities, all finite-root affine inversion strings, the length-additive translation factorization and five reduced braid witnesses. This command reads neither the FC catalogue nor discovery data. Its certificate `research/en_affine_referee/proved-affine-certificate.json` records the base matrix by rows, column slopes, inversion data, translation lengths and braid witnesses; it is shipped and compared with the regenerated file. FC catalogue closure and all 21 base-cover exclusions (`cmd/affine-fc-covers`) remain additional corroboration.
* The E10 all-parameter reflection proof (`cmd/e10-all-k`), including real-root lowering witnesses and the reduced base word.
* The complete E6 Kazhdan--Lusztig table with every leading coefficient (`cmd/e6-mu-table`): maximum mu = 10 over all 51,840 elements, attained by 8 pairs; 1,386 pairs with mu >= 2 (histogram 2:400, 3:556, 4:310, 5:108, 6:4, 10:8), none with a fully commutative lower endpoint; 6,431 pairs with a fully commutative lower endpoint and nonzero mu, all equal to 1; and the example with fully commutative upper endpoint and mu = 5.
* The affine E8 odd-gap polynomials (`cmd/e9-quotient-kl`): Kazhdan--Lusztig polynomials in the parabolic quotient by W_{R(w)} for the length-33 reflection r_beta, beta = (1,2,3,3,2,2,1,1,2), and its two length-34 extensions, with their eligible fully commutative lower endpoints (ideals of 5,826,496, 6,943,680 and 6,651,712 elements; 364,156, 216,990 and 207,866 cosets). The engine is first validated against a naive recursion on random elements of A4, D4 and E6.
* The Kazhdan--Lusztig polynomials P_{x_6,w_6} = 1 + 6q + 11q^2 + 6q^3 + q^4 + q^5 and P_{x_8^D,w_8^D} = P_{e,w_8^D} = 1 + 12q + 59q^2 + 154q^3 + 233q^4 + 221q^5 + 147q^6 + 70q^7 + 20q^8 + 2q^9 of Gern's elements (`cmd/d8-gern-kl`), in the signed-permutation model by default and also in the geometric model with `-full`.
* The affine D4 witness P_{x,xcx} = 1 + 3q + 2q^2, mu = 2, by direct subwords and R-polynomial reciprocity (`cmd/affine-d4`).
* The E8 classification repeated along four chains of parabolic subgroups (`cmd/terminals-recursive -chains-certificate`): the chain through E7 used in the finite proof, chains through D7 and A7, and a second ordering through E7 (adding generators 0,1,7,2,3,4,5,6). Along each chain, the program enumerates the same 2,160 right-terminal and 64 terminal elements. The certificate records every chain's coset and candidate counts.
* The Kazhdan--Lusztig polynomial of the length-15 pair of Table 1 computed from scratch by the KL recursion on the 3,184-element lower ideal of b in ambient E7 (b = 132543621324356, x = 1356), in ambient E8 (b = 132543721324357, x = 1357) and in the D6 parabolic (`cmd/d6-ambient`): 1 + 6q + 11q^2 + 6q^3 + q^4 + q^5, mu = 1, interval of 1,676 elements with the stated rank vector.
* A third check of the Table 1 data from the E6/E7 integer-matrix certificates (`cmd/terminal-data-check`): the Poincare polynomials against the length distributions, the commuting terminal counts 22 and 36 as independent sets of the diagrams, the transcription of the four E7 bad terminal words and the identification of the D6 pair with the signed-permutation model.
* The affine E8 reflection family r_{beta+k delta}, beta = (1,2,3,3,2,2,1,1,2) (`cmd/affine-reflection-family`): the translation T with T(beta) = beta + delta, the inversion count 33 + 58k from the 240 finite roots, and for k = 0..10 that the reflection is terminal with full support and that both one-generator extensions of length 34 + 58k have a full-support fully commutative bottom element.
* The E10 seed of the all-k reflection proof regenerated by exhausting the terminal reflections of E10 whose root pairings on a maximum independent set are bounded by 2 (`cmd/cartan-candidates -rank 10 -max-entry 2 -max-only`): eight real terminal roots, each certified by reduction to a simple root, including the root beta_0 that `cmd/e10-all-k` then reads from the regenerated file.
* The finite counterchecks of the rank-uniform E_{4r+1} construction for r = 3..30 (`cmd/uniform-construction -max-r 30`): the signed-coordinate seed certificate, the translation identities, and for k = 0, 1, 2, 10 that beta_{r,k} is a norm-2 positive root whose reflection has the stated right descent set and is right terminal.

In default mode, `cmd/d6-certificate` verifies the supplied D6 recurrence dependency certificate. It recalculates the recurrence correction sets, Bruhat comparisons, all target descents, and the target R-polynomial reciprocity identity. Its arithmetic library contains no KL evaluator. The target polynomial is `1+6q+11q^2+6q^3+q^4+q^5`; `cmd/d8-gern-kl` recomputes the same value independently.

The default also compares the newly computed finite classifications with the supplied E6/E7 full-group snapshots (`cmd/verify-exceptional`). It does not rerun those four searches unless `-full` is selected. The full mode regenerates both root-index snapshots (`cmd/enumerate-bad`) and both integer-matrix snapshots (`cmd/matrix-search`) and requires each to equal the shipped file.

## Source independence and dependencies

The programs are the packages under `cmd/`; shared arithmetic is in the four packages under `internal/`. Every `cmd/` package imports the Go standard library and at most one `internal/` package; `cmd/proofs`, the runner, imports none of them. The data files are copied unchanged under `payload/`, with their original repository layout. No discovery searches or capped interval probes are included or run.

* `cmd/terminals-flat` and `cmd/terminals-recursive` both import `internal/parabolic`: those two algorithms share arithmetic, so they are not independent implementations (the original recursive engine included the flat engine's source).
* `cmd/e8-d7` uses a separate implementation and parabolic decomposition and imports no internal package. The exact integer verifier `cmd/verify-e8` imports neither enumeration engine. The proof runner reruns the complete enumeration engines before the output-file checks.
* `cmd/fc-catalogue` is self-contained. `cmd/affine-fc-covers` imports `internal/affine`, as does `cmd/affine-proof`: the E9 closure verifier is independent of the FC enumerator, but shares affine arithmetic with the affine verifier.
* `cmd/e10-all-k` imports `internal/indefinite` as an arithmetic library (the original library's E13 discovery-snapshot-reading main is not ported; its E13 input is not needed or shipped).
* `cmd/d6-certificate` imports `internal/d6` for its signed-permutation, length, Bruhat, lower-ideal, R-polynomial and polynomial primitives. The checker and the certificate's generator share these primitives; `internal/d6` contains no KL evaluator. The independent recomputation of the D6 value is `cmd/d8-gern-kl`.
* `cmd/uniform-verify` imports neither the construction nor the group enumerators. `referee.txt` (an AI-generated internal note) and the manuscript explain the symbolic proof; arithmetic samples supplement it.
* `cmd/finite-descents` imports no enumeration engine. Its transcribed table words are additionally compared with the regenerated terminal matrices by the runner.
* The six review programs `cmd/e6-mu-table`, `cmd/e9-quotient-kl`, `cmd/d8-gern-kl`, `cmd/fc-maxima`, `cmd/terminal-structure` and `cmd/uniform-covers` import nothing from the rest of the module; each uses its own group arithmetic (integer matrices in the simple-root basis, signed permutations, or height vectors) and asserts the values stated in the manuscript before the runner compares its output with the shipped certificate.
* `cmd/verify-exceptional`, `cmd/verify-outputs`, `cmd/enumerate-bad`, `cmd/matrix-search` and `cmd/affine-d4` are likewise self-contained, as are the five programs added in v0.4.0: `cmd/d6-ambient`, `cmd/terminal-data-check`, `cmd/affine-reflection-family`, `cmd/cartan-candidates` and `cmd/uniform-construction` (the four-chain E8 certificate is written by `cmd/terminals-recursive` itself). Small routines such as Cartan matrices and reflection actions are intentionally duplicated across self-contained programs.

Where the earlier programs used threads (`cmd/d8-gern-kl`, `cmd/fc-maxima`), the Go programs use goroutines with a deterministic merge; their output does not depend on scheduling. All arithmetic is exact: 64-bit integers with overflow guards, or arbitrary-precision rationals where the originals used Python fractions. Three of the review programs also accept debugging-only flags that were used for the differential tests of the port against the earlier programs (`e6-mu-table -pairs FILE`, `e9-quotient-kl -mode pairs -pairs FILE`, `d8-gern-kl -pairs FILE` with optional `-pairsgeo`): they print the Kazhdan--Lusztig polynomials of the pairs listed in the file and never write a certificate; the runner does not use them.

`INPUTS.json` records every collected file, its original repository path, its purpose, and its source hash, including the manuscript snapshot's hash. The runner compares that hash with the one recorded by the uniform verifier. One discovery-result JSON file (the E9 eligible rows) contains only the base row used by the supplementary cover check; its original source hash and selection rule are recorded. The E10 seed is shipped complete. The runner regenerates it with `cmd/cartan-candidates` and compares it before `cmd/e10-all-k` reads it. The runner generates the large FC catalogues. It verifies the shipped D6 dependency certificate. `MANIFEST.json` contains hashes of every other archive file, including every Go source file. Checksum comparisons verify file integrity within the package. Source authentication and proof review are separate tasks.

The terminal engines and verifiers write detailed certificates in the isolated working tree. `expected-summary.json` records the counts, polynomials, finite table support matchings, and family outcomes. The earlier mathematical values are retained, and v0.6.0 adds the affine braid count and all-parameter cover-exclusion checks. The `-finite` summary contains only finite-theorem inputs. Results for all coefficients of arbitrary KL polynomials, unresolved exploratory families, and historical priority are outside the verification's scope.

## Build from the original working tree

The local repository's `supplement/cmd/package` collects an explicit inventory, records the current manuscript snapshot, and writes a deterministic ZIP and staging directory under `supplement/dist/`. Its output can be rebuilt from the repository root with:

```sh
go run -C supplement ./cmd/package
```

The Go module is in the repository's `supplement/` directory; `-C supplement` runs the command there. The packager writes the archive under `supplement/dist/`.

After a source or manuscript change, the rebuilt archive and manifest differ. The ZIP uses sorted paths, fixed timestamps (2000-01-01) and fixed permissions (0644). The earlier Python and C++ programs remain in release v0.2.0 and in the repository history (commit d997526); the repository itself contains Go only. Releases at `https://github.com/tygern/kl/releases` have fixed tags and downloadable assets. They have no DOI deposit.
