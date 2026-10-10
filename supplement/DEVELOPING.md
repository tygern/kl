# Go module of the proof supplement

The Go module `github.com/tygern/kl/supplement` contains the enumeration engines, certificate verifiers, six review programs, proof runner and archive packager for `results/exceptional-leading.tex`. Go 1.22 or later suffices; the module uses only the standard library.

The Go programs replaced the Python and C++ programs in release v0.3.0. The earlier programs remain in release v0.2.0 and in git history (commit d997526); the repository has contained Go only since v0.4.0. Every certificate from the Go port was compared with the earlier programs' outputs.

Tracked files are sources only. Binaries are never tracked: `go build -C supplement -o bin/ ./cmd/...` writes them to `supplement/bin/` (ignored), and the proof runner builds its own copies into `runs/<mode>-<random>/bin/` inside an extracted archive.

## Build, test, lint

```sh
cd supplement
gofmt -l .          # must print nothing
go vet ./...
go build ./...
go test ./...       # small sanity tests, a few seconds
```

## Commands

Each directory under `cmd/` contains one program. Commands use `flag` for options and take no positional arguments. They read inputs relative to the current working directory, using the repository's layout, and write certificates to the original programs' paths or to stdout for capture by the runner. On a failed assertion, the command exits with a non-zero status. Certificate fields exclude timing, host, compiler and absolute paths.

| Command | Replaces | Invocation | Output |
|---|---|---|---|
| `cmd/fc-catalogue` | `research/en_independent/fc_catalogue.cpp` | `fc-catalogue -rank N` (N = 6..9) | stdout, captured to `research/en_independent/eN-fc.json` |
| `cmd/terminals-flat` | `research/en_e8/parabolic_terminals.cpp` | `terminals-flat -rank N` (N = 6, 7, 8) | stdout, captured to `research/en_e8/eN-validation.json` (N < 8) or `e8-terminals.json` |
| `cmd/terminals-recursive` | `research/en_e8/recursive_terminals.cpp` | `terminals-recursive -rank N`; `terminals-recursive -rank 8 -chains-certificate research/en_e8/e8-recursive-chains.json` runs the four E8 parabolic chains | stdout, captured to `research/en_e8/eN-recursive.json`; the chains certificate (snapshot compared) |
| `cmd/e8-d7` | `research/en_independent/e8_d7_cosets.cpp` | `e8-d7` | stdout, captured to `research/en_independent/e8-d7-terminals.json` |
| `cmd/enumerate-bad` | `research/broad_exceptional/enumerate_bad.cpp` | `enumerate-bad -rank N` (N = 6, 7; full mode) | stdout, captured to `research/broad_exceptional/eN_bad.json` |
| `cmd/matrix-search` | `research/verify_e6_independent.py`, `research/verify_e7_matrices.cpp` | `matrix-search -rank 6` writes `results/e6-independent-certificate.json`; `-rank 7` prints `results/e7-independent-certificate.json` (full mode) | snapshot compared |
| `cmd/verify-exceptional` | `research/verify_exceptional.py` | `verify-exceptional` | `results/exceptional-audit.json` |
| `cmd/verify-outputs` | `research/en_e8/verify_outputs.py` | `verify-outputs`; `verify-outputs -matrices N` prints the generated terminal matrices of `eN-recursive.json` | stdout |
| `cmd/verify-e8` | `research/en_independent/verify_e8.py` | `verify-e8` | `results/e8-independent-audit.json` |
| `cmd/finite-descents` | `research/ai-review-notes/check_finite_descents.py` | `finite-descents` | `research/ai-review-notes/finite-descents.json` |
| `cmd/d6-certificate` | `computations/verify_certificate.py` (+ `sparse_kl.py`, `coxeter.py`) | `d6-certificate` | JSON summary on stdout |
| `cmd/uniform-verify` | `research/en_uniform/referee_verify.py` | `uniform-verify` | `research/en_uniform/referee-certificate.json` (records the manuscript hash) |
| `cmd/affine-proof` | `supplement/affine_proof.py` | `affine-proof` | `research/en_affine_referee/proved-affine-certificate.json` |
| `cmd/affine-fc-covers` | `research/en_affine_referee/verify_fc_catalogue.py` | `affine-fc-covers` | `research/en_affine_referee/fc-cover-certificate.json` |
| `cmd/e10-all-k` | `research/en_affine_referee/verify_indefinite_e10.py` | `e10-all-k` | `research/en_affine_referee/e10-certificate.json` |
| `cmd/e6-mu-table` | `research/review_checks/e6_mu_table.cpp` | `e6-mu-table -type E -rank 6 -out results/e6-mu-table.json` | snapshot compared |
| `cmd/e9-quotient-kl` | `research/review_checks/e9_quotient_kl.cpp` | `e9-quotient-kl -mode certify -out results/e9-odd-gap-certificate.json` | snapshot compared |
| `cmd/d8-gern-kl` | `research/review_checks/d8_gern_kl.cpp` | `d8-gern-kl -out results/d8-gern-certificate.json -threads 0 -models sp` (`both` in full mode) | compared on `ranks` and `status` |
| `cmd/fc-maxima` | `research/review_checks/fc_maxima.cpp` | `fc-maxima -out results/fc-maxima-certificate.json -ranks 6,7,8,9,10,11,12,13` | snapshot compared |
| `cmd/terminal-structure` | `research/review_checks/terminal_structure.py` | `terminal-structure` | `results/terminal-structure-certificate.json` |
| `cmd/uniform-covers` | `research/review_checks/uniform_family_checks.py` | `uniform-covers [-full]` | `results/uniform-family-certificate.json` |
| `cmd/affine-d4` | `research/verify_affine_d4_r.py` | `affine-d4` | `results/affine-d4-independent.json` |
| `cmd/d6-ambient` | review code `kl.cpp`, `kl_d6.py` (history, d997526) | `d6-ambient -out results/d6-ambient-certificate.json` | snapshot compared |
| `cmd/terminal-data-check` | `research/exceptional_referee/verify_terminal_data.py` | `terminal-data-check` | stdout, captured to `research/exceptional_referee/checks.json` (snapshot compared) |
| `cmd/affine-reflection-family` | `research/en_families/affine_reflection_family.py` | `affine-reflection-family` | `research/en_families/affine_reflection_family.json` (snapshot compared) |
| `cmd/cartan-candidates` | `research/en_families/cartan_candidates.py` | `cartan-candidates -rank 10 -max-entry 2 -max-only` (the other committed variants: `-rank 10`, `-rank 11`, `-rank 12`, `-rank 13`, `-rank 13 -max-entry 2 -max-only`) | `research/en_families/cartan_E{N}_m{M}_max{0|1}.json` (E10 m2 max1 snapshot compared; `-out` overrides) |
| `cmd/uniform-construction` | `research/en_uniform/construction.py` | `uniform-construction -max-r 30` | `research/en_uniform/construction.json` (snapshot compared) |
| `cmd/proofs` | `supplement/run_proofs.py` | `go run ./cmd/proofs [-finite \| -full]` from an extracted archive | `runs/<mode>-<random>/{work,logs,bin,steps.json,summary.json,status.json}` |
| `cmd/package` | `supplement/build_package.py` | `go run -C supplement ./cmd/package` from the repository root | `supplement/dist/` (staging directory, zip, build-summary.json) |

Not ported: `research/broad_exceptional/exact_e.py` (historical provenance only; `verify-exceptional` records `enumerate-bad` as the snapshot generator). The `Replaces` column names the programs at their v0.2.0 paths; they are no longer in the tree (release v0.2.0 and commit d997526 have them).

Debugging-only flags (not used by the runner, never writing a certificate): `e6-mu-table -pairs FILE`, `e9-quotient-kl -mode pairs -pairs FILE` and `d8-gern-kl -pairs FILE` (with `-pairsgeo` to add the geometric model) print the Kazhdan--Lusztig polynomial of each pair listed in `FILE`, one per line. They exist for differential tests against the earlier engines (the port was tested on random pairs in A4, D4, D5, D6 and E6 with zero mismatches); `e9-quotient-kl -mode validate` is the original program's self-validation against a naive recursion.

## Shared code and separate implementations

Appendix A and the archive README list the following package imports and shared routines.

1. The root-index engine (`cmd/enumerate-bad`) and the integer-matrix engine (`cmd/matrix-search`) share no code: each imports no `internal` package.
2. The flat and recursive parabolic engines share arithmetic: `internal/parabolic`, imported by `cmd/terminals-flat` and `cmd/terminals-recursive` only. `cmd/e8-d7` and `cmd/verify-e8` are self-contained (the latter has its own integer-column arithmetic and inversion-count lengths).
3. `cmd/fc-catalogue` is self-contained. `internal/affine` is shared by `cmd/affine-proof` and `cmd/affine-fc-covers`. `internal/indefinite` is used by `cmd/e10-all-k`.
4. `cmd/d6-certificate` uses `internal/d6` for SparseCoxeter primitives and polynomial arithmetic. This library contains no KL evaluator.
5. `cmd/uniform-verify` is self-contained.
6. The six review programs (`cmd/e6-mu-table`, `cmd/e9-quotient-kl`, `cmd/d8-gern-kl`, `cmd/fc-maxima`, `cmd/terminal-structure`, `cmd/uniform-covers`) are each self-contained.
7. `cmd/verify-exceptional`, `cmd/verify-outputs`, `cmd/finite-descents`, `cmd/affine-d4`, `cmd/d6-ambient`, `cmd/terminal-data-check`, `cmd/affine-reflection-family`, `cmd/cartan-candidates`, `cmd/uniform-construction`, `cmd/proofs` and `cmd/package` are self-contained.

Small routines (Cartan matrices, reflection actions) are intentionally duplicated across self-contained commands. Where the C++ used threads (`d8-gern-kl`, `fc-maxima`) the Go programs use goroutines with a deterministic merge. All integer arithmetic is `int64` with overflow guards; exact rationals (`math/big`) replace Python fractions.

## Running verification

Build the archive, extract it outside the repository, then run the verification command from the archive root.

```sh
go run -C supplement ./cmd/package                      # from the repository root; writes supplement/dist/
unzip supplement/dist/exceptional-leading-proof.zip -d /some/scratch/dir
cd /some/scratch/dir/exceptional-leading-proof
go run ./cmd/proofs                             # default mode, about 30 s
go run ./cmd/proofs -finite                     # finite theorem only
go run ./cmd/proofs -full                       # full mode, about 90-100 s
```

The runner checks `MANIFEST.json`, copies `payload/` to `runs/<mode>-<random>/work/`, then builds the selected commands into `runs/<mode>-<random>/bin/` with `go build`. It executes the commands and compares regenerated snapshots with the shipped files, ignoring only `seconds` and `elapsed_seconds`. In default and full modes, it checks that the manuscript hash recorded by `uniform-verify` matches the one in `INPUTS.json` and that `summary.json` matches `expected-summary.json`. It also checks that the printed table words generate the regenerated terminal matrices (`verify-outputs -matrices N`). `status.json` records status, mode, seconds, `go_version`, step count, relative work path, the three boolean checks and the compared snapshots.

`-finite` and `-full` are mutually exclusive. Finite mode rebuilds the flat and recursive E6/E7/E8 terminal classifications, the independent E8/D7 classification, all four E8 parabolic chains, the printed descent/support matching data, and both the D6 recurrence check and the fresh ambient E7/E8/D6 polynomial calculation. It compares the regenerated terminal matrices with the printed table and with the saved full-group E6/E7 snapshots; only full mode reruns those historical full-group searches. Finite mode checks the manuscript hash directly against `INPUTS.json` and compares only the `finite`, `E8_D7`, `D6`, `finite_matchings`, `terminal_data_checks`, `E8_chains` and `D6_ambient` fields of `expected-summary.json`. Its saved summary contains only these recomputed outcomes. It builds no affine or indefinite-family program and does not run the additional complete E6 KL table.

With a `go.work` file at the repository root, `go run ./cmd/proofs` fails inside archives extracted below the repository, including the staging directory. The runner's `go build` steps use `GOWORK=off` and an empty `GOFLAGS`. For an archive below a `go.work` file, invoke the initial command as `GOWORK=off go run ./cmd/proofs [-finite | -full]`.

To regenerate a committed certificate, run its command from the repository root. The repository and runner's working directory have the same layout; the command rewrites the committed file in place.

```sh
go build -C supplement -o bin/ ./cmd/...                # binaries in supplement/bin/ (ignored by git)
supplement/bin/uniform-verify                           # rewrites research/en_uniform/referee-certificate.json
supplement/bin/matrix-search -rank 6                    # rewrites results/e6-independent-certificate.json
supplement/bin/fc-catalogue -rank 8 > research/en_independent/e8-fc.json   # stdout programs
```

`cmd/package` reads the archive README from `supplement/README.md`, the data files from their repository paths, and this module's sources; it computes `expected-summary.json` from the committed certificates with the same assertions as the runner and writes a deterministic zip (sorted paths, timestamps 2000-01-01, mode 0644, Deflate).
