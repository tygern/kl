# Provenance

Tyson Gern directed this work in October 2026 and takes responsibility for
the manuscript. Large language models assisted with literature searches,
constructions, proof drafts, code, editing, and internal checks.

## Development history

- OpenAI's GPT-6.1 Sol and GPT-6 Astra assisted with the initial research,
  enumeration programs, constructions, and manuscript drafts.
- Anthropic's Claude Fable 5.1 produced an internal review on 7 October
  2026 and assisted with subsequent revisions and additional computations.
- Claude Sonnet 5.5 assisted with the Go port on 8 October; its outputs
  were compared against the earlier Python and C++ programs.
- The 9 October exposition revision expanded the central proofs and
  construction explanations, made the computational dependencies explicit,
  and shortened the review-history material. The stated main results and
  mathematical certificate values were retained.

- On 10 October, three GPT-6 Astra subagents independently reviewed finite
  correctness, infinite families, and novelty. The current revision fixes
  the FC-cover proof's opening quantifier, replaces the affine catalogue
  dependency by five braid witnesses and a translation factorization,
  prints compact terminal descriptions, and adds finite-only reproduction.
  The internal reports and fresh validation records are retained under
  `research/review-2026-10-10*`. Their audit scripts were preserved as
  independent Go commands under `tools/cmd/review-finite` and
  `tools/cmd/review-infinite`; these reviews are not peer review.

The original research prompts are retained in `prompts/`. Historical notes
and reviews are under `research/`; they describe earlier snapshots and are
not an alternative specification of the current results. Earlier Python,
C++ and review programs are available in release v0.2.0 and the repository
history at commit `d997526`. Release v0.4.0 preserves the manuscript before
the 9 October exposition revision. Release v0.5.0 contains that revision,
the matching PDF, and the proof-study companion with worked examples,
exercises, and solution sketches.

## What the checks establish

The current proof supplement regenerates the finite classifications and
computational inputs, compares outputs, and checks a stable mathematical
summary. Its README identifies shared code and separate implementations.
Appendix A of the manuscript distinguishes computations used in proofs
from additional checks and supplementary results.

Successful execution establishes reproducibility of those computations.
It does not by itself establish the correctness of every mathematical
argument, the completeness of the literature search, or historical
priority. AI-generated review statements are internal working notes, not
peer review or a record of the author's personal verification of each step.
The previous review's recommendations and categorical verdicts should be
read in that historical context.

## Snapshots and integrity

Build the supplement from the current checkout to include its manuscript
snapshot. `INPUTS.json` records the collected sources, and `MANIFEST.json`
hashes the archive contents. Each proof run keeps its logs, summary and
status. The PDF build manifest records the source and rendered PDF hashes.

Older manifests may refer to sources before later edits or the initial
history squash (`ab58ab2`). They are historical records and should not be
expected to match the current tree. Checksums identify snapshots and detect
changes; they do not authenticate mathematical claims.
