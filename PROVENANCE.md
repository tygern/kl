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
- The 9 October exposition revision contains longer central proofs and
  construction explanations, explicit computational dependency lists, and
  shorter review-history material. The stated main results and mathematical
  certificate values were retained.

- On 10 October, three GPT-6 Astra subagents independently reviewed finite
  correctness, infinite families, and novelty. The revisions correct
  the FC-cover proof's opening quantifier, replace the affine catalogue
  dependency with five braid witnesses and a translation factorization,
  include compact terminal descriptions, and add finite-only reproduction.
  The internal reports and fresh validation records are retained under
  `research/review-2026-10-10*`. Their audit scripts were preserved as
  independent Go commands under `tools/cmd/review-finite` and
  `tools/cmd/review-infinite`; these reviews are not peer review.

Release v0.6.1 adds two Dynkin diagrams and prose revisions throughout the
manuscript, supporting documentation, code comments, and diagnostic messages.
The style audit in `results/style-audit-2026-10-10.json` records file coverage,
retained historical wording, and inspection limits before release preparation.
The mathematical statements and computational results are unchanged.

The original research prompts are retained in `prompts/`. Historical notes
and reviews under `research/` describe earlier snapshots; the manuscript
and current supplement state the current results. Earlier Python,
C++ and review programs are available in release v0.2.0 and the repository
history at commit `d997526`. Release v0.4.0 preserves the manuscript before
the 9 October exposition revision. Release v0.5.0 contains that revision,
the matching PDF, and the proof-study companion with worked examples,
exercises, and solution sketches.

## Computational checks

The current proof supplement regenerates the finite classifications and
computational inputs, compares outputs, and checks a stable mathematical
summary. Its README identifies shared code and separate implementations.
Appendix A of the manuscript lists computations used in proofs separately
from additional checks and supplementary results.

A successful run reproduces those computations. Mathematical correctness,
literature coverage, and historical priority require separate review.
AI-generated review statements are internal working notes. They have not
undergone peer review and do not record the author's personal verification
of each step. The previous review's recommendations and categorical
verdicts concern the snapshot reviewed at that time.

## Snapshots and integrity

Build the supplement from the current checkout to include its manuscript
snapshot. `INPUTS.json` lists the collected sources, and `MANIFEST.json`
contains hashes of the archive contents. Each proof run keeps its logs,
summary and status. The PDF build manifest records the source and rendered
PDF hashes.

Older manifests may refer to sources before later edits or the initial
history squash (`ab58ab2`). Their hashes correspond to those historical
snapshots. Checksum comparisons detect file changes; mathematical claims
require proof review.
