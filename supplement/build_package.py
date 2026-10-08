#!/usr/bin/env python3
"""Collect an explicit, minimal proof supplement (release v0.2.0); no network or installations."""
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

VERSION='v0.2.0'
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
DIST=HERE/'dist'
STAGE=DIST/'exceptional-leading-proof'
# The internal review notes directory was renamed from research/journal-review
# to research/ai-review-notes; the archive always uses the new path.
RENAMED={'research/ai-review-notes/':'research/journal-review/'}
FILES={
 'output/pdf/exceptional-leading.pdf':'Rendered manuscript snapshot supplied alongside the TeX; not executed by the proof runner.',
 'results/exceptional-leading.tex':'Manuscript snapshot; the uniform verifier records its SHA256, which the runner compares with the value recorded here at build time.',
 'research/en_e8/parabolic_terminals.cpp':'Flat finite terminal engine; included as arithmetic library by the recursive engine.',
 'research/en_e8/recursive_terminals.cpp':'Recursive complete finite parabolic pruning engine.',
 'research/en_e8/verify_outputs.py':'Arbitrary-precision cross-comparison of all finite terminal outputs and baseline snapshots.',
 'research/en_independent/e8_d7_cosets.cpp':'Separate complete E8/D7 terminal enumeration.',
 'research/en_independent/fc_catalogue.cpp':'Generates complete FC catalogues of E6-E9; no catalogue snapshots shipped.',
 'research/en_independent/verify_e8.py':'Independent integer-matrix closure, inversion-length, descent and endpoint verifier; hashes JSON inputs without their timing fields.',
 'research/verify_exceptional.py':'Compares the four historical E6/E7 full-group snapshots (root-index and integer-matrix enumerations) and identifies the signed D6 pair.',
 'research/verify_e6_independent.py':'Full E6 independent integer-matrix enumeration (rerun with --full).',
 'research/verify_e7_matrices.cpp':'Full E7 independent integer-matrix enumeration (rerun with --full).',
 'research/broad_exceptional/enumerate_bad.cpp':'Full E6/E7 root-index enumeration (rerun with --full).',
 'research/broad_exceptional/exact_e.py':'Historical E6 search program whose hash is recorded by verify_exceptional; not executed by the runner.',
 'research/broad_exceptional/e6_bad.json':'Full E6 root-index snapshot written by enumerate_bad.cpp 6; regenerated and compared with --full.',
 'research/broad_exceptional/e7_bad.json':'Full E7 root-index snapshot written by enumerate_bad.cpp 7; regenerated and compared with --full.',
 'results/e6-independent-certificate.json':'Full E6 integer-matrix snapshot written by verify_e6_independent.py; regenerated and compared with --full.',
 'results/e7-independent-certificate.json':'Full E7 integer-matrix snapshot written by verify_e7_matrices.cpp; regenerated and compared with --full.',
 'research/en_uniform/referee_verify.py':'Independent symbolic uniform verifier and finite matrix diagnostics; records the manuscript hash.',
 'research/en_uniform/referee.txt':'AI-generated internal review note (not peer review) explaining the all-rank argument referenced by the uniform verifier.',
 'research/en_affine_referee/verify_families.py':'Row-matrix affine arithmetic and all-k inversion/terminal verifier; exploratory main not run.',
 'research/en_affine_referee/verify_fc_catalogue.py':'E9 arbitrary-precision FC closure and complete base cover verification.',
 'research/en_affine_referee/verify_indefinite.py':'Arithmetic library imported by the E10 verifier; E13 main not run.',
 'research/en_affine_referee/verify_indefinite_e10.py':'All-k E10 root/reflection proof verifier.',
 'research/en_affine_referee/proved-affine-certificate.json':'Shipped certificate of the proved affine family written by affine_proof.py (matrix of c_0 by rows, column slopes, inversion data, eligible lower endpoint); regenerated and compared.',
 'research/ai-review-notes/check_finite_descents.py':'Independent table descent and support matching verifier (from the AI-generated internal review notes).',
 'research/ai-review-notes/finite-descents.json':'Expected finite table matching certificate; regenerated and compared.',
 'computations/verify_certificate.py':'D6 saved recurrence certificate checker; shares length, Bruhat and R-polynomial primitives with the generator.',
 'computations/sparse_kl.py':'D6 signed permutation, length, Bruhat, lower ideal, and R-polynomial primitives.',
 'computations/coxeter.py':'Polynomial arithmetic imported by the D6 verifier and sparse primitives.',
 'results/d6-recurrence-certificate.json':'All 24,245 saved D6 dependency records, checked without calling a KL evaluator.',
 'research/review_checks/e6_mu_table.cpp':'Complete E6 Kazhdan-Lusztig table with every mu value: maximum 10 attained by 8 pairs, histogram of mu >= 2, and mu in {0,1} for every FC lower endpoint (manuscript Section 5).',
 'research/review_checks/e9_quotient_kl.cpp':'Parabolic-quotient KL engine for the length-33 affine E8 reflection r_beta and its two length-34 extensions; polynomials of the odd-gap eligible bottoms (Section 5); validated against a naive recursion.',
 'research/review_checks/d8_gern_kl.cpp':'KL polynomials on the lower ideals of Gern\'s w_6 and w_8 in two independent models of D_n: the D6 value 1+6q+11q^2+6q^3+q^4+q^5 and the D8 polynomial of the archived transfer note (geometric model with --full).',
 'research/review_checks/fc_maxima.cpp':'Independent FC enumeration of E6-E13 by height vectors: counts and maximum lengths (55, 66, 78, 92 for E10-E13; Section 4 remark).',
 'research/review_checks/terminal_structure.py':'Gern-plus-two matrix identities, exhaustive D5/D6/D7 enumeration, layered palindromes, orthogonal reflections, weak-order chain and w_0(J) factorizations of w_4, w_6, w_7, w_8 (Section 3 remarks).',
 'research/review_checks/uniform_family_checks.py':'Finite checks of the cover lemma ingredients for b_{r,k}: descents and non-FC status of bs, sbs = r_{s beta}, lengths 8r^2+3+116k, no FC Bruhat covers (Section 4).',
 'research/verify_affine_d4_r.py':'Independent affine D4 witness P_{x,xcx} = 1+3q+2q^2, mu = 2, by direct subwords and R-polynomial reciprocity (Section 5).',
 'results/affine-d4-independent.json':'Shipped affine D4 certificate; regenerated and compared.',
 'results/e6-mu-table.json':'Shipped E6 mu-table certificate; regenerated and compared.',
 'results/e9-odd-gap-certificate.json':'Shipped affine E8 odd-gap certificate; regenerated and compared.',
 'results/d8-gern-certificate.json':'Shipped D6/D8 Gern-element certificate (both models); regenerated and compared.',
 'results/fc-maxima-certificate.json':'Shipped FC maxima certificate for E6-E13; regenerated and compared.',
 'results/terminal-structure-certificate.json':'Shipped structure certificate for the exceptional terminals; regenerated and compared.',
 'results/uniform-family-certificate.json':'Shipped uniform-family certificate (default pairs); regenerated and compared.',
}


def source(path):
    """Repository file for an archive path, tolerating the pre-rename location of the review notes."""
    candidate=ROOT/path
    if candidate.is_file():return candidate
    for new,old in RENAMED.items():
        if path.startswith(new) and (ROOT/(old+path[len(new):])).is_file():
            print(f'note: {path} taken from the pre-rename location {old}')
            return ROOT/(old+path[len(new):])
    raise SystemExit('Missing source file: '+path)
def read(path):return json.loads(source(path).read_text())
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,value):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2)+'\n')


def expected():
    """Stable mathematical outcomes; must agree with the summary built by run_proofs.py."""
    result={'finite':{}}
    for n in (6,7,8):
        d=read(f'research/en_e8/e{n}-recursive.json')
        result['finite'][f'E{n}']={k:d[k] for k in ('group_order','right_terminal_count','commuting_terminals','noncommuting_terminal_count','fc_count')}
        result['finite'][f'E{n}']['bad_lengths']=sorted(row['length'] for row in d['bad'])
    d=read('research/en_independent/e8-d7-terminals.json')
    result['E8_D7']={k:d[k] for k in ('cosets','parabolic_right_terminals','candidates_tested','terminal_count')}
    result['D6']={'records':24245,'root_polynomial':[1,6,11,6,1,1],
                  'root_right_descent_checks':4,'root_R_reciprocity':True,'evaluator_called':False}
    result['uniform']={'symbolic_assertions':True,'checked_ranks':list(range(3,41))}
    result['affine']={'length':[27,92],'terminal_for_all_k':True}
    d=read('research/en_affine_referee/fc-cover-certificate.json')
    result['affine'].update({k:d[k] for k in ('FC_count','maximum_FC_length','base_Bruhat_cover_count','FC_base_covers')})
    d=read('research/en_affine_referee/e10-certificate.json')
    result['E10']={k:d[k] for k in ('rank','base_length','T_minus_identity_cube_zero','full_support_for_all_k','maximum_independent_size','eligible_mu')}
    result['finite_matchings']=read('research/ai-review-notes/finite-descents.json')
    d=read('results/e6-mu-table.json')
    result['E6_mu']={k:d[k] for k in ('group_order','fully_commutative_count','max_mu','mu_histogram','pairs_with_mu_at_least_2',
                                       'fully_commutative_lower_endpoint','fully_commutative_upper_endpoint','pairs_attaining_max_mu',
                                       'largest_mu_with_fully_commutative_upper_endpoint')}
    d=read('results/e9-odd-gap-certificate.json')
    result['E9_odd_gap']=[{k:e[k] for k in ('name','word_string','length','L','R','involution','terminal','quotient_size','lower_ideal_size','eligible_fully_commutative_bottoms')} for e in d['elements']]
    d=read('results/d8-gern-certificate.json')
    result['D6_D8_Gern']={rank:{k:v[k] for k in ('w_length','x_length','lower_ideal_size','interval_size','P_x_w_ascending','P_e_w_ascending','mu_x_w')} for rank,v in d['ranks'].items()}
    d=read('results/fc-maxima-certificate.json')
    result['FC_maxima']={rank:{k:v[k] for k in ('fully_commutative_count','maximum_length')} for rank,v in d['ranks'].items()}
    d=read('results/terminal-structure-certificate.json')
    result['terminal_structure']={'rows_not_from_type_D':{k:v['rows_not_from_type_D'] for k,v in d['gern_plus_two'].items() if isinstance(v,dict)},
                                  'type_D_noncommuting_terminals':{k:v['lengths'] for k,v in d['type_D_enumeration'].items()},
                                  'negated_roots':{k:v['negated_positive_roots'] for k,v in d['structure'].items()},
                                  'layers':{k:v['layers'] for k,v in d['structure'].items()},
                                  'chain_E8':[[r['lower'],r['upper'],r['prefix'] and r['suffix']] for r in d['right_weak_order_chain']['E8']['relations']],
                                  'inner_longest_parabolic_lengths':{k:v['max_inner_longest_parabolic_length'] for k,v in d['inner_longest_parabolic_factors'].items()},
                                  'w0_factorizations':{k:v['w0_factorization'] for k,v in d['inner_longest_parabolic_factors'].items() if 'w0_factorization' in v}}
    d=read('results/uniform-family-certificate.json')
    result['uniform_family']={'checked_pairs':d['checked_pairs'],'lengths':[r['length'] for r in d['pairs']],
                              'length_formula_holds':d['length_formula_holds_for_all_checked_pairs'],
                              'bruhat_covers':[[c['r'],c['k'],c['distinct_covers'],c['fully_commutative_covers']] for c in d['bruhat_covers']]}
    d=read('results/affine-d4-independent.json')
    result['affine_D4']={k:d[k] for k in ('lengths','lower_ideal_size','polynomial','mu')}
    # Values stated in the manuscript, asserted here so that a changed certificate cannot pass unnoticed.
    assert result['E6_mu']['max_mu']==10 and result['E6_mu']['mu_histogram']['10']==8 and result['E6_mu']['fully_commutative_lower_endpoint']=={'nonzero_mu_pairs':6431,'max_mu':1}
    assert [e['quotient_size'] for e in result['E9_odd_gap']]==[364156,216990,207866] and [e['lower_ideal_size'] for e in result['E9_odd_gap']]==[5826496,6943680,6651712]
    assert result['D6_D8_Gern']['D8']['P_x_w_ascending']==[1,12,59,154,233,221,147,70,20,2] and result['D6_D8_Gern']['D6']['P_x_w_ascending']==[1,6,11,6,1,1]
    assert {k:v['maximum_length'] for k,v in result['FC_maxima'].items() if k in ('E10','E11','E12','E13')}=={'E10':55,'E11':66,'E12':78,'E13':92}
    assert result['affine_D4']['polynomial']==[1,3,2] and result['affine_D4']['mu']==2
    for cert in ('results/e6-mu-table.json','results/e9-odd-gap-certificate.json','results/d8-gern-certificate.json','results/fc-maxima-certificate.json',
                 'results/terminal-structure-certificate.json','results/uniform-family-certificate.json'):
        assert read(cert)['status']=='passed',cert
    assert len(read('results/uniform-family-certificate.json')['pairs'])==24, 'ship the default-mode uniform certificate'
    return result


def main():
    # Only replace this builder's generated staging directory. Runs live in
    # extracted archives, so no proof-review logs are erased here.
    if STAGE.exists():shutil.rmtree(STAGE)
    STAGE.mkdir(parents=True)
    inventory=[]
    for archive_path,purpose in sorted(FILES.items()):
        origin=source(archive_path)
        target=STAGE/'payload'/archive_path;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(origin,target)
        inventory.append({'archive_path':str(target.relative_to(STAGE)),
                          'source_path':archive_path,'source_sha256':digest(origin),
                          'transformation':'unchanged','purpose':purpose})
    selections=[
      ('research/en_families/e9_full_support_eligible.json',
       lambda data:[next(row for row in data if row['length']==27)],
       'Retain only the length-27 base word row used by the proved affine family.'),
      ('research/en_families/cartan_E10_m2_max1.json',
       lambda data:{'real_terminal_roots':[next(row for row in data['real_terminal_roots'] if row['beta']==[3,7,10,9,7,6,4,3,1,6])]},
       'Retain only the E10 beta0 row and reduced-word witness used by verify_indefinite_e10.'),
    ]
    for path,select,purpose in selections:
        target=STAGE/'payload'/path;save(target,select(read(path)))
        inventory.append({'archive_path':str(target.relative_to(STAGE)),
                          'source_path':path,'source_sha256':digest(ROOT/path),
                          'transformation':purpose,'purpose':'Required fixed proof seed; discovery code and other results excluded.'})
    for name,destination in [('run_proofs.py','run_proofs.py'),('affine_proof.py','payload/affine_proof.py'),('README.md','README.md')]:
        target=STAGE/destination;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(HERE/name,target)
        inventory.append({'archive_path':destination,'source_path':'supplement/'+name,
                          'source_sha256':digest(HERE/name),'transformation':'unchanged',
                          'purpose':'Supplement packaging, proof-only entry point, or documentation.'})
    license_file=ROOT/'LICENSE'
    if not license_file.is_file():raise SystemExit('LICENSE (MIT) is missing from the repository root; it must be shipped in the archive.')
    shutil.copyfile(license_file,STAGE/'LICENSE')
    inventory.append({'archive_path':'LICENSE','source_path':'LICENSE','source_sha256':digest(license_file),
                      'transformation':'unchanged','purpose':'MIT License covering all code in this archive (copyright 2026 Tyson Gern).'})
    save(STAGE/'INPUTS.json',{'version':VERSION,'files':inventory,
                              'scope':f'Proof-only snapshot prepared for release {VERSION} of https://github.com/tygern/kl.',
                              'provenance':'Programs and notes were drafted by large language model systems (OpenAI GPT-6.1 Sol and GPT-6 Astra; review certificates adapted from an independent review by Anthropic Claude Fable 5.1) under the direction of the author, Tyson Gern, who takes responsibility for the content.'})
    save(STAGE/'expected-summary.json',expected())
    hashes={str(p.relative_to(STAGE)):digest(p) for p in sorted(STAGE.rglob('*')) if p.is_file()}
    save(STAGE/'MANIFEST.json',{'algorithm':'SHA256','scope':'Every distributed input file except this manifest itself.','sha256':hashes})
    archive=DIST/'exceptional-leading-proof.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for p in sorted(STAGE.rglob('*')):
            if not p.is_file():continue
            info=zipfile.ZipInfo('exceptional-leading-proof/'+str(p.relative_to(STAGE)),date_time=(2000,1,1,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            info.create_system=3;info.external_attr=0o100644<<16
            z.writestr(info,p.read_bytes(),compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)
    summary={'version':VERSION,'archive':str(archive.relative_to(ROOT)),'bytes':archive.stat().st_size,
             'sha256':digest(archive),'files':len(hashes)+1,'staging':str(STAGE.relative_to(ROOT)),
             'manuscript_sha256':digest(ROOT/'results/exceptional-leading.tex')}
    save(DIST/'build-summary.json',summary)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
