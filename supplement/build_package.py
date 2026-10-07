#!/usr/bin/env python3
"""Collect an explicit, minimal proof supplement; no network or installations."""
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
DIST=HERE/'dist'
STAGE=DIST/'exceptional-leading-proof'
FILES={
 'output/pdf/exceptional-leading.pdf':'Rendered manuscript snapshot supplied alongside the TeX; not executed by the proof runner.',
 'results/exceptional-leading.tex':'Manuscript snapshot; uniform verifier records its SHA256; contains all-rank proof.',
 'research/en_e8/parabolic_terminals.cpp':'Flat finite terminal engine; included as arithmetic library by recursive engine.',
 'research/en_e8/recursive_terminals.cpp':'Recursive complete finite parabolic pruning engine.',
 'research/en_e8/verify_outputs.py':'Arbitrary-precision cross-comparison of all finite terminal outputs and baseline snapshots.',
 'research/en_independent/e8_d7_cosets.cpp':'Separate complete E8/D7 terminal enumeration.',
 'research/en_independent/fc_catalogue.cpp':'Generates complete FC catalogues; no catalogue snapshots shipped.',
 'research/en_independent/verify_e8.py':'Independent integer-matrix closure, inversion-length, descent and endpoint verifier.',
 'research/verify_exceptional.py':'Historical E6/E7 snapshot comparison, signed D6 identification, and optional full-group driver.',
 'research/verify_e6_independent.py':'Optional full E6 independent integer-matrix enumeration.',
 'research/verify_e7_matrices.cpp':'Optional full E7 independent integer-matrix enumeration.',
 'research/broad_exceptional/enumerate_bad.cpp':'Optional full E6/E7 root-index enumeration.',
 'research/broad_exceptional/exact_e.py':'Historical source provenance hashed by verify_exceptional; no execution required in default/full runner.',
 'research/broad_exceptional/e6_bad.json':'Historical full E6 root-index snapshot; rebuilt with --full.',
 'research/broad_exceptional/e7_bad.json':'Historical full E7 root-index snapshot; rebuilt with --full.',
 'results/e6-independent-certificate.json':'Historical full E6 matrix snapshot; rebuilt with --full.',
 'results/e7-independent-certificate.json':'Historical full E7 matrix snapshot; rebuilt with --full.',
 'research/en_uniform/referee_verify.py':'Independent symbolic uniform verifier and finite matrix diagnostics.',
 'research/en_uniform/referee.txt':'Explanation of the all-rank proof referenced by the uniform verifier.',
 'research/en_affine_referee/verify_families.py':'Row-matrix affine arithmetic and all-k inversion/terminal verifier; exploratory main not run.',
 'research/en_affine_referee/verify_fc_catalogue.py':'E9 arbitrary-precision FC closure and complete base cover verification.',
 'research/en_affine_referee/verify_indefinite.py':'Arithmetic library imported by E10 verifier; E13 main not run.',
 'research/en_affine_referee/verify_indefinite_e10.py':'All-k E10 root/reflection proof verifier.',
 'research/journal-review/check_finite_descents.py':'Independent table descent and support matching verifier.',
 'research/journal-review/finite-descents.json':'Expected finite table matching certificate; regenerated and compared.',
 'computations/verify_certificate.py':'D6 saved recurrence certificate checker.',
 'computations/sparse_kl.py':'D6 signed permutation, length, Bruhat, lower ideal, and R-polynomial primitives.',
 'computations/coxeter.py':'Polynomial arithmetic imported by D6 verifier and sparse primitives.',
 'results/d6-recurrence-certificate.json':'All 24,245 saved D6 dependency records, checked without calling a KL evaluator.',
}


def read(path):return json.loads((ROOT/path).read_text())
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,value):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2)+'\n')


def expected():
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
    result['finite_matchings']=read('research/journal-review/finite-descents.json')
    return result


def main():
    # Only replace this builder's generated staging directory. Runs live in
    # extracted archives, so no proof-review logs are erased here.
    if STAGE.exists():shutil.rmtree(STAGE)
    STAGE.mkdir(parents=True)
    inventory=[]
    for source,purpose in sorted(FILES.items()):
        target=STAGE/'payload'/source;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(ROOT/source,target)
        inventory.append({'archive_path':str(target.relative_to(STAGE)),
                          'source_path':source,'source_sha256':digest(ROOT/source),
                          'transformation':'unchanged','purpose':purpose})
    selections=[
      ('research/en_families/e9_full_support_eligible.json',
       lambda data:[next(row for row in data if row['length']==27)],
       'Retain only the length-27 base word row used by the proved affine family.'),
      ('research/en_families/cartan_E10_m2_max1.json',
       lambda data:{'real_terminal_roots':[next(row for row in data['real_terminal_roots'] if row['beta']==[3,7,10,9,7,6,4,3,1,6])]},
       'Retain only the E10 beta0 row and reduced-word witness used by verify_indefinite_e10.'),
    ]
    for source,select,purpose in selections:
        target=STAGE/'payload'/source;save(target,select(read(source)))
        inventory.append({'archive_path':str(target.relative_to(STAGE)),
                          'source_path':source,'source_sha256':digest(ROOT/source),
                          'transformation':purpose,'purpose':'Required fixed proof seed; discovery code and other results excluded.'})
    for source,destination in [('run_proofs.py','run_proofs.py'),('affine_proof.py','payload/affine_proof.py'),('README.md','README.md')]:
        target=STAGE/destination;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(HERE/source,target)
        inventory.append({'archive_path':destination,'source_path':'supplement/'+source,
                          'source_sha256':digest(HERE/source),'transformation':'unchanged',
                          'purpose':'Supplement packaging, proof-only entry point, or documentation.'})
    # Output parents absent from a minimal input payload must be present before
    # original scripts write them. Empty directories are made by the runner's
    # generated FC output operations where needed.
    save(STAGE/'INPUTS.json',{'files':inventory,'scope':'Proof-only snapshot for https://github.com/tygern/kl/releases/tag/v0.1.0.'})
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
    summary={'archive':str(archive.relative_to(ROOT)),'bytes':archive.stat().st_size,
             'sha256':digest(archive),'files':len(hashes)+1,'staging':str(STAGE.relative_to(ROOT))}
    save(DIST/'build-summary.json',summary)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
