"""Inverse diagnostic reusing frozen scalar-mode data; no spectral recomputation.

Read the existing 512-mode spectrum and sum rule. Apply the reviewed
single-response heavy-residue bounds, with a positive omitted-tail allowance.
All numbers are benchmark estimates: stored roots/residues are floating point,
not validated interval enclosures. This is not a 2T-parent or survey result.
"""
from pathlib import Path
import csv,hashlib,json
ROOT = Path(__file__).resolve().parent
INPUT = ROOT/'finite_residue_5d_benchmark.json'
MODES = ROOT/'finite_residue_5d_spectrum.csv'
OUT = Path(__file__).with_name('gapped_residue_inverse_diagnostic_results.json')


def run():
    reference=json.loads(INPUT.read_text())
    with MODES.open() as stream:
        modes=list(csv.DictReader(stream))
    masses=[float(row['mass_squared']) for row in modes]
    weights=[float(row['normalized_weight']) for row in modes]
    assert len(modes)==reference['contract']['mode_count']
    assert all(m>0 for m in masses) and all(w>0 for w in weights)
    assert all(a<b for a,b in zip(masses,masses[1:]))
    missing=max(0.0,1-sum(weights))
    assert abs(missing-reference['sum_rule']['truncation_shortfall'])<1e-13
    light=masses[0];gap=reference['heavy_sector']['critical_gap_squared']
    assert light<gap<=masses[1]
    truth=reference['heavy_sector']['heavy_weight_exact']
    rows=[]
    for fraction in (0.01,0.05,0.20):
        q=fraction*gap
        partial=sum(w*q/(q+m) for w,m in zip(weights,masses))
        tail=missing*q/(q+masses[-1])
        x=q/(q+light);xh=q/(q+gap)
        base={'q_over_gap_squared':fraction,'q':q,'response_lower':partial,
              'response_upper':partial+tail,'omitted_response_allowance':tail,
              'light_response_unit_weight':x,'heavy_threshold_response':xh,
              'fitted_light_only_response':weights[0]*x,
              'heavy_response_upper_relative_to_light':(partial+tail-weights[0]*x)/(weights[0]*x)}
        for uncertainty in (0.0,0.005,0.01,0.03):
            lower_s=partial/(1+uncertainty)
            upper_s=(partial+tail)/(1-uncertainty)
            lo=max(0.0,1-upper_s/x)
            hi=min(1.0,(x-lower_s)/(x-xh))
            assert 0<=lo<=truth<=hi<=1
            rows.append(dict(base,normalization_relative_halfwidth=uncertainty,
                             residue_lower=lo,residue_upper=hi))
    out={'status':'numerically_verified inverse diagnostic using frozen scalar benchmark; not interval certified',
         'input_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (INPUT,MODES)},
         'contract':reference['contract'],'light_mass_squared':light,'heavy_gap_squared':gap,
         'true_heavy_weight':truth,'truncated_weight_shortfall':missing,'rows':rows,
         'limitations':['known light pole and declared heavy gap','normalization interval supplied externally to the inverse problem',
         'tail allowance treats frozen numerical modes as inputs; no root/residue rounding enclosure',
         'conditional scalar-response inference only; no covariant 2T completion, survey confidence or radiative protection']}
    OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'true_heavy_weight':truth,'rows':[{k:v for k,v in row.items() if k in ('q_over_gap_squared','normalization_relative_halfwidth','residue_lower','residue_upper','heavy_response_upper_relative_to_light')} for row in rows]},indent=2))


if __name__=='__main__':
    run()
