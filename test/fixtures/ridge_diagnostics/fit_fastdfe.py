"""fastDFE GammaExp on the tgp SFS projected to 10 and 20 individuals.

Fits fastDFE (default p_b bound 0.5, and capped at GRAPES' 0.1), then evaluates the
fastDFE likelihood at the Julia and GRAPES GammaExpo optima (mapped: S_d = -negGmean,
b = negGshape, p_b = pos_prop, S_b = posGmean, eps = syno_orientation_error).
"""
import sys, csv, json
import numpy as np
import fastdfe as fd

d = sys.argv[1]
div = np.loadtxt(f"{d}/div.tsv")
Ln, Ls = div[2], div[3]
rows = list(csv.DictReader(open(f"{d}/params.tsv"), delimiter="\t"))
out = []
for k in (10, 20):
    s = np.loadtxt(f"{d}/sfs_k{k}.tsv")
    pn, ps = s[:, 1], s[:, 2]
    sfs_sel = fd.Spectrum([Ln - pn.sum(), *pn, 0.0])
    sfs_neut = fd.Spectrum([Ls - ps.sum(), *ps, 0.0])
    for cap in (0.5, 0.1):
        inf = fd.BaseInference(sfs_neut=sfs_neut, sfs_sel=sfs_sel, model=fd.GammaExpParametrization(),
                               n_runs=10, do_bootstrap=False, parallelize=False, seed=0,
                               bounds=dict(p_b=(0, cap)))
        inf.run()
        p = inf.params_mle
        rec = dict(k=k, source=f"fastdfe_cap{cap}", lnL_fastdfe=inf.likelihood, alpha_fastdfe=inf.alpha,
                   S_d=p["S_d"], b=p["b"], p_b=p["p_b"], S_b=p["S_b"], eps=p["eps"])
        if cap == 0.5:
            for r in rows:
                if int(r["k"]) != k:
                    continue
                q = {"all": dict(S_d=-float(r["negGmean"]), b=float(r["negGshape"]), p_b=float(r["pos_prop"]),
                                 S_b=float(r["posGmean"]), eps=float(r["syno_orientation_error"]))}
                out.append(dict(k=k, source=f"{r['source']}_params", lnL_fastdfe=inf.evaluate_likelihood(q),
                                alpha_fastdfe=inf.get_alpha(q["all"]), **{kk: vv for kk, vv in q["all"].items()}))
        out.append(rec)
json.dump(out, open(f"{d}/fastdfe_results.json", "w"), indent=1)
for o in out:
    print({kk: (round(vv, 5) if isinstance(vv, float) else vv) for kk, vv in o.items()})
