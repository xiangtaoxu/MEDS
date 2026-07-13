import netCDF4, numpy as np, glob, csv, sys

VARS = ['cas_temp_site','gpp_rate_fast','le_flux_fast','h_flux_fast','rnet_fast','soil_temp_top_site']
UMOL_TO_GC = 12.011e-6   # gC per umol CO2

def load_run(prefix):
    """Return dict (month,day,hour)->{var:val} concatenated over all daily FAST files."""
    files = sorted(glob.glob(f'out_integ/{prefix}-F-*.nc'))
    if not files: return None
    rec = {}
    for f in files:
        d = netCDF4.Dataset(f)
        mo = np.array(d.variables['month'][:]); dy = np.array(d.variables['day'][:]); hr = np.array(d.variables['hour'][:])
        cols = {v: np.array(d.variables[v][:]) for v in VARS if v in d.variables}
        for i in range(len(hr)):
            rec[(int(mo[i]),int(dy[i]),int(hr[i]))] = {v: cols[v][i] for v in cols}
        d.close()
    return rec

def metrics(ref, cand):
    keys = sorted(set(ref) & set(cand))
    out = {'n': len(keys)}
    for v in VARS:
        a = np.array([ref[k][v]  for k in keys if v in ref[k]])
        b = np.array([cand[k][v] for k in keys if v in cand[k]])
        if len(a)==0: continue
        diff = b - a
        out[v] = dict(rmse=float(np.sqrt(np.mean(diff**2))),
                      bias=float(np.mean(diff)),
                      maxerr=float(np.max(np.abs(diff))))
    # month-total GPP (umol/m2/s -> gC/m2 over the hours; each record is an hour mean)
    ga = np.array([ref[k]['gpp_rate_fast']  for k in keys])
    gb = np.array([cand[k]['gpp_rate_fast'] for k in keys])
    tot_ref = ga.sum()*3600*UMOL_TO_GC; tot_cand = gb.sum()*3600*UMOL_TO_GC
    out['gpp_month_gC'] = (tot_ref, tot_cand)
    out['gpp_month_bias_pct'] = (100*(tot_cand-tot_ref)/tot_ref) if abs(tot_ref)>1e-9 else float('nan')
    return out

# wall times
wall = {}
try:
    for r in csv.DictReader(open('integ/timings.csv')):
        wall[r['prefix']] = float(r['wall_sec'])
except FileNotFoundError:
    pass

RUN_TAGS = ['ref','ark_a_30','split_15','split_600','split_900','split_1800',
            'arkf_600','arkf_900','arkf_1800','arka_600','arka_900','arka_1800']

for month in ['jan','jun']:
    ref = load_run(f'i_{month}_ref')
    if ref is None: print(f"[{month}] REF missing"); continue
    print(f"\n{'='*118}\nMONTH = {month.upper()}   (REF = ARK adaptive @15s, {len(ref)} hourly records)\n{'='*118}")
    hdr = f"{'run':14s}{'wall_s':>8s} | {'casT_RMSE':>9s}{'soilT_RMSE':>11s} | {'GPP_RMSE':>9s}{'GPPmo_bias%':>12s} | {'LE_RMSE':>8s}{'H_RMSE':>8s}{'Rn_RMSE':>8s}"
    print(hdr); print('-'*len(hdr))
    for tag in RUN_TAGS:
        p = f'i_{month}_{tag}'
        cand = load_run(p)
        w = wall.get(p, float('nan'))
        if tag=='ref':
            print(f"{tag:14s}{w:8.1f} | {'(truth)':>9s}")
            continue
        if cand is None: print(f"{tag:14s}{w:8.1f} | MISSING"); continue
        m = metrics(ref, cand)
        def g(v,f): return m[v][f] if v in m else float('nan')
        print(f"{tag:14s}{w:8.1f} | {g('cas_temp_site','rmse'):9.4f}{g('soil_temp_top_site','rmse'):11.4f} | "
              f"{g('gpp_rate_fast','rmse'):9.4f}{m.get('gpp_month_bias_pct',float('nan')):12.3f} | "
              f"{g('le_flux_fast','rmse'):8.3f}{g('h_flux_fast','rmse'):8.3f}{g('rnet_fast','rmse'):8.3f}")
