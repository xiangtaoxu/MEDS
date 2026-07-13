import re, os
tmpl = open('meds_config_census.toml').read()

MONTHS = {'jan': ('2024-01-01','2024-02-01'), 'jun': ('2024-06-01','2024-07-01')}
# (tag, integrator, adaptive, fixed_substep, dt_fast_s)
RUNS = [
    ('ref',        'ark',   'true',  1, 15),     # REF (truth)
    ('ark_a_30',   'ark',   'true',  1, 30),     # self-convergence check vs REF
    ('split_15',   'split', 'true',  1, 15),     # cross-scheme check vs REF
    ('split_600',  'split', 'true',  1, 600),
    ('split_900',  'split', 'true',  1, 900),
    ('split_1800', 'split', 'true',  1, 1800),
    ('arkf_600',   'ark',   'false', 1, 600),    # ARK fixed (1 coupled march/dt_fast)
    ('arkf_900',   'ark',   'false', 1, 900),
    ('arkf_1800',  'ark',   'false', 1, 1800),
    ('arka_600',   'ark',   'true',  1, 600),    # ARK adaptive
    ('arka_900',   'ark',   'true',  1, 900),
    ('arka_1800',  'ark',   'true',  1, 1800),
]

def gen(month, tag, integ, adaptive, fsub, dt):
    fis = 3600 // dt
    assert 3600 % dt == 0 and 86400 % dt == 0
    s = tmpl
    st, en = MONTHS[month]
    prefix = f"i_{month}_{tag}"
    s = s.replace('start_time = "2024-07-05"', f'start_time = "{st}"')
    s = s.replace('end_time   = "2024-07-06"', f'end_time   = "{en}"')
    s = re.sub(r'dt_fast\s*=\s*"[^"]*"', f'dt_fast             = "{dt}s"', s)
    s = re.sub(r'time_integrator\s*=\s*"[^"]*"', f'time_integrator     = "{integ}"', s)
    s = re.sub(r'ark_adaptive\s*=\s*\w+', f'ark_adaptive        = {adaptive}', s)
    # insert ark_fixed_substep after ark_adaptive line
    s = s.replace(f'ark_adaptive        = {adaptive}',
                  f'ark_adaptive        = {adaptive}\nark_fixed_substep   = {fsub}')
    s = re.sub(r'fast_interval_steps = \d+', f'fast_interval_steps = {fis}', s)
    s = re.sub(r'output_prefix\s*=\s*"[^"]*"', f'output_prefix         = "{prefix}"', s)
    s = re.sub(r'^prefix\s*=\s*"[^"]*"', f'prefix        = "{prefix}"', s, flags=re.M)
    s = re.sub(r'^dir\s*=\s*"out"', 'dir           = "out_integ"', s, flags=re.M)
    s = re.sub(r'write_state\s*=\s*true', 'write_state           = false', s)
    s = re.sub(r'cohort_max\s*=\s*4096', 'cohort_max    = 64', s)
    # daily ON, monthly/annual OFF
    s = s.replace('[output.monthly]\nenabled    = true', '[output.monthly]\nenabled    = false')
    open(f'integ/{prefix}.toml','w').write(s)
    return prefix

order=[]
for m in MONTHS:
    for r in RUNS:
        order.append(gen(m, *r))
open('integ/run_order.txt','w').write('\n'.join(order)+'\n')
print(f"generated {len(order)} configs into integ/")
