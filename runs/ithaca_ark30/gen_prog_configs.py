#!/usr/bin/env python3
"""Generate fully-prognostic (leaf+wood) SPLIT/Picard configs at dt_fast=30s for
Jan (cold/snow) and Jun (hot). Output = the FAST netCDF tier at 15-min cadence
(fast_interval_steps = 30 * 30s = 900s); air temp + per-cohort leaf temp + height
ride the tier, so tallest-cohort selection is pure post-processing. Base = census."""
import re

tmpl = open('meds_config_census.toml').read()
MONTHS = {'jan': ('2024-01-01', '2024-02-01'), 'jun': ('2024-06-01', '2024-07-01')}

def gen(month):
    st, en = MONTHS[month]
    s = tmpl
    s = s.replace('start_time = "2024-07-05"', f'start_time = "{st}"')
    s = s.replace('end_time   = "2024-07-06"', f'end_time   = "{en}"')
    s = re.sub(r'dt_fast\s*=\s*"[^"]*"', 'dt_fast             = "30s"', s)
    s = re.sub(r'time_integrator\s*=\s*"[^"]*"', 'time_integrator     = "split"', s)
    s = re.sub(r'integration_scheme\s*=\s*"[^"]*"', 'integration_scheme  = "picard"', s)
    # fully prognostic leaf + wood
    s = re.sub(r'(fast_probe_file\s*=\s*"[^"]*")',
               r'\1\nleaf_energy_model   = "prognostic"\nwood_energy_model   = "prognostic"', s)
    s = re.sub(r'write_state\s*=\s*true', 'write_state           = false', s)
    # FAST netCDF output tier @ 15 min (30 * 30s); energy group carries air/cas/soil/leaf/height.
    s = re.sub(r'fast_interval_steps\s*=\s*\d+', 'fast_interval_steps = 30', s)
    s = re.sub(r'^carbon_fluxes\s*=\s*true', 'carbon_fluxes = false', s, flags=re.M)
    s = re.sub(r'^water_fluxes\s*=\s*true',  'water_fluxes  = false', s, flags=re.M)
    s = re.sub(r'^dir\s*=\s*"out"', 'dir           = "out_nc"', s, flags=re.M)
    s = re.sub(r'^\s*prefix\s*=\s*"[^"]*"', f'prefix        = "prog_{month}"', s, flags=re.M)
    s = re.sub(r'output_prefix\s*=\s*"[^"]*"', f'output_prefix         = "prog_{month}"', s)
    # only the FAST stream: disable daily / monthly (annual already off)
    s = s.replace('[output.daily]\nenabled    = true',   '[output.daily]\nenabled    = false')
    s = s.replace('[output.monthly]\nenabled    = true', '[output.monthly]\nenabled    = false')
    fn = f'meds_config_prog_{month}.toml'
    open(fn, 'w').write(s)
    return fn

for m in MONTHS:
    print('wrote', gen(m))
