#!/usr/bin/env python3
"""Plot the fully-prognostic (leaf+wood) SPLIT/Picard 30s runs from the FAST netCDF
tier: air (forcing), tallest-cohort leaf, CAS, and surface-soil temperature, for Jan
and Jun. Tallest-cohort selection is done HERE (post-processing) via argmax(height).
Two panels per month: whole month + a representative 5-day diurnal zoom."""
import glob
import datetime as dt
import numpy as np
import netCDF4
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

K = 273.15
SERIES = [
    ('air',  'Air temperature (forcing)', '#111111', 1.4),
    ('leaf', 'Tallest-cohort leaf',       '#2ca02c', 1.1),
    ('cas',  'Canopy air space (CAS)',    '#ff7f0e', 1.1),
    ('soil', 'Surface soil',              '#8c564b', 1.3),
]

def load(month):
    """Concatenate the daily FAST files; pick the tallest cohort's leaf temp per record."""
    times, air, cas, soil, leaf = [], [], [], [], []
    for fn in sorted(glob.glob(f'out_nc/prog_{month}-F-*.nc')):
        d = netCDF4.Dataset(fn)
        yr, mo, dy = d['year'][:], d['month'][:], d['day'][:]
        hh, mi, ss = d['hour'][:], d['minute'][:], d['second'][:]
        ncoh = np.asarray(d['n_cohort'][:], dtype=int)
        ltemp = np.asarray(d['leaf_temp_cohort_fast'][:])   # (time, cohort)
        hgt   = np.asarray(d['height_cohort_fast'][:])      # (time, cohort)
        a = np.asarray(d['air_temp_fast'][:])
        c = np.asarray(d['cas_temp_site'][:])
        s = np.asarray(d['soil_temp_top_site'][:])
        for t in range(len(a)):
            n = max(1, ncoh[t])
            itall = int(np.argmax(hgt[t, :n]))              # tallest = max height among live cohorts
            times.append(dt.datetime(int(yr[t]), int(mo[t]), int(dy[t]), int(hh[t]), int(mi[t]), int(ss[t])))
            air.append(a[t]); cas.append(c[t]); soil.append(s[t]); leaf.append(ltemp[t, itall])
        d.close()
    order = np.argsort(times)
    times = [times[i] for i in order]
    out = dict(air=np.array(air)[order] - K, cas=np.array(cas)[order] - K,
               soil=np.array(soil)[order] - K, leaf=np.array(leaf)[order] - K)
    return times, out

def plot_month(month, month_label, zoom_start, png):
    t, data = load(month)
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(13, 8.5))
    for ax, title, window in ((ax0, f'{month_label} 2024 — whole month', None),
                              (ax1, f'{month_label} 2024 — {zoom_start:%b %-d}–{zoom_start + dt.timedelta(days=5):%-d} (diurnal detail)',
                               (zoom_start, zoom_start + dt.timedelta(days=5)))):
        for key, lab, col, lw in SERIES:
            ax.plot(t, data[key], color=col, lw=lw, label=lab)
        ax.axhline(0.0, color='#4488cc', lw=0.8, ls=':', alpha=0.7)   # freezing line
        ax.set_ylabel('Temperature (°C)'); ax.set_title(title, fontsize=11); ax.grid(alpha=0.25)
        if window:
            ax.set_xlim(window)
            ax.xaxis.set_major_locator(mdates.DayLocator())
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %-d'))
            ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6)))
        else:
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %-d'))
    ax0.legend(ncol=4, loc='upper center', bbox_to_anchor=(0.5, 1.28), frameon=False, fontsize=10)
    fig.suptitle('MEDS fully-prognostic leaf+wood energy (split/Picard, dt=30 s, FAST netCDF @15 min) — Ithaca NY',
                 y=0.995, fontsize=12, fontweight='bold')
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(png, dpi=130)
    rng = lambda k: (min(data[k]), max(data[k]))
    print('wrote', png, '|', ', '.join(f'{k}=[{rng(k)[0]:.1f},{rng(k)[1]:.1f}]' for k in ('air', 'leaf', 'cas', 'soil')))

plot_month('jan', 'January', dt.datetime(2024, 1, 15), 'out_nc/prog_jan.png')
plot_month('jun', 'June',    dt.datetime(2024, 6, 15), 'out_nc/prog_jun.png')
