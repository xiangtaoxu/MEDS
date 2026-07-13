#!/bin/bash
cd /home/xiangtao/projects/MEDS/runs/ithaca_ark30
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1
export LD_LIBRARY_PATH=$HOME/miniforge3/envs/common/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
BIN=../../build-ifx/meds_main
echo "prefix,wall_sec,exit" > integ/timings.csv
while read p; do
  [ -z "$p" ] && continue
  rm -f out_integ/${p}-F-*.nc out_integ/${p}-D-*.nc
  t0=$(date +%s.%N)
  $BIN integ/${p}.toml >integ/log_${p}.txt 2>/dev/null
  rc=$?; t1=$(date +%s.%N)
  w=$(echo "$t1-$t0"|bc)
  echo "${p},${w},${rc}" >> integ/timings.csv
  echo "  done ${p}: wall=${w}s exit=${rc} files=$(ls out_integ/${p}-F-*.nc 2>/dev/null|wc -l)"
done < integ/run_order.txt
echo "MATRIX_COMPLETE"
