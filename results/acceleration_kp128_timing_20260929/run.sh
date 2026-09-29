#!/usr/bin/env bash
set -u
cd /home/zm/CCM_CE/OpenSourceCodes_new || exit 1
# User stopped the experiment after 75 completed paired trials.
/usr/local/bin/matlab -batch "exp_Acceleration(75,'results/acceleration_kp128_timing_20260929',isfile('results/acceleration_kp128_timing_20260929/acceleration.mat'),128,26092900)" -logfile /home/zm/CCM_CE/OpenSourceCodes_new/results/acceleration_kp128_timing_20260929/run.log
run_status=$?
printf '\nMATLAB exit status: %s\n' "$run_status"
exit "$run_status"
