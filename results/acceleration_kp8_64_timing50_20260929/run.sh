#!/usr/bin/env bash
set -u
cd /home/zm/CCM_CE/OpenSourceCodes_new || exit 1
/usr/local/bin/matlab -batch "exp_Acceleration(50,'results/acceleration_kp8_64_timing50_20260929',isfile('results/acceleration_kp8_64_timing50_20260929/acceleration.mat'),[8,16,32,64],26092900)" -logfile /home/zm/CCM_CE/OpenSourceCodes_new/results/acceleration_kp8_64_timing50_20260929/run.log
run_status=$?
printf '\nMATLAB exit status: %s\n' "$run_status"
exit "$run_status"
