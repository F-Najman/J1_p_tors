#!/bin/sh
# Runs all the computations of Section 6 and writes the logs to logs/.
# Each log ends with "Claim n successfully verified" lines; an assertion failure stops the script.
set -e
cd "$(dirname "$0")"
mkdir -p logs
magma -b congruence_module.m > logs/congruence_module.log 2>&1
magma -b orders.m            > logs/orders.log            2>&1
magma -b two_rank.m          > logs/two_rank.log          2>&1
magma -b level15.m           > logs/level15.log           2>&1
grep -h "successfully verified" logs/*.log
