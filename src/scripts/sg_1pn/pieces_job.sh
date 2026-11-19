#!/bin/bash
#SBATCH --partition=maxcpu --constraint="7402|75F3|7H12" --nodes=1 --time=0-06:00:00
#SBATCH --job-name=sg_1pn_pieces
# All pieces of the complete 1PN correction of the self-gravity amplitude (see README.md), for the points of a job
# file with lines "S1 S2 S3 alpha" (the format of ongoing_tests/sg_relativistic/scan_gen.jl).
#   sbatch --export=ALL,WORK=<work dir>,JOBS=<job file>,SRC_TESTS=<repository>/src/ongoing_tests pieces_job.sh
#   then:  python3 complete_1pn.py <work dir>          (NW = workers per set, default 16; JULIA = julia binary)
# 1. Newtonian amplitudes with Kerr bound states (scan_gen.jl newton_only), five times: as in the rate tables (n0), in
#    the harmonic radius (shift), and with the weights of the 1PN source (w40: 4 V_B where the field acts; w46: also
#    6 V_B at the source; wrec: also the recoil weight 5/2 V_B and retardation in the acceleration of the hole).
# 2. Gravitomagnetism and retardation (sg_1pn_gen.jl) and the induction term of the recoil (sg_1pn_recoil.jl), with
#    hydrogenic states.
source ~/.bashrc
: ${WORK:?set WORK}; : ${JOBS:?set JOBS}
S=${SRC_TESTS:?set SRC_TESTS to <repository>/src/ongoing_tests}
H=$S/../scripts/sg_1pn; JULIA=${JULIA:-$HOME/julia-1.8.0/bin/julia}; NW=${NW:-16}
export OPENBLAS_NUM_THREADS=1
mkdir -p $WORK; cp $JOBS $WORK/jobs.txt; cd $WORK
run() { # <dir> <environment assignments ...>
  d=$1; shift; mkdir -p $d/logs
  for w in $(seq 1 $NW); do ( env "$@" JULIA_NUM_THREADS=2 $JULIA $S/sg_relativistic/scan_gen.jl jobs.txt $d newton_only > $d/logs/w$w.log 2>&1 ) & sleep 1; done
}
run n0    SGX_NEWTON_SHIFT=0
run shift SGX_NEWTON_SHIFT=1
run w40   SGX_NEWTON_SHIFT=1 SGX_NEWTON_WX=4
run w46   SGX_NEWTON_SHIFT=1 SGX_NEWTON_WX=4 SGX_NEWTON_WY=6
run wrec  SGX_NEWTON_SHIFT=1 SGX_NEWTON_WX=4 SGX_NEWTON_WY=6 SGX_NEWTON_WRX=2.5 SGX_NEWTON_RECRET=1
wait
mkdir -p pn ind
python3 $H/complete_1pn.py $WORK --tasks > tasks_1pn.txt
cd $S/sg_vs_si
grep '^pn '  $WORK/tasks_1pn.txt | cut -d' ' -f2- | xargs -P 40 -L 1 bash -c "$JULIA sg_1pn_gen.jl \$0 \$1 \$2 \$3 \$4 \$5 $WORK/pn 2>&1 | tail -n 2" > $WORK/pn.log 2>&1 &
grep '^ind ' $WORK/tasks_1pn.txt | cut -d' ' -f2- | xargs -P 20 -L 1 bash -c "$JULIA sg_1pn_recoil.jl \$0 \$1 \$2 \$3 $WORK/ind 2>&1 | tail -n 2" > $WORK/ind.log 2>&1 &
wait; echo "end $(date)"
