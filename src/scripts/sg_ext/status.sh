#!/bin/bash
# Progress of the extended-SG rate campaign and runs.
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext
cd $X
date
squeue -u $USER -h -o "%.12i %.12j %.2t %.11M %R" | grep sgx
for p in p0 p1 p2; do echo "chunks $p: $(ls done 2>/dev/null | grep -c "chunk_${p}_") done of $(ls tasks | grep -c "chunk_${p}_")"; done
echo "locks (chunks in progress): $(ls locks | grep -c chunk_)   failed-twice: $(for f in fail/*; do [ -e "$f" ] && [ "$(wc -l < $f)" -ge 2 ] && echo $f; done | wc -l)"
echo "channels with points: $(ls rates_out | wc -l) of $(wc -l < sgx_channels.txt)"
echo "FAILED points in chunk logs: $(cat logs/chunks/*.log 2>/dev/null | grep -c FAILED)"
[ -e CHAIN_STOPPED.txt ] && cat CHAIN_STOPPED.txt
ls runs 2>/dev/null | head
