#!/bin/bash
# Pack the relativistic GW tables (rate_sve/gw_kerr/*.dat, built by build_gw_kerr_tables.jl) into the archive
# that is tracked in git, rate_sve/gw_kerr_tables.zip, and write its SHA256 next to it. The loose tables are
# not tracked; unpack with `cd src/rate_sve && unzip -o gw_kerr_tables.zip`.
#   bash scripts/pack_gw_kerr_tables.sh
set -e
cd "$(dirname "$0")/../rate_sve"
n=$(ls gw_kerr/*.dat | wc -l)
rm -f gw_kerr_tables.zip
ls gw_kerr/*.dat | LC_ALL=C sort | zip -q -X -9 gw_kerr_tables.zip -@
sha256sum gw_kerr_tables.zip | cut -d' ' -f1 > gw_kerr_tables.zip.sha256
echo "packed $n tables into $(pwd)/gw_kerr_tables.zip ($(du -h gw_kerr_tables.zip | cut -f1)); sha256 $(cat gw_kerr_tables.zip.sha256)"
