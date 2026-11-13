#!/usr/bin/env bash
# Tests with spaces in every path (inputs, reference, output): single and paired, flat gz.
cd "$(dirname "$0")"
for mode in single paired; do
  bash run_test.sh $mode spaces flat gz
  echo "--- verify"; bash verify.sh ${mode}_spaces_flat_gz | tail -2
  echo "--- leftover links:"; ls /tmp/Rtmp*/tuxnet_links 2>/dev/null | head -2
  echo "======"
done
