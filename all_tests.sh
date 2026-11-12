#!/usr/bin/env bash
cd "$(dirname "$0")"
for t in "single win flat gz" "paired win flat gz" "single local folders"; do
  bash run_test.sh $t; echo "----- verify:"; bash verify.sh $(echo $t | awk '{n=$1"_"$2"_"$3; if($4)n=n"_"$4; print n}') | tail -3; echo "=========="
done
