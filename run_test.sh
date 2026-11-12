#!/usr/bin/env bash
# usage: run_test.sh single|paired [local|win] [folders|flat] [gz]
#   Synthetic end-to-end test. Data in /root/tuxtest/<name>; outputs in <name>/out.
#   win : copy FASTQ inputs to /mnt/c/tuxtest_inputs_<name> (a Windows-drive mount,
#         like an external hard drive) and read them from there.
#   flat: all FASTQ in one folder, groups inferred from names; gz: gzip the FASTQs.
set -e
mode=${1:-single}; where=${2:-local}; layout=${3:-folders}; gz=${4:-}
name=${mode}_${where}_${layout}${gz:+_gz}
export PATH=/root/mamba/envs/tuxnet/bin:$PATH
src=$(cd "$(dirname "$0")" && pwd)
d=/root/tuxtest/$name; rm -rf $d; mkdir -p $d
flags="--layout $layout"; [ "$mode" = paired ] && flags="$flags --paired"; [ -n "$gz" ] && flags="$flags --gz"
Rscript $src/make_synthetic.R --out $d $flags --seed 7
inputs=$d/inputs
if [ "$where" = win ]; then
  inputs=/mnt/c/tuxtest_inputs_$name; rm -rf "$inputs"; cp -r $d/inputs "$inputs"
fi
gtf=$d/genes.gtf; fasta=$d/genome.fa; outdir=$d/out
if [ "$where" = spaces ]; then      # spaces in every path (inputs, reference, output)
  sp="$d/my data"; mkdir -p "$sp/in puts" "$sp/ref files"; cp -r $d/inputs/. "$sp/in puts/"
  cp $d/genes.gtf "$sp/ref files/my genes.gtf"; cp $d/genome.fa "$sp/ref files/my genome.fa"
  inputs="$sp/in puts"; gtf="$sp/ref files/my genes.gtf"; fasta="$sp/ref files/my genome.fa"; outdir="$sp/out put"
fi
cat > $d/cfg.json <<J
{"inputs_folder":"$inputs","output_folder":"$outdir","gtf_file":"$gtf","fasta_file":"$fasta",
 "read_type":"$mode","threads":4,"adapter_fasta":"$src/IlluminaAdaptorSeq.fasta"}
J
set +e
Rscript $src/run_cli.R $d/cfg.json > $d/run.log 2>&1; rc=$?
grep -E "samples in|Input size|WARNING|ERROR|Done|failed|STEP" $d/run.log
echo "[$name] EXIT=$rc"
echo "inputs dir listing:"; ls "$inputs" | head -8
echo "samples_used.csv:"; cat "$outdir"/samples_used.csv 2>/dev/null | head -8
