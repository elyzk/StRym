#!/usr/bin/env bash
# usage: verify.sh <mode> : map cuffdiff loci (XLOC) to truth genes by coordinates
d=/root/tuxtest/$1
awk -F'\t' '
FILENAME==ARGV[1]{ if(FNR>1){de[$1]=$3; fc[$1]=$4}; next }
FILENAME==ARGV[2]{ match($9,/gene_id "[^"]+"/); g=substr($9,RSTART+9,RLENGTH-10);
   if(!(g in lo)||$4<lo[g])lo[g]=$4; if($5>hi[g])hi[g]=$5; ch[g]=$1; next }
FNR>1{ split($4,a,/[:-]/); found="?";
   for(g in lo) if(ch[g]==a[1] && a[2]<=hi[g] && a[3]>=lo[g]){found=g}
   printf "%s\t%s\t%s\t%s\t%.2f\t%.3g\t%s\n",$2,found,de[found],fc[found],$10,$13,$14 }' \
  $d/truth.tsv $d/genes.gtf "$([ -d "$d/my data" ] && echo "$d/my data/out put" || echo "$d/out")/diff_out/gene_exp.diff" | sort -k2 \
  | awk -F'\t' 'BEGIN{print "XLOC\ttrue_gene\ttruth\ttrueFC\tlog2FC\tq\tsig"}{print; t=($3!="none"); s=($7=="yes"); tp+=t&&s; fn+=t&&!s; fp+=!t&&s; tn+=!t&&!s}
   END{printf "\nTP=%d FN=%d FP=%d TN=%d\n",tp,fn,fp,tn}'
