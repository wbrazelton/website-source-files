#! /bin/sh
#SBATCH -J WAR.cooccur
#SBATCH --cpus-per-task 1
#SBATCH --mem 2GB
#SBATCH --time 10:0
#SBATCH --error cooccur.err

# Functions
input_check() {
  if [[ ! -s $1 ]]; then
    if [[ ! -s ${1}.gz ]]; then
        echo "$(date +"%b %-d %k:%M:%S") ${SITE}: error: cannot find $1." >>${log};
        exit 1;
    else
      local input="${1}.gz"
    fi
  else
    local input=$1
  fi
  echo $input
}

# Command line arguments
args=$*
SITE=${args[0]}

# Path variables
project="WAR"
project_dir="~/projects/${project}"
results="${project_dir}/results"
data="${project_dir}/data"
logs="${results}/logs"
analysis="${results}/analysis"
assembly="${results}/assembly"
qc="${results}/qc"
stats="${results}/stats"
log="${logs}/cooccur.log"

outdir="${analysis}/${SITE}"
contigs=$(input_check "${assembly}/${SITE}/${SITE}.contigs.renamed.fa")
mge_gff=$(input_check "${outdir}/${SITE}.contigs.annotated.mgefam.gff")
amr_gff=$(input_check "${outdir}/${SITE}.contigs.annotated.amrfinder.gff")
sets=$(input_check "${data}/sets.csv")
#set1=$(input_check "${data}/set1.csv")
#set2=$(input_check "${data}/set2.csv")

out_gff="${outdir}/${SITE}.contigs.annotated.arg-mge.gff"
occur_table="${outdir}/${SITE}.cooccur.table.tsv"
occur_dist="${outdir}/${SITE}.cooccur.dist.tsv"

# Output check
if [[ ! -d $outdir ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${SITE}: creating output directory ${outdir}." >>${log};
  mkdir --parents ${outdir};
fi

# Combine feature annotations
if [[ ! -s $out_gff ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${SITE}: combining feature annotations." >>${log};
  combine_features --precedence --preserve --stranded --out ${out_gff} ${mge_gff} ${amr_gff};
else
  echo "$(date +"%b %-d %k:%M:%S") ${SITE}: feature annotations already combined." >>${log};
fi

# Find co-located ARGs and TEs
if [[ -s $out_gff ]]; then
  if [[ ! -s $occur_dist || ! -s $occur_table ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${SITE}: finding co-located ARGs and TEs." >>${log};
#    colocate_features --attr gene --out-dist ${occur_dist} --out ${occur_table} --fasta ${contigs} --set ${set1} --set2 ${set2} ${out_gff};
    colocate_features --symmetric --attr gene --out-dist ${occur_dist} --out ${occur_table} --fasta ${contigs} --set ${sets} ${out_gff};

  else
    echo "$(date +"%b %-d %k:%M:%S") ${SITE}: searching for co-located ARGs and TEs already performed." >>${log};
  fi
fi
