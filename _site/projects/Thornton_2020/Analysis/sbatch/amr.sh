#! /bin/sh
#SBATCH -J WAR.final
#SBATCH --cpus-per-task 1
#SBATCH --mem 20GB

# Command line arguments
args=$*
SITE=${args[0]}

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

log="${logs}/final.log"
site_dir="${analysis}/${SITE}"
outdir="${analysis}/final/"
assem=$(input_check "${assembly}/${SITE}/${SITE}.contigs.renamed.fa")
ARG=$(input_check "${site_dir}/${SITE}.amrfinder.b6")
CARD=$(input_check "${site_dir}/${SITE}.contigs.card.screened.b6")
MGE=$(input_check "${site_dir}/${SITE}.contigs.mgefam.csv")
map_file=$(input_check "/home/cthornton/NCBIfam-AMR.mod.csv")

# Program parameters
cpus=1

mkdir --parents ${outdir}

# Create associative array mapping samples to sites
dna_samples="${data}/WAR.dna_samples.csv"
samples=()
while read line; do
  if [[ ${line:0:1} == "#" ]]; then
    continue  #ignore header
  fi
  site_code=$(echo $line | cut -d',' -f10);
  if [[ $SITE != $site_code ]]; then
    continue
  fi
  sample=$(echo $line | cut -d',' -f1);
  last_char="${sample: -1}";
  if [[ $(( $last_char % 2 )) -eq 0 ]]; then  #odd samples are sterivex
    continue;
  fi
  samples+=( $sample )
done < $dna_samples

# Find all AMR-bearing contigs and CDS annotated as an AMR
shared=($(awk '{print $1}' ${CARD} | grep -F -f - ${ARG} | awk '{print $1}'))
exception=($(awk '{print $1}' ${CARD} | grep -v -F -f - ${ARG} | awk '{print $2}' | grep -w -F -f - ${map_file} | grep "exception" | awk '{print $1}' | grep -F -f - ${ARG} | awk '{print $1}'))
allele=($(awk '$3 != "-" {print $1}' ${ARG}))
merged_arg=( "${shared[@]}" "${exception[@]}" "${allele[@]}" )

mge=($(awk '!/^#/ {print $0}' ${MGE} | sort --parallel=1 --buffer-size=4G -k1,1 -k6,6gr | sort -u --merge -k1,1 | grep -v -E "IstB_IS21_ATP" | awk '{print $1}'))

# Extract AMRs from FASTA file of assembled contigs
seed_arg="${outdir}/${SITE}.contigs.arg.fna"
if [[ ! -s $seed_arg ]]; then
  echo "${merged_arg[@]}" | tr ' ' '\n' | awk '{split($1, contig, "_"); print contig[1]"_"contig[2]}' | filterbyname.sh ow=t include=t names=/dev/stdin in=${assem} out=${seed_arg};
fi

# Extract AMRs from CDS abundance estimations
for sample_id in ${samples[@]}; do
  # Find abundance file for sample
  props="${site_dir}/${sample_id}.coverage.prop.csv";
  counts="${site_dir}/${sample_id}.coverage.counts.csv";
  if [[ ! -s $props ]]; then
    if [[ ! -s ${props}.gz ]]; then
      echo "$(date +"%b %-d %k:%M:%S") ${SITE}: unable to locate $props for applying feature abundances." >>${log};
      continue;
    else
      gunzip ${props}.gz;
      sleep 2s;
    fi
  fi

  # Extract ARGs and MBRGs from abundance tables
  arg_prop="${outdir}/${sample_id}.arg.prop.csv";
  mge_prop="${outdir}/${sample_id}.mge.prop.csv";
  if [[ ! -s $mge_prop || ! -s $arg_prop ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${SITE}: ${sample_id}: Extracting ARGs and MGEs from abundance table." >>${log};
    echo "${merged_arg[@]}" | tr ' ' '\n' | grep -w -F -f - ${props} >${arg_prop};
    echo "${mge[@]}" | tr ' ' '\n' | grep -w -F -f - ${props} >${mge_prop};
    sleep 2s;
  fi

done
