#! /bin/sh
#SBATCH -J WAR.count_s
#SBATCH --cpus-per-task 2
#SBATCH -p highmem
#SBATCH --mem 60G
#SBATCH --time 1-0

# Command line arguments
args=$*
site_code=${args[0]}

# Path variables
project="~/projects/WAR"
data="${project}/data"
results="${project}/results"
figures="${results}/figures"
assembly="${results}/assembly"
analysis="${results}/analysis"
logs="${results}/logs"
stats="${results}/stats"

jsondir="/srv/databases/internal/json"

assem_dir="${assembly}/${site_code}"
outdir="${analysis}/${site_code}"

gff_annot="${outdir}/${site_code}.contigs.annotated.mge.gff";
log="${logs}/coverage.log"

# Input check
if [[ ! -s $gff_annot ]]; then
  if [[ -s ${gff_annot}.gz ]]; then
    gff_annot="${gff_annot}.gz";
  else
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: no annotated features found in ${gff_annot} to count coverage of." >>${log};
    exit 1;
  fi
else
  gzip --best $gff_annot;
  gff_annot="${gff_annot}.gz";
fi

# Array of samples associated with site code
dna_samples="${data}/WAR.dna_samples.csv"
samples=()
while read line; do
  if [[ ${line:0:1} == "#" ]]; then
    continue  #ignore header
  fi
  site=$(echo $line | cut -d',' -f10);
  if [[ $site_code != $site ]]; then
    continue
  fi
  sample=$(echo $line | cut -d',' -f1);
  last_char="${sample: -1}";
  if [[ $(( $last_char % 2 )) -eq 0 ]]; then  #odd samples are sterivex
    continue;
  fi
  samples+=( $sample )
done < $dna_samples

# Associative array of total DNA extracted by sample
extractions="${data}/WAR.dna_extraction.csv"
declare -A dna_extracted
while read line; do
  if [[ ${line:0:1} == "#" ]]; then
    continue  #ignore comments
  fi
  sample_id=$(echo $line | cut -d',' -f1);
  last_char="${sample_id: -1}";
  if [[ $(( $last_char % 2 )) -eq 0 ]]; then  #odd samples are sterivex
    continue;
  fi
  for sample_name in ${samples[@]}; do
    if [[ $sample_name == $sample_id ]]; then
      edna=$(echo $line | cut -d',' -f19);
      fdna=$(echo $line | cut -d',' -f21);
      if [[ $edna != $fdna ]]; then  #sample diluted because concentration out of range of measuring device
	tdna=$(echo $fdna \* 47 / 5 | bc -l);  #total DNA of undiluted sample
      else
	tdna=$fdna;  #sample was undiluted
      fi
      dna_extracted[$sample_name]=$tdna;
    fi
  done
done < $extractions

# Calculate feature abundances and normalize for feature length and dataset size
for sample_id in ${samples[@]}; do

  mapfile="${assem_dir}/${sample_id}.mapped.sorted.derep.bam";
  map_stats="${stats}/${sample_id}.mapping_stats.derep.csv";
  if [[ ! -s $mapfile ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: warning: unable to obtain coverage information from ${mapfile}." >>${log};
    continue;
  else
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: sorting ${sample_id} BAM file by name." >>${log};
    map_temp="${assem_dir}/${sample_id}.mapped.sorted-name.derep.bam";
    samtools sort -l 0 -m 54G -n ${mapfile} >$map_temp;
    sleep 10s;
    map_sort_stats="${stats}/${sample_id}.mapping_stats.sort.csv";
    samtools stats ${map_temp} 1>${map_sort_stats};
    sleep 10s;
    total_mapped=$(grep "reads mapped:" ${map_stats} | cut -f 3);
    final_mapped=$(grep "reads mapped:" ${map_sort_stats} | cut -f 3);
    if [[ $total_mapped != $final_mapped ]]; then
      echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: error: mapped reads post-sorting ${final_mapped} is different from the expected ${total_mapped}." >>${log};
      exit 1;
    else
      mv $map_temp $mapfile;
    fi
  fi

  units=("prop" "fpk")
  outpref="${outdir}/${sample_id}.coverage.mge";
  if [[ ! -s ${outpref}.prop.csv && ! -s ${outpref}.fpk.csv ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: estimating feature abundances for sample ${sample_id}." >>${log};
    count_features --nonunique --gzip --qual 2 --mode union --attr Alias --type CDS --sorting name --format bam --units $(IFS=", "; shift; echo "${units[*]}") --outpref $outpref $mapfile $gff_annot 2>>${logs}/${site_code}.coverage.mge.log;
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: finished estimating feature abundances for sample ${sample_id}." >>${log};
  else
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: estimation of feature abundances for sample ${sample_id} has already been performed ... skipping." >>${log};
  fi

done
