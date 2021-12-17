#! /bin/sh
#SBATCH -J WAR.derep
#SBATCH --cpus-per-task 2
#SBATCH --mem 80G
#SBATCH --time 360
#SBATCH --partition highmem

# Command line arguments
args=$*
site_code=${args[0]}

# Path variables
project="~/projects/WAR"
data="${project}/data"
results="${project}/results"
assembly="${results}/assembly"
analysis="${results}/analysis"
stats="${results}/stats"
logs="${results}/logs"
qc="${results}/qc"
log="${logs}/derep.log"
assem_dir="${assembly}/${site_code}";
renamed="${assem_dir}/${site_code}.contigs.renamed.fa";
index="${assem_dir}/${site_code}.contigs";

# Data checks
if [[ ! -d $assem_dir || ! -s $renamed ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: assembly doesn't exist." >>${log};
  exit 1
fi

# Create associative array mapping samples to sites
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

for sample_id in ${samples[@]}; do
  mapfile="${assem_dir}/${sample_id}.mapped.sorted.bam";
  if [[ ! -f $mapfile ]]; then
    continue
  fi

  map_derep="${assem_dir}/${sample_id}.mapped.sorted.derep.bam";
  metrics="${stats}/${sample_id}.replicates.txt";
  if [[ ! -f $map_derep ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: removing replicate mapped reads." >>${log};
    picard -Xmx65g MarkDuplicates COMPRESSION_LEVEL=9 USE_JDK_INFLATER=true USE_JDK_DEFLATER=true REMOVE_DUPLICATES=false INPUT=${mapfile} OUTPUT=${map_derep} METRICS_FILE=${metrics} 2>${logs}/${sample_id}.derep.log;
    sleep 30s;
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: finished dereplicating mapped reads." >>${log};
  fi

  if [[ ! -s ${map_derep} ]]; then
      echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: failed to dereplicate mapped reads." >>${log};
      exit 1;
  fi

  map_stats="${stats}/${sample_id}.mapping_stats.csv";
  map_derep_stats="${stats}/${sample_id}.mapping_stats.derep.csv";
  if [[ ! -s ${map_derep_stats} ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: generating mapping statistics on dereplicated mapped reads." >>${log};
    samtools stats ${map_derep} 1>${map_derep_stats};
    sleep 30s;
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: finished generating mapping statistics on dereplicated mapped reads." >>${log};
  fi

  echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: verifying that dereplication of mapped reads was successful." >>${log};
  if [[ ! -s ${map_stats} || ! -s ${map_derep_stats} ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: error: cannot read mapping statistics file." >>${log};
    exit 1;
  fi

  total_mapped=$(grep "reads mapped:" ${map_stats} | cut -f 3)
  final_mapped=$(grep "reads mapped:" ${map_derep_stats} | cut -f 3);
  if [[ $total_mapped != $final_mapped ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: error: mapped reads post-dereplication ${final_mapped} is different from the expected ${total_mapped}." >>${log};
    exit 1;
  else
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: finished generating mapping statistics." >>${log};
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: removing old mapping file." >>${log};
    rm $mapfile;
    rm ${mapfile}.bai;
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: finished removing old mapping file." >>${log};
  fi

  echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: indexing dereplicated mapped reads." >>${log};
  samtools index ${map_derep};
  echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: finished indexing of dereplicated mapped reads." >>${log};

done
