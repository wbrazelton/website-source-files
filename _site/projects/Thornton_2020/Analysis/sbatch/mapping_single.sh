#! /bin/sh
#SBATCH -J war.mapping_single
#SBATCH --cpus-per-task 12
#SBATCH --mem 80000
#SBATCH --time 480
#SBATCH -p highmem

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
log="${logs}/mapping.log"
assem_dir="${assembly}/${site_code}";
renamed="${assem_dir}/${site_code}.contigs.renamed.fa";
index="${assem_dir}/${site_code}.contigs";

# Program paramters
cpus=12

# Data checks
if [[ ! -d $assem_dir || ! -s $renamed ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: assembly doesn't exist." >> ${log};
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
  fq="${qc}/${sample_id}.interleaved.atrim.decontam.qtrim.fq.gz";

  mapfile="${assem_dir}/${sample_id}.mapped.sorted.bam";
  if [[ ! -f $mapfile ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: starting read remapping." >> ${log};

    # Estimate insert-size distribution
    ishist="${stats}/${sample_id}.inserts.hist";
    if [[ ! -s $ishist ]]; then
      echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: calculating insert-size distribution." >> ${log};
      bbmap.sh -Xmx54g threads=${cpus} interleaved=true in=${fq} ihist=${ishist} reads=500000 ref=${renamed} nodisk;
      sleep 30s;
      echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: finished calculating insert-size distribution." >> ${log};
    fi

    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: estimating insert-size min and max from ${ishist}." >> ${log};
    ismean=$(grep "#Mean" ${ishist} | sed -n "/#Mean\t/s/#Mean\t//p" | xargs printf '%.*f\n' 0);
    issd=$(grep "#STDev" ${ishist} | sed -n "/#STDev\t/s/#STDev\t//p" | xargs printf '%.*f\n' 0);
    ismin=$(value=$(expr ${ismean} - ${issd} \* 3); if [[ ${value} -lt 0 ]]; then echo 0; else echo $value; fi);
    ismax=$(expr ${ismean} + ${issd} \* 3);
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: finished estimating insert-size min and max." >> ${log};

    # Map short reads to the assembly
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: mapping short reads to ${renamed}." >> ${log};
    bowtie2 -q --phred33 --very-sensitive -x ${index} -p $(( $cpus - 1 )) -I ${ismin} -X ${ismax} --interleaved ${fq} 2>${logs}/${sample_id}.mapping.log | samtools sort -@ 1 -m 56G -l 9 -O bam -T ${sample_id} -o ${mapfile} - 2>${logs}/${sample_id}.sort.log;
    sleep 30s;
    if [[ ! -s ${mapfile} ]]; then
      echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: failed to map short reads to the assembly. Please check the log for details." >> ${log};
      exit 1;
    fi
    samtools index ${mapfile};
    sleep 30s;
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: finished mapping short reads to ${renamed}." >> ${log};

    # Mapping statistics
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: generating mapping statistics." >> ${log};
    samtools stats ${mapfile} 1>${stats}/${sample_id}.mapping_stats.csv;
    sleep 30s;
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: verifying that mapping was successful." >> ${log};
    total_mapped=$(grep "raw total sequences:" ${stats}/${sample_id}.mapping_stats.csv | cut -f 3)
    total_qc=$(grep "Result:" ${logs}/${sample_id}.qtrim.log | cut -f 2 | awk '{print $1}')
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: ${total_mapped} reads successfully mapped out of ${total_qc}." >> ${log};
    if [[ $total_mapped != $total_qc ]]; then
      echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: error: read count post-mapping ${total_mapped} is less than the expected ${total_qc}." >> ${log};
      exit 1;
    else
      echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: finished generating mapping statistics." >> ${log};
    fi

    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: read remapping completed." >> ${log};
  else
    echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: short reads have already been mapped to the assembly ... skipping." >> ${log};
  fi
done
