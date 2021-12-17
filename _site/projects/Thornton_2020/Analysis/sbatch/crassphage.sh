#! /bin/sh
#SBATCH -J war.crAssphage
#SBATCH --cpus-per-task 2
#SBATCH -p batch
#SBATCH --mem 200
#SBATCH --time 120

# Command line arguments
args=$*
sample_id=${args[0]}

# Path variables
project="~/projects/WAR"
data="${project}/data"
results="${project}/results"
analysis="${results}/analysis"
logs="${results}/logs"
log="${logs}/crassphage.log"

forward="${data}/${sample_id}.01.forward.fastq.gz";
reverse="${data}/${sample_id}.01.reverse.fastq.gz";

outdir="${analysis}/crAssphage"
if [[ ! -d $outdir ]]; then
  mkdir --parents ${outdir};
fi

if [[ ! -s ${forward} || ! -s ${reverse} ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: missing one or more files of the read pair.">> ${log};
  exit 1;
fi

echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: mapping short reads to crAssphage genome." >>${log};
bowtie2 -x ${data}/crAssphage -1 ${forward} -2 ${reverse} --very-sensitive -q 2>${outdir}/${sample_id}.crassphage_mapping.log | samtools view -h -b -F 4 -o ${outdir}/${sample_id}.crassphage.bam -;
echo "$(date +"%b %-d %k:%M:%S") ${sample_id}: finished mapping short reads to crAssphage genome." >>${log};
