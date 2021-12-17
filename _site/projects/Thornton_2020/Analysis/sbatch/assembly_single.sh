#! /bin/sh
#SBATCH -J assembly_single
#SBATCH --cpus-per-task 16
#SBATCH -p highmem
#SBATCH --mem 107520
#SBATCH --time 1200

# Command line arguments
args=$*
site_code=${args[0]}

# Path variables
project="~/projects/WAR"
data="${project}/data"
results="${project}/results"
assembly="${results}/assembly"
analysis="${results}/analysis"
logs="${results}/logs"
qc="${results}/qc"
stats="${results}/stats"
log="${logs}/assembly.log"

# Program parameters
cpus=16
kstart=27
kend=127
kstep=20

# Create array of samples by associated site
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

# Combine samples by site
singles=();
pairs=();
for sample_id in ${samples[@]}; do
  pair="${qc}/${sample_id}.interleaved.atrim.decontam.qtrim.fq.gz";
  single="${qc}/${sample_id}.singles.atrim.decontam.qtrim.fq.gz";
  if [[ ! -s $pair ]]; then
    continue;  #don't include failed samples
  fi
  pairs+=( $pair );
  singles+=( $single );
done
pair_combined=$(joined=$(printf ",%s" "${pairs[@]}"); echo ${joined:1});
single_combined=$(joined=$(printf ",%s" "${singles[@]}"); echo ${joined:1});

# Combined assembly
assem_dir="${assembly}/${site_code}";
orig="${assem_dir}/${site_code}.contigs.fa";
renamed="${assem_dir}/${site_code}.contigs.renamed.fa";
if [[ ! -d $assem_dir ]]; then
  if [[ ! -z $pair_combined && ! -z $single_combined ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: starting assembly of short reads with Megahit." >> ${log};
    megahit -t ${cpus} --memory 0.40 --12 ${pair_combined} --read ${single_combined} --out-prefix ${site_code} --out-dir ${assem_dir} --min-contig-len 200 --k-min ${kstart} --k-max ${kend} --k-step ${kstep};
    sleep 30s;
    if [[ ! -f "${assem_dir}/done" ]]; then
      echo "$(date +"%b %-d %k:%M:%S") ${site_code}: assembly failed. See assembly log for details" >> ${log};
      exit 1;
    fi 
    if [[ ! -s $orig ]]; then
      echo "$(date +"%b %-d %k:%M:%S") ${site_code}: assembly failed. No contigs were created." >> ${log};
      exit 1;
    fi
  else
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: assembly failed. No short reads available to assemble." >> ${log};
    exit 1;
  fi
else
  if [[ ! -f $orig ]]; then
    if [[ ! -f $renamed ]]; then
      echo "$(date +"%b %-d %k:%M:%S") ${site_code}: restarting an incomplete assembly run." >> ${log};
      megahit --continue --out-dir ${assem_dir};
      sleep 30s;
    fi
  fi
fi

# Rename headers
int_contig="${assem_dir}/intermediate_contigs/";
index="${assem_dir}/${site_code}.contigs";
if [[ ! -s $renamed ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: simplifying contig names for downstream use." >> ${log};
  anvi-script-reformat-fasta ${orig} -o ${renamed} -l 0 --simplify-names --report-file ${assem_dir}/${site_code}.header_map.tsv;
  sleep 30s;
  if [[ ! -s $renamed ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: assembly failed. ${renamed} is empty." >> ${log};
    exit 1;
  else
    # Assembly statistics
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: generating basic assembly statistics." >> ${log};
    stats.sh in=${renamed} format=2 >${logs}/${site_code}.assembly_stats.txt;
    sleep 30s;

    # Discard intermediates
    n_orig=$(stats.sh in=${orig} format=2 | grep "Main genome contig total:" | cut -f 2);
    n_renamed=$(grep "Main genome contig total:" ${logs}/${site_code}.assembly_stats.txt | cut -f 2);
    if [[ $n_orig == $n_renamed ]]; then
      echo "$(date +"%b %-d %k:%M:%S") ${site_code}: removing original assembly." >> ${log};
      rm ${orig};
      echo "$(date +"%b %-d %k:%M:%S") ${site_code}: removing intermediate files in ${int_contig}." >> ${log};
      rm -R ${int_contig};
    else
      echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: original (${n_orig}) and renamed (${n_renamed}) assemblies do not have the same number of contigs." >> ${log};
      exit 1;
    fi

    # Index contigs for short read mapping
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: generating bowtie2 index for short read mapping." >> ${log};
    bowtie2-build ${renamed} ${index} 1>${logs}/${site_code}.bt_index.log 2>${logs}/${site_code}.bt_index.err;
  fi
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: finished assembly of short reads." >> ${log};
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: assembly of short reads already completed ... skipping." >> ${log};
  continue
fi
