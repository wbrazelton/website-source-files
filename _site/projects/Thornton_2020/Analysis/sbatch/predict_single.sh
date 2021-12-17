#! /bin/sh
#SBATCH -J WAR.predict_s
#SBATCH --cpus-per-task 1
#SBATCH --mem 6G
#SBATCH --time 180
#SBATCH -p batch

# Command line arguments
args=$*
site_code=${args[0]}

# Path variables
project="~/projects/WAR"
results="${project}/results"
assembly="${results}/assembly"
analysis="${results}/analysis"
logs="${results}/logs"
qc="${results}/qc"
stats="${results}/stats"

blastdir="/srv/databases/internal/blast/"
fastadir="/srv/databases/internal/fasta/"
jsondir="/srv/databases/internal/json/"
csvdir="/srv/databases/internal/csv/"
hmmdir="/srv/databases/internal/hmm/"
cmdir="/srv/databases/internal/cm/"

assem_dir="${assembly}/${site_code}"
outdir="${analysis}/${site_code}"

contig="${assem_dir}/${site_code}.contigs.renamed.fa"
log="${logs}/predict.log"

# Program parameters
cpus=8
evalue=0.00001

echo "$(date +"%b %-d %k:%M:%S") ${site_code}: starting functional profiling of the assembly." >>${log};

# Input check
if [[ ! -s $contig && ! -s ${contig}.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: no contigs found in ${contig} to profile." >>${log};
  exit 1
fi

# Output check
if [[ ! -d $outdir ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: creating output directory ${outdir}/${site_code}." >>${log};
  mkdir --parents ${outdir}/${site_code};
fi

# Index the assembly for searching
mkblast_log="${logs}/${site_code}.contigs.indexed.log";
if [[ ! -f ${mkblast_log} ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: indexing assembly for searching." >>${log};
  makeblastdb -dbtype nucl -max_file_sz 2GB -logfile ${mkblast_log} -title ${site_code} -in ${contig} -out ${assem_dir}/${site_code}.contigs;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: assembly indexing already completed ... skipping." >>${log};
fi

# Search for putative CRISPR arrays
crispr="${outdir}/${site_code}.contigs.crispr.gff";
if [[ ! -f ${crispr} && ! -f ${crispr}.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: searching for putative CRISPR arrays." >>${log};
  minced -minNR 3 -gffFull ${contig} ${crispr};
  sleep 30s;
  gzip --best ${crispr};
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: searching for CRISPRs already completed ... skipping." >>${log};
fi

# Use motif detection to find putative protein-coding DNA sequences
cds="${outdir}/${site_code}.contigs.cds.faa";
gff="${outdir}/${site_code}.contigs.cds.gff";
if [[ ! -f ${cds}.gz && ! -f ${cds} || ! -f ${gff}.gz && ! -f ${gff} ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: predicting protein-coding genes from the assembly." >>${log};
  if [[ -s ${contig}.gz ]]; then
      gzip -dc ${contig}.gz | prodigal -q -f gff -p meta -i /dev/stdin -o ${gff} -a ${cds};
  else
      prodigal -q -f gff -p meta -i ${contig} -o ${gff} -a ${cds};
  fi
  sleep 30s;
  if [[ ! -s $cds || ! -s $gff ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: no putative protein-coding genes could be identified in the assembly." >>${log};
    exit 1;
  else
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: removing invalid characters from the protein sequences and simplifying headers." >>${log};
    clean_cds="${outdir}/${site_code}.contigs.cds.clean.faa"
    awk '/^>/{split($0, a, " "); print a[1]; next} {gsub(/[-\*_\.]/,""); print}' ${cds} >${clean_cds}
    sed -i 's/\x0//g' ${clean_cds}  # NULL bytes: common undesirable artifact produced by prodigal for unknown reasons
    mv ${clean_cds} ${cds};
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: linking sequence IDs to GFF entries using the GFF's ID attribute." >>${log};
    clean_gff="${outdir}/${site_code}.contigs.cds.clean.gff"
    awk -F'\t' '/^#/{print $0; next} {sub(/ID\=[^_]*_/, "ID="$1"_", $9); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' ${gff} >${clean_gff};
    mv ${clean_gff} ${gff};
  fi
  gzip --best ${cds};
  gzip --best ${gff};
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: protein-coding gene prediction has already been performed ... skipping." >>${log};
fi
