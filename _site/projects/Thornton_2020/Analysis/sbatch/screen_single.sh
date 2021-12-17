#! /bin/sh
#SBATCH -J WAR.screen_s
#SBATCH --cpus-per-task 3
#SBATCH -p batch
#SBATCH --mem 16G
#SBATCH --time 180

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

jsondir="/srv/databases/internal/json/"
csvdir="/srv/databases/internal/csv/"

assem_dir="${assembly}/${site_code}"
outdir="${analysis}/${site_code}"

log="${logs}/screen.log"
contig="${assem_dir}/${site_code}.contigs.renamed.fa"

# Analysis parameters
cpus=1
ident_assem=80
cov_assem=25

# Create associative array of reference database names
declare -A dbs
dbs=(["amr_homologs"]="CARD homolog model" ["bmr"]="BacMet database" ["Pc"]="IntegronFinder cassette promotor database" ["vf"]="virulence factor database" ["ices"]="ICEberg database" ["transposons"]="ACLAME database")
variant_dbs=(["amr_overexpression"]="CARD overexpression model" ["amr_protein_variants"]="CARD variant model" ["amr_rRNA_variants"]="CARD rRNA variant model")

# Filter hits to the reference sequence databases
kegg="${outdir}/${site_code}.contigs.kegg.b6"
kegg_out="${outdir}/${site_code}.contigs.kegg.screened.b6"
if [[ ! -s $kegg_out && ! -s ${kegg_out}.gz ]]; then
  if [[ ! -s $kegg ]]; then
    if [[ ! -s ${kegg}.gz ]]; then
      echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: no hits found in ${kegg} to screen." >>${log};
      exit 1;
    else
      echo "$(date +"%b %-d %k:%M:%S") ${site_code}: screening hits to KEGG." >>${log};
      gzip -dc ${kegg}.gz | sort --parallel=1 --buffer-size=14G -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u --merge -k1,1 | gzip --best - >${kegg_out}.gz;
      echo "$(date +"%b %-d %k:%M:%S") ${site_code}: finished screening hits to KEGG." >>${log};
    fi
  else
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: screening hits to KEGG." >>${log};
    sort --parallel=1 --buffer-size=14G -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u --merge -k1,1 | gzip --best - >${kegg_out}.gz;
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: finished screening hits to KEGG." >>${log};
  fi
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: screening hits to KEGG has already been performed ... skipping." >>${log};
fi

for reference in "${!dbs[@]}"; do
  outfile="${outdir}/${site_code}.contigs.${reference}.screened.b6";
  if [[ ! -f $outfile ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: screening hits to the ${dbs[$reference]}." >>${log};
    infile="${outdir}/${site_code}.contigs.${reference}.b11";
    if [[ ! -s $infile ]]; then
      if [[ ! -s ${infile}.gz ]]; then
        echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: ${infile} cannot be found." >>${log};
        exit 1;
      else
        gzip -dc ${infile}.gz | blast_formatter -archive /dev/stdin -outfmt 6 -max_target_seqs 1 | screen_features --defaults --best --length ${cov_assem} --identity ${ident_assem} --out ${outfile} /dev/stdin 2>>${logs}/${site_code}.screen.log;
      fi
    else
      blast_formatter -archive ${infile} -outfmt 6 -max_target_seqs 1 | sort --parallel=1 --buffer-size=10G -k1,1 -k12,12gr -k11,11g -k3,3gr | screen_features --defaults --best --length ${cov_assem} --identity ${ident_assem} --out ${outfile} /dev/stdin 2>>${logs}/${site_code}.screen.log;
    fi
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: finished screening hits to the ${dbs[$reference]}." >>${log};
  else
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: screening hits to the ${dbs[$reference]} has already been performed ... skipping." >>${log};
  fi
done

#for reference in "${!variant_dbs[@]}"; do
#  outfile="${outdir}/${site_code}.contigs.${reference}.screened.b6";
#  if [[ ! -f $outfile ]]; then
#    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: screening hits to the ${variant_dbs[$reference]}." >>${log};
#    infile="${outdir}/${site_code}.contigs.${reference}.b11";
#    if [[ ! -s $infile ]]; then
#      if [[ ! -s ${infile}.gz ]]; then
#        echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: ${infile} cannot be found." >>${log};
#        exit 1;
#      else
#        gzip -dc ${infile}.gz | blast_formatter -archive /dev/stdin -max_target_seqs 1 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" | screen_features --specifiers "qaccver,saccver,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qseq,sseq" --mapping ${jsondir}/CARD.mapping.json --evalue ${evalue} --length ${cov_assem} --identity ${ident_assem} --alleles snp --defaults --out ${outfile} /dev/stdin 2>>${logs}/${site_code}.screen.log;
#      fi
#    else
#      blast_formatter -archive ${infile} -max_target_seqs 1 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" | screen_features --specifiers "qaccver,saccver,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qseq,sseq" --mapping ${jsondir}/CARD.mapping.json --evalue ${evalue} --length ${cov_assem} --identity ${ident_assem} --alleles snp --defaults --out ${outfile} /dev/stdin 2>>${logs}/${site_code}.screen.log;
#    fi
#    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: finished screening hits to the ${variant_dbs[$reference]}." >>${log};
#  else
#    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: screening hits to the ${variant_dbs[$reference]} has already been performed ... skipping." >>${log};
#  fi
#done
