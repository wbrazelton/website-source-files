#! /bin/sh
#SBATCH -J WAR.annotate_s
#SBATCH --cpus-per-task 2
#SBATCH --mem 30G
#SBATCH --time 20
#SBATCH --partition batch

# Command line arguments
args=$*;
site_code=${args[0]};

# Path variables
project="~/projects/WAR";
results="${project}/results";
assembly="${results}/assembly";
analysis="${results}/analysis";
logs="${results}/logs";

blastdir="/srv/databases/internal/blast/";
fastadir="/srv/databases/internal/fasta/";
jsondir="/srv/databases/internal/json/";
csvdir="/srv/databases/internal/csv/";
hmmdir="/srv/databases/internal/hmm/";
cmdir="/srv/databases/internal/cm/";

assem_dir="${assembly}/${site_code}";
outdir="${analysis}/${site_code}";

log="${logs}/annotate.log";
gff="${outdir}/${site_code}.contigs.cds.gff";
gff_crispr="${outdir}/${site_code}.contigs.crispr.gff";

#screened=("${outdir}/${site_code}.contigs.amr_homologs.screened.b6" "${outdir}/${site_code}.contigs.amr_overexpression.screened.b6" "${outdir}/${site_code}.contigs.amr_protein_variants.screened.b6" "${outdir}/${site_code}.contigs.amr_rRNA_variants.screened.b6" "${outdir}/${site_code}.contigs.bmr.screened.b6" "${outdir}/${site_code}.contigs.kegg.screened.b6");
screened_blast=("${outdir}/${site_code}.contigs.amr_homologs.screened.b6" "${outdir}/${site_code}.contigs.bmr.screened.b6" "${outdir}/${site_code}.contigs.kegg.screened.b6");
screened_hmmer=("${outdir}/${site_code}.contigs.intI.b6" "${outdir}/${site_code}.contigs.tnp.screened.b6");
mappings=("${jsondir}BacMet-EXP.mapping.json" "${jsondir}CARD.mapping.json" "${jsondir}KEGG-prokaryotes.mapping.json");

# Input check
for i in "${!screened_blast[@]}"; do
  input=${screened_blast[$i]};
  if [[ ! -s $input ]]; then
    if [[ ! -s ${input}.gz ]]; then
      echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: cannot find ${input}." >>${log};
      exit 1;
    else
      screened_blast[$i]="${input}.gz";
    fi
  fi
done

for input in "${mappings[@]}"; do
  if [[ ! -s $input ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: cannot find ${input}." >>${log};
    exit 1;
  fi
done

if [[ ! -s $gff ]]; then
  if [[ ! -s ${gff}.gz ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: no amino acid sequences in ${gff} to annotate." >>${log};
    exit 1;
  else
    gff="${gff}.gz";
  fi
else
  gzip --best $gff;
  gff="${gff}.gz";
fi

if [[ ! -f $gff_crispr ]]; then
  if [[ -f ${gff_crispr}.gz ]]; then
    gff_crispr="${outdir}/${site_code}.contigs.crispr.gff.gz";
  else
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: no CRISPR sequences in ${gff_crispr} to annotate." >>${log};
    exit 1;
  fi
fi

# Output check
if [[ ! -d $outdir ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: creating output directory ${outdir}." >>${log};
  mkdir --parents ${outdir};
fi

# Annotate metagenome using BLAST results
gff_res="${outdir}/${site_code}.contigs.annotated.gff";
if [[ ! -s $gff_res && ! -s ${gff_res}.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: merging results of BLAST homology searches and annotating features." >>${log};
  gzip -dc $gff | annotate_features --clear --type CDS --product "hypothetical protein" --mapping $(IFS=", "; shift; echo "${mappings[*]}") --fields gene,Dbxref,Ontology_term,product,organism,database --conflict quality --out ${gff_res}.gz --discarded ${logs}/${site_code}.discarded.tsv --b6 $(IFS=", "; shift; echo "${screened_blast[*]}") /dev/stdin 2>${logs}/${site_code}.annotate.log;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: merging and annotation of BLAST homology searches has already been completed ... skipping." >>${log};
fi

# Annotate metagenome using HMMER results
gff_mge="${outdir}/${site_code}.contigs.annotated.mge.gff";
if [[ ! -s $gff_mge && ! -s ${gff_mge}.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: merging results of HMM homology searches and annotating features." >>${log};
  gzip -dc $gff | annotate_features --clear --filter --type CDS --conflict quality --out ${gff_mge}.gz --discarded ${logs}/${site_code}.discarded.tsv --b6 $(IFS=", "; shift; echo "${screened_hmmer[*]}") /dev/stdin 2>>${logs}/${site_code}.annotate.log;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: merging and annotation of HMM homology searches has already been completed ... skipping." >>${log};
fi

# Merge predicted features
gff_merged="${outdir}/${site_code}.contigs.annotated.merged.gff";
if [[ ! -s $gff_merged && ! -s ${gff_merged}.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: merging feature annotations." >>${log};
  combine_features --preserve --precedence --out ${gff_merged}.gz ${gff_crispr} ${gff_mge}.gz ${gff_res}.gz 2>${logs}/${site_code}.combine.log;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: merging of feature annotations has already been performed ... skipping." >>${log};
fi
