#! /bin/sh
#SBATCH -J WAR.profile_s
#SBATCH --cpus-per-task 4
#SBATCH --mem 16G
#SBATCH --time 720
#SBATCH --partition batch

# Command line arguments
args=$*
site_code=${args[0]}

# Path variables
project="~/projects/WAR"
results="${project}/results"
assembly="${results}/assembly"
analysis="${results}/analysis"
logs="${results}/logs"

blastdir="/srv/databases/internal/blast/"
fastadir="/srv/databases/internal/fasta/"
jsondir="/srv/databases/internal/json/"
csvdir="/srv/databases/internal/csv/"
hmmdir="/srv/databases/internal/hmm/"
cmdir="/srv/databases/internal/cm/"

assem_dir="${assembly}/${site_code}"
outdir="${analysis}/${site_code}"

contig="${assem_dir}/${site_code}.contigs.renamed.fa"
cds="${outdir}/${site_code}.contigs.cds.faa"
log="${logs}/profile.log"

# Program parameters
cpus=4
evalue=0.00001

echo "$(date +"%b %-d %k:%M:%S") ${site_code}: starting functional profiling of the assembly." >>${log};

# Input check
if [[ ! -s $contig ]]; then
  if [[ ! -s ${contig}.gz ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: no contigs found to profile." >>${log};
    exit 1
  else
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: uncompressing input contigs." >>${log};
    gunzip ${contig}.gz;  #uncompress the file
  fi
fi

if [[ ! -s $cds ]]; then
  if [[ ! -s ${cds}.gz ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: error: no amino acid sequences in ${cds} to profile." >>${log};
    exit 1
  else
    echo "$(date +"%b %-d %k:%M:%S") ${site_code}: uncompressing input protein sequences." >>${log};
    gunzip ${cds}.gz;  #uncompress the file
  fi
fi

# Output check
if [[ ! -d $outdir ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: creating output directory ${outdir}/${site_code}." >>${log};
  mkdir --parents ${outdir}/${site_code};
fi

# Resistome analysis using assembled contigs
echo "$(date +"%b %-d %k:%M:%S") ${site_code}: identifying genes involved in resistance to antibiotics, metals, and biocides." >>${log};

## Homologs of protein-coding genes known to confer resistance to antibiotics
homologs="${outdir}/${site_code}.contigs.amr_homologs.b11"
if [[ ! -f $homologs && ! -f ${homologs}.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: searching for homologs of protein-coding genes known to confer resistance to antibiotics." >>${log};
  blastp -num_threads ${cpus} -max_hsps 1 -query $cds -db ${blastdir}/CARD-protein_homolog_model -outfmt 11 | gzip --best >${homologs}.gz;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: homology searching for antimicrobial resistance genes has already been performed ... skipping." >>${log};
fi

## Homologs of antibiotic target rRNA genes with resistance-conferring SNPs
rna_variants="${outdir}/${site_code}.contigs.amr_rRNA_variants.b11";
if [[ ! -f ${rna_variants} && ! -f ${rna_variants}.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: searching for homologs of rRNA genes known to confer resistance to antibiotics through mutation." >>${log};
  blastn -num_threads ${cpus} -max_hsps 1 -query ${contig} -db ${blastdir}/CARD-rRNA_gene_variant_model -outfmt 11 | gzip --best >${rna_variants}.gz;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: homology searching for antimicrobial resistance-conferring rRNA gene mutations has already been performed ... skipping." >>${log};
fi

## Homologs of protein-coding genes known to confer resistance to antibiotics through mutation
variants="${outdir}/${site_code}.contigs.amr_protein_variants.b11"
if [[ ! -f $variants && ! -f ${variants}.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: searching for homologs of protein-coding genes known to confer resistance to antibiotics through mutation." >>${log};
  blastp -num_threads ${cpus} -max_hsps 1 -query $cds -db ${blastdir}/CARD-protein_variant_model -outfmt 11 | gzip --best >${variants}.gz;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: homology searching for antimicrobial resistance-conferring protein-coding gene mutations has already been performed ... skipping." >>${log};
fi

## Homologs of protein-coding genes known to confer resistance to antibiotics through overexpression
overexpression="${outdir}/${site_code}.contigs.amr_overexpression.b11"
if [[ ! -f $overexpression && ! -f ${overexpression}.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: searching for homologs of protein-coding genes known to confer resistance to antibiotics through overexpression." >>${log};
  blastp -num_threads ${cpus} -max_hsps 1 -query $cds -db ${blastdir}/CARD-protein_overexpression_model -outfmt 11 | gzip --best >${overexpression}.gz;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: homology searching for overexpressed protein-coding genes has already been performed ... skipping." >>${log};
fi

## Homologs of genes known to confer resistance to biocides and metals:
bmr="${outdir}/${site_code}.contigs.bmr.b11"
if [[ ! -f $bmr && ! -f ${bmr}.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: searching for homologs of genes known to confer resistance to biocides and metals." >>${log};
  blastp -num_threads ${cpus} -max_hsps 1 -evalue ${evalue} -outfmt 11 -query $cds -db ${blastdir}/BacMet-EXP | gzip --best >${bmr}.gz;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: homology searching for metals and biocides resistance genes has already been performed ... skipping." >>${log};
fi

# Mobile genetic elements
echo "$(date +"%b %-d %k:%M:%S") ${site_code}: identifying genomic elements involved in the mobility of genes." >>${log};

## Transposases
tnp="${outdir}/${site_code}.contigs.tnp";
if [[ ! -f ${tnp}.csv && ! -f ${tnp}.csv.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: searching for homologs of known transposases." >>${log};
  hmmsearch --cpu $(( ${cpus} - 1)) --noali -E ${evalue} --tblout ${tnp}.csv ${hmmdir}/transposase.hmm ${cds} | gzip --best >${tnp}.log.gz;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: transposase search already completed ... skipping." >>${log};
fi

## Integrons
echo "$(date +"%b %-d %k:%M:%S") ${site_code}: identifying potential integrons by searching for primary components." >>${log};
### Integron-integrases
inti="${outdir}/${site_code}.contigs.intI";
if [[ ! -f ${inti}.csv && ! -f ${inti}.csv.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: searching for homologs of known integron integrases." >>${log};
  hmmsearch --cpu $(( ${cpus} - 1)) --noali -E ${evalue} --tblout ${inti}.csv ${hmmdir}/integron_integrase.hmm ${cds} | gzip --best >${inti}.log.gz;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: integron integrase search already completed ... skipping." >>${log};
fi

### Gene cassette promotors
pc="${outdir}/${site_code}.contigs.Pc.b11";
if [[ ! -f $pc && ! -f ${pc}.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: searching for potential cassette promoters." >>${log};
  blastn -num_threads ${cpus} -max_hsps 1 -evalue ${evalue} -query ${fastadir}/variants_Pc_intI1.fst -db ${assem_dir}/${site_code}.contigs -outfmt 11 | gzip --best >${pc}.gz;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: search for cassette promoters already completed ... skipping." >>${log};
fi

### ATTC recombination sites
attc="${outdir}/${site_code}.contigs.attC";
if [[ ! -f ${attc}.csv && ! -f ${attc}.csv.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: scanning contigs for attC recombination sites." >>${log};
  cmsearch --cpu $(( ${cpus} - 1)) --noali -E ${evalue} --tblout ${attc}.csv ${cmdir}/attc_4.cm ${contig} | gzip --best >${attc}.log.gz;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: scanning for attC recombination sites already completed ... skipping." >>${log};
fi
echo "$(date +"%b %-d %k:%M:%S") ${site_code}: finished identifying potential integrons in the assembly." >>${log};

## Transposons, plasmids, ICEs, and phages

### Transposons
#transposons="${outdir}/${site_code}.contigs.transposons.b11";
#if [[ ! -f $transposons && ! -f ${transposons}.gz ]]; then
#  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: identifying transposons and other mobile genetic elements." >>${log};
#  gzip -dc ${cds}.gz | blastp -num_threads ${cpus} -max_hsps 1 -evalue ${evalue} -outfmt 11 -out /dev/stdout -query /dev/stdin -db ${blastdir}/ACLAME | gzip --best >${transposons}.gz;
#else
#  echo "$(date +"%b %-d %k:%M:%S") ${site_code}:  identification fo transposons already completed ... skipping." >>${log};
#fi

# Virulence factors
vf="${outdir}/${site_code}.contigs.vf.b11";
if [[ ! -f $vf && ! -f ${vf}.gz ]]; then
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: identifying genes involved in pathogen virulence." >>${log};
  blastp -num_threads ${cpus} -max_hsps 1 -evalue ${evalue} -outfmt 11 -out /dev/stdout -query $cds -db "${blastdir}/VFDB-setA ${blast_dir}/VFDB-setB" | gzip --best >${vf}.gz;
else
  echo "$(date +"%b %-d %k:%M:%S") ${site_code}: identification of genes involved in pathogen virulence already completed ... skipping." >>${log};
fi

echo "$(date +"%b %-d %k:%M:%S") ${site_code}: finished functional profiling of the assembly." >>${log};

# Final Steps
echo "$(date +"%b %-d %k:%M:%S") ${site_code}: recompressing input files." >>${log};
gzip --best $contig;
gzip --best $cds;
