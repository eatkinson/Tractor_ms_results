###shortened technique.

##Pre-step: run RFmixv2, or v1 with converted output

# Remove all INFO fields and all FORMAT fields except for GT and PL from the input VCF, if necessary
use Bcftools
for i in {1..22}; do \
qsub -l h_vmem=50G wrapper1.sh /broad/software/free/Linux/redhat_6_x86_64/pkgs/bcftools/bcftools_1.1/bin/bcftools annotate -x INFO,^FORMAT/GT AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr$i.recode.vcf > AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr$i.noINFO.vcf ; done

##Run Tractor
reuse  Python-2.7
for i in {1..21}; do \
python /humgen/atgu1/fs03/eatkinso/scripts/Tractor/UnkinkMSPfile-Flags-switchpoints.py --msp UKBB_AfEur_QCed_lipids.rfmix.chr$i ; done

for i in {1..22}; do \
qsub wrapper1.sh /humgen/atgu1/fs03/eatkinso/scripts/Tractor/UnkinkGenofile-Flags.py --switches  UKBB_AfEur_QCed_lipids.rfmix.chr$i.switches.txt --genofile AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr$i.noINFO \
; done

for i in {1..22}; do \
qsub -l h_vmem=100G -N dosage wrapper1.sh python /humgen/atgu1/fs03/eatkinso/scripts/ExtractTracts-Flags-dosage-minor1.py --msp UKBB_AfEur_QCed_lipids.rfmix.chr$i --vcf AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr$i.noINFO1


#export dosage files if want to run jointly on Google cloud
cat AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr*.noINFO1.anc1.dosage.txt > temp
cat DosHeader temp > UKBB_AfEur_QCed_lipids.autosomes.anc1.dosage.txt
bgzip UKBB_AfEur_QCed_lipids.autosomes.anc1.dosage.txt

cat AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr*.noINFO1.anc0.dosage.txt > temp
cat DosHeader temp > UKBB_AfEur_QCed_lipids.autosomes.anc0.dosage.txt
bgzip UKBB_AfEur_QCed_lipids.autosomes.anc0.dosage.txt
#download and upload to google cloud to run


##export as plink if want to use plink to run GWAS
use PLINK2
#for original file
for i in {1..22}; do \
qsub -l h_vmem=100G wrapper1.sh plink --vcf AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr$i.noINFO1.vcf --make-bed --out AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr$i.noINFO1  ;done

#and the partial genomes for the two ancestries
for i in {1..22}; do \
qsub -l h_vmem=100G wrapper1.sh plink --vcf AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr$i.noINFO1.anc0.vcf --vcf-half-call haploid --out AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr$i.noINFO1.anc0 --make-bed  ;done

for i in {1..22}; do \
qsub -l h_vmem=100G wrapper1.sh plink --vcf AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr$i.noINFO1.anc1.vcf --vcf-half-call haploid --out AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr$i.noINFO1.anc1 --make-bed ;done

##merge chromosomes together into one file for easier use
#needed to remove tri-allelic sites first
qsub -l h_vmem=100G wrapper1.sh plink --bfile AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr1.noINFO2.anc0 --merge-list autosomefiles1.anc0 --make-bed --out UKBB_AfEur_QCed_lipids.auto.anc0 --noweb 

qsub -l h_vmem=100G wrapper1.sh plink --bfile AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr1.noINFO2.anc1 --merge-list autosomefiles1.anc1 --make-bed --out UKBB_AfEur_QCed_lipids.auto.anc1 --noweb 

qsub -l h_vmem=100G wrapper1.sh plink --bfile AdmixedAfrEur_Results_UKBB_AfEur_QCed_lipids.chr1.noINFO2 --merge-list autosomefiles1 --make-bed --out UKBB_AfEur_QCed_lipids.auto --noweb 

#run GWAS on the original file
cat phenolist | while read line ; do \
qsub -l h_vmem=100G -l h_rt=48:00:00 -N GWAS$line wrapper1.sh /broad/software/free/Linux/redhat_6_x86_64/pkgs/plink_1.07/bin/plink --noweb --bfile UKBB_AfEur_QCed_lipids.auto --pheno AdmixedAfAm_lipid_phenos --linear --covar AdmixedAfAm_lipidcovs --pheno-name $line --allow-no-sex --covar-name isFemale,age,dilutionFactor --hide-covar --out $line.CI --ci 0.95 ; done

##run GWAS on the two anc halves
cat phenolist | while read line ; do \
qsub -l h_vmem=100G -N GWAS$line.anc0 -l h_rt=40:00:00 wrapper1.sh /broad/software/free/Linux/redhat_6_x86_64/pkgs/plink_1.07/bin/plink --noweb --bfile UKBB_AfEur_QCed_lipids.auto.anc0 --pheno AdmixedAfAm_lipid_phenos --linear --covar AdmixedAfAm_lipidcovs --pheno-name $line --allow-no-sex --covar-name isFemale,age,dilutionFactor --hide-covar --out $line.anc0.CI --ci 0.95 ; done

cat phenolist | while read line ; do \
qsub -l h_vmem=100G -N GWAS$line.anc1 -l h_rt=40:00:00 wrapper1.sh /broad/software/free/Linux/redhat_6_x86_64/pkgs/plink_1.07/bin/plink --noweb --bfile UKBB_AfEur_QCed_lipids.auto.anc1 --pheno AdmixedAfAm_lipid_phenos --linear --covar AdmixedAfAm_lipidcovs --pheno-name $line --allow-no-sex --covar-name isFemale,age,dilutionFactor --hide-covar --out $line.anc1.CI --ci 0.95 ; done

#meta-analyse the anc0 and anc1 files
cat phenolist | while read line ; do \
qsub -l h_vmem=100G -l h_rt=10:00:00 wrapper1.sh /broad/software/free/Linux/redhat_6_x86_64/pkgs/plink_1.07/bin/plink --noweb --meta-analysis  $line.anc0.CI.assoc.linear $line.anc1.CI.assoc.linear + qt --out $line ;done


#plot up results - move to R. Wrote up a script to plot this up with qqman package - PlotManQQ.R

cat phenolist1 | while read line ; do \
PlotManQQ.R 


