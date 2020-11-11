#!/usr/bin/env python
# coding: utf-8

# # Initialize

# In[1]:


import argparse
import hail as hl
import numpy as np
hl.init()


# In[2]:


from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()


# # Load in the filtered and phenotype annotated genotype data for the individuals
# ### Note: the example file is in Matrix Table format, the native Hail format. VCF formats may also be imported with 'hl.import_vcf'. See documentation for more details: https://hail.is/

# In[3]:


#key columns by sample and rows by rsID for easier merge with dosage data
mt = hl.read_matrix_table('gs://ukb-diverse-pops/AdmixedAfrEur/Results/UKBB_AfEur_QCed_lipids2.mt').key_rows_by('locus')


# In[4]:


#get out just the sample pheno info I want
mt = mt.select_cols(mt.pheno.isMale, mt.pheno.age, mt.pheno.dilutionFactor, mt.pheno.TC, mt.pheno.HDLC, mt.pheno.LDLC, mt.pheno.TG, mt.sqc1.PC1, mt.sqc1.PC2,
                mt.sqc1.PC3, mt.sqc1.PC4, mt.sqc1.PC5, mt.sqc1.PC6, mt.sqc1.PC7, mt.sqc1.PC8, mt.sqc1.PC9, 
                     mt.sqc1.PC10, mt.sqc1.PC11, mt.sqc1.PC12, mt.sqc1.PC13, mt.sqc1.PC14, mt.sqc1.PC15, mt.sqc1.PC16, mt.sqc1.PC17, mt.sqc1.PC18, mt.sqc1.PC19, mt.sqc1.PC20)


# In[5]:


#clear out the genotype information that we no longer need
mt = mt.select_entries()


# #export sample covariates for comparison running in plink
mt.cols().export('gs://ukb-diverse-pops/AdmixedAfrEur/Results/AfrEur_samples_covariates.tsv')

# # Load in the dosage files from Tractor
# ## First implementing the GWAS version that includes the full VCF but with dosage calls per ancestry

# In[6]:


row_fields={'CHROM': hl.tstr, 'POS': hl.tint, 'ID': hl.tstr, 'REF': hl.tstr, 'ALT': hl.tstr} 
anc0dos = hl.import_matrix_table('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids.autosomes.anc0.dosage_v1.txt.gz', 
                                 force_bgz=True, row_fields=row_fields, row_key=[], min_partitions=32) 
anc0dos = anc0dos.key_rows_by().drop('row_id')
anc0dos = anc0dos.key_rows_by(locus=hl.locus(anc0dos.CHROM, anc0dos.POS))     


# In[7]:


row_fields={'CHROM': hl.tstr, 'POS': hl.tint, 'ID': hl.tstr, 'REF': hl.tstr, 'ALT': hl.tstr} 
anc1dos = hl.import_matrix_table('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids.autosomes.anc1.dosage_v1.txt.gz', 
                                 force_bgz=True, row_fields=row_fields, row_key=[], min_partitions=32) 
anc1dos = anc1dos.key_rows_by().drop('row_id')
anc1dos = anc1dos.key_rows_by(locus=hl.locus(anc1dos.CHROM, anc1dos.POS))     


# In[8]:


#save these temporary files to relieve  memory burden
anc1dos = anc1dos.checkpoint('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids.autosomes.anc1.dosage_v1.mt')


# In[9]:


anc0dos = anc0dos.checkpoint('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids.autosomes.anc0.dosage_v1.mt')


# In[10]:


anc0dos.show()


# In[11]:


#join the dosage files to the genotype data. 
#Specifically, annotating the samples with their info for how many copies of the minor allele per ancestry were seen
mt = mt.annotate_entries(anc0dos = anc0dos[mt.locus, mt.s], anc1dos = anc1dos[mt.locus, mt.s])


# In[13]:
##write out this file for future use

mt.write('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids_dosages.mt')


# In[14]:


##add in the haplotype counts per site per indiv also

#read in haplotype counts for anc0, African
row_fields={'CHROM': hl.tstr, 'POS': hl.tint, 'ID': hl.tstr} 
hapcounts0 = hl.import_matrix_table('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids.auto.anc0.hapcount.txt.gz', 
                                 force_bgz=True, row_fields=row_fields, row_key=[], min_partitions=32) 
hapcounts0 = hapcounts0.key_rows_by().drop('row_id')
hapcounts0 = hapcounts0.key_rows_by(locus=hl.locus(hapcounts0.CHROM, hapcounts0.POS))   


# In[26]:


withinPCs.describe()


# In[27]:


withinPCs.show(1)


# In[21]:


hapcounts0.show(1)


# In[22]:


#annotate the dosage mt with within sample PCs and haplotype counts
mt = mt.annotate_entries(hapcounts0 = hapcounts0[mt.locus, mt.s])


# In[30]:


mt = mt.annotate_cols(withinPCs = withinPCs[mt.s])


# In[31]:


mt.describe()


# In[32]:


#export this new file with the haplotype counts and within sample PCs added in
mt.write('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids_dosages_hapcount.mt')


# In[8]:


anc0dos1 = anc0dos.annotate_cols(pheno_info=mt.cols()[anc0dos.col_id])


# In[10]:


anc1dos1 = anc1dos.annotate_cols(pheno_info=mt.cols()[anc1dos.col_id])


# # Run linear regression across pheno using the Tractor method

# In[3]:

##load mt back in for memory savings
mt = hl.read_matrix_table('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids_dosages_hapcount.mt')


# In[37]:


mt.show(1)


# In[5]:


admixFrac = hl.import_table('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKB_admixed.1kG.2.Q.txt').key_by('IID')


# In[20]:


mt.describe()


# In[7]:


admixFrac.show(2)
#will need to convert fractions to integers when run the GWAS


# In[8]:


#annotate in with the ADMIXTURE fractions to run instead of PCs
mt = mt.annotate_cols(admixFrac = admixFrac[mt.s])



#run all 3 lipid phenotypes in batch with the ADMIXTURE fract as global ancestry covariate
phenonames = ['TC', 'HDLC', 'LDLC']

#scale up to run on batch of phenotypes, global UKB PCs
mt = mt.annotate_rows(results=hl.struct(**{pheno: hl.agg.linreg(mt[pheno], 
                                                                [1.0, mt.hapcounts0.x, mt.anc0dos.x, mt.anc1dos.x, mt.isMale, mt.age, mt.dilutionFactor, hl.float(mt.admixFrac.AFR)]) 
                                           for pheno in phenonames}))


# In[19]:


#export the final file with the annotated results for just ADMIXTURE rather than PCs
mt.write('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids_dosages_admixfrac.mt')


# In[21]:


#then can read the mt file back in and plot up the manhattan and QQ plots
mt = hl.read_matrix_table('gs://ukb-diverse-pops/AdmixedAfrEur/DosageFiles/UKBB_AfEur_QCed_lipids_dosages_admixfrac.mt')


# Note - for plotting, will need to use workers rather than premptible workers

# In[23]:


#plot up TC as a manhattan plot, anc 0
p_TC0 = hl.plot.manhattan(mt.results.TC.p_value[2], title='Admixed Afr-Eur UKBB, TC, anc0', collect_all=False, significance_line=5e-08) 
#colors=["#030303", "#7F7F7F"])
show(p_TC0)


# In[24]:


#plot up TC manhattan plot, anc 1
p_TC1 = hl.plot.manhattan(mt.results.TC.p_value[3], title='Admixed Afr-Eur UKBB, TC, anc1', collect_all=False, significance_line=5e-08) 
#colors=["#030303", "#7F7F7F"])
show(p_TC1)


# In[25]:


#make a QQ plot for TC anc0
p = hl.plot.qq(mt.results.TC.p_value[2], title="QQ plot, TC, anc0")
show(p)


# In[26]:


#make a QQ plot for TC anc1
p = hl.plot.qq(mt.results.TC.p_value[3], title="QQ plot, TC, anc1")
show(p)


# In[33]:


#plot up HDLC manhattan plot, anc 0 then anc1
p = hl.plot.manhattan(mt.results.HDLC.p_value[2], title='Admixed Afr-Eur UKBB, HDLC, anc0', collect_all=False, significance_line=5e-08) 
show(p)


# In[ ]:


p = hl.plot.manhattan(mt.results.HDLC.p_value[3], title='Admixed Afr-Eur UKBB, HDLC, anc1', collect_all=False, significance_line=5e-08) 
#colors=["#030303", "#7F7F7F"])
show(p)


# In[35]:


#plot up LDLC manhattan plot, anc 0 then anc1
p = hl.plot.manhattan(mt.results.LDLC.p_value[2], title='Admixed Afr-Eur UKBB, LDLC, anc0', collect_all=False, significance_line=5e-08) 
show(p)


# In[36]:


p = hl.plot.manhattan(mt.results.LDLC.p_value[3], title='Admixed Afr-Eur UKBB, LDLC, anc1', collect_all=False, significance_line=5e-08) 
#colors=["#030303", "#7F7F7F"])
show(p)


# In[29]:


#make a QQ plots for HDLC and LDLC
p = hl.plot.qq(mt.results.HDLC.p_value[3], title="QQ plot, HDLC, anc1")
show(p)


# In[37]:


p = hl.plot.qq(mt.results.HDLC.p_value[2], title="QQ plot, HDLC, anc0")
show(p)


# In[38]:


p = hl.plot.qq(mt.results.LDLC.p_value[3], title="QQ plot, LDLC, anc1")
show(p)


# In[39]:


p = hl.plot.qq(mt.results.LDLC.p_value[2], title="QQ plot, LDLC, anc0")
show(p)


# In[48]:


mt1 = mt.rows()


# In[49]:


mt1.describe()


# In[50]:


#export results as a text file to plot up offline if need be
mt1.export('gs://ukb-diverse-pops/AdmixedAfrEur/Results/mt.results.tsv.bgz')
