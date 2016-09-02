setwd('~/d/sci/src/freq_filter')
options(stringsAsFactors=FALSE)
source('../exac_2015/exac_constants.R')

# load ExAC dataset
if (!("exac" %in% ls(globalenv()))) {
  exac = load_exac_data(reload=T)
}

clinvar = read.table('../clinvar/output/clinvar.tsv',sep='\t',header=T,quote='',comment.char='')
clinvar$pos_id = paste(clinvar$chrom, formatC(clinvar$pos,width=9,flag='0'), clinvar$ref, clinvar$alt, sep='_')
clinvar$exac_ac_adj = exac$ac_adj[match(clinvar$pos_id, exac$pos_id)]
clinvar$exac_ac_adj[is.na(clinvar$exac_ac_adj)] = 0
# "Literature only" refers to ClinVar submissions that have no submitter except OMIM and/or GeneReviews
clinvar$lit_only = nchar(gsub('OMIM','',gsub('GeneReviews','',gsub(';','',clinvar$all_submitters)))) == 0

hcm = grepl("cardiomyopathy",clinvar$all_traits,ignore.case=TRUE) & grepl("hypertrophic",clinvar$all_traits,ignore.case=TRUE)

sum(hcm) # 1063
sum(hcm & clinvar$pathogenic==1) # 631
sum(hcm & clinvar$exac_ac_adj > 9 & clinvar$pathogenic==1) # 50
sum(hcm & clinvar$exac_ac_adj > 9 & clinvar$pathogenic==1 & clinvar$all_pmids != '') # 50
sum(hcm & clinvar$exac_ac_adj > 9 & clinvar$pathogenic==1 & clinvar$lit_only) # 2
sum(hcm & clinvar$exac_ac_adj > 9 & clinvar$pathogenic==1 & clinvar$conflicted==0) # 28
sum(hcm & clinvar$exac_ac_adj < 9 & clinvar$pathogenic==1) # 577
sum(hcm & clinvar$exac_ac_adj < 9 & clinvar$pathogenic==1 & clinvar$conflicted==0) # 573

# create a table of these variants to browse in Rstudio
hcm_prev = clinvar[hcm & clinvar$exac_ac_adj > 9 & clinvar$pathogenic==1 & clinvar$conflicted==0,]

# write out a table of HCM variants per Nicky's request
write.table(clinvar[hcm,],'data/clinvar_hcm.tsv',sep='\t',row.names=F,col.names=T,quote=F)

pcd = grepl("primary",clinvar$all_traits,ignore.case=TRUE) & grepl("ciliary",clinvar$all_traits,ignore.case=TRUE) & grepl("dyskinesia",clinvar$all_traits,ignore.case=TRUE)
sum(pcd) # 246
sum(pcd & clinvar$pathogenic==1 & clinvar$conflicted==0) # 128
sum(pcd & clinvar$pathogenic==1 & clinvar$conflicted==0 & clinvar$exac_ac_adj > 154) # 1

nccm = grepl("cardiomyopathy",clinvar$all_traits,ignore.case=TRUE) & grepl("compaction",clinvar$all_traits,ignore.case=TRUE)
sum(nccm) # 85
sum(nccm & clinvar$pathogenic==1 & clinvar$conflicted==0) # 42
sum(nccm & clinvar$pathogenic==1 & clinvar$conflicted==0 & clinvar$exac_ac_adj > 18) # 3

to_look_at = clinvar$pathogenic==1 & clinvar$conflicted==0 & ((nccm & clinvar$exac_ac_adj > 18) | (pcd & clinvar$exac_ac_adj > 154))
write.table(clinvar[to_look_at,],'data/clinvar_pcd_nccm_high_af.tsv',sep='\t',row.names=F,col.names=T,quote=F)
