setwd('~/d/sci/src/freq_filter')
options(stringsAsFactors=FALSE)
source('../exac_2015/exac_constants.R')

# load ExAC dataset
if (!("exac" %in% ls(globalenv()))) {
  exac = load_exac_data(reload=F)
}

# forward method by James Ware:
# given a desired AF to filter on for Mendelian analyses, what is the max AC that would be expected in a dataset?
find_max_ac = function(af,an,ci=.95) {
  if (af == 0) {
    return (0)
  } else {
    quantile_limit = ci # ci for one-sided, 1-(1-ci)/2 for two-sided
    max_ac = qpois(quantile_limit,an*af)
    return (max_ac)
  }
}


# backward method: given an observed AC, what is the highest AF filter for which one should remove this variant?
# This function will accept an AC_Adj and an AN_Adj and will return the highest AF filter for which you would want to reject
# this variant as a potential Mendelian causal allele.# for uniroot to work, lower has to be rarer than a singleton, and upper has to be higher than fixed
# Important notes about parameters:
# 1. for uniroot to work, lower has to be rarer than a singleton, and upper has to be higher than fixed (100% AF), 
#    hence the AF here runs from 0.1 alleles in ExAC to an impossible AF of 200%
# 2. you need a tight tolerance (tol) to get accurate results from uniroot at low AF, so we use 1e-7 as default. 
# 3. we also specify a precision (1e-6) for incrementing to find the upper boundary of the range of AFs for which
#    the 95%CI AC is less than the observed AC
find_af_filter = function(ac, an, ci=.95, lower=(.1/(2*60706)), upper=2, tol=1e-7, precision=1e-6) { 
  if (is.na(ac) | is.na(an) | ac == 0 | an == 0 | ac == 1) {
    return (0.0)
  } else {
    quantile_limit = ci # ci for one-sided, 1-(1-ci)/2 for two-sided
    # this has been tested and with the default parameters here, uniroot seems to never fail for an AC, AN values
    # in ExAC. still, it's good to have this tryCatch in here for debugging or in case anyone wants to override
    # the defaults. instead of just giving the error/warning message it also tells you which AC and AN values
    # it failed on.
    attempt_uniroot = tryCatch({
      uniroot_result = uniroot(f = function(af,ac,an) { return (ac - 1 - qpois(p=quantile_limit,lambda=an*af)) },lower=lower,upper=upper,ac=ac,an=an,tol=tol)
    }, warning = function(w) {
      print(paste("ac= ",as.character(ac),", an= ",as.character(an)," warning = ",as.character(w),sep=''))
      return (0.0)
    }, error = function(e) {
      print(paste("ac= ",as.character(ac),", an= ",as.character(an)," error = ",as.character(e),sep=''))
      return (0.0)
    }, finally = {
    })
    max_af = round(uniroot_result$root,-log10(precision)) # round to nearest millionth
    while(find_max_ac(af=max_af,an=an) < ac) {
      max_af = max_af + precision # increase by millionths until you find the upper bound - the highest AF for which 95%CI AC is still less than observed AC
    }
    max_af = max_af - precision # back off one unit from the AF that violated the 95%CI AC < obs AC condition
    return (max_af)
  }
}

# This finds the af_filter for, say, a vector of AC at 10,000 differnet alleles and a corresponding vector of 10,000 AN values
# it's basically just a wrapper for mapply, which is trivial but makes more readable the code below in find_highest_af_filter.
find_af_filter_vectorized = function(ac_vector, an_vector, ci=.95, lower=(.1/(2*60706)), upper=2, tol=1e-7) { 
  return (mapply(find_af_filter, ac_vector, an_vector, ci=ci, lower=lower, upper=upper, tol=tol))
}

# In some cases, the highest af_filter value will _NOT_ come from the popmax population, because
# sometimes the population with higher AN_Adj will give more filtering power even if its AF is slightly lower
# therefore, it is important to be able to take the highest af_filter value from a vector of population AC and AN
# this function will accept a list of vectors of AC and of AN. for instance, ac_list may be a list of 5 vectors, which are
# each 10,000 value of ac_afr, ac_amr, ac_eas, ac_nfe, and ac_sas. an_list is the corresponding list of an_ vectors.
find_highest_af_filter = function(ac_list,an_list,ci=.95) {
  len_ac = length(list(ac_list))
  len_an = length(list(an_list))
  stopifnot(len_ac==len_an) # check that we have the same number of AC and AN columns
  # make a list of af_filter vectors for each population
  # SIMPLIFY=FALSE makes it return a list instead of concatenating the vectors
  af_filter_list = mapply(find_af_filter_vectorized,ac_list,an_list,SIMPLIFY=FALSE) 
  # now take pmax of the population-specific af_filters, returning ONE vector
  # containing the highest af_filter value across populations
  highest_af_filter = do.call(pmax, af_filter_list) 
  highest_af_filter = numeric(length(af_filter_list[[1]]))
  highest_af_filter_pop = character(length(af_filter_list[[1]]))
  for (i in 1:length(highest_af_filter)) { # couldn't figure out a way to do both things with do.call so just looping
    to_compare = mapply("[[",af_filter_list,i)
    idx = which.max(to_compare)
    highest_af_filter[i] = to_compare[idx] # get max filter
    highest_af_filter_pop[i] = names(to_compare)[idx] # get the pop that gave the max filter
  }
  highest_af_filter_pop = gsub("ac_","",highest_af_filter_pop)
  highest_af_filter_pop[highest_af_filter==0] = 'none'
  return (list(af_filter=highest_af_filter, af_filter_pop=highest_af_filter_pop))
}

# QC test #1: do we get the same answer forward and backward?
an = 2*60706
find_af_filter(find_max_ac(af=.001,an=an),an=an)
find_af_filter(find_max_ac(af=.01,an=an),an=an)

# QC test 2. do we get roughly the same answer across a whole range of AC?
an = 2*60706 # for now assume full genotyping
afs = c((1:9)*10^-6,(1:9)*10^-5,(1:9)*10^-4,(1:9)*10^-3,1*10^-2)
plot(NA,NA,xlim=range(afs),ylim=c(1,300),xlab='AF',ylab='95%CI max AC',log='xy',axes=FALSE)
forward_ac = numeric(length(afs))
backward_af = numeric(length(afs))
for (i in 1:length(afs)) {
  forward_ac[i] = find_max_ac(afs[i],an=an)
  backward_af[i] = find_af_filter(ac=forward_ac[i],an=an)
}
points(x=afs,y=forward_ac,col='red',type='s',lwd=10)
points(x=backward_af,y=forward_ac,col='blue',type='s',lwd=3)
legend('topleft',c('forward (James)','backward (Eric)'),col=c('red','blue'),lwd=c(10,3))
axis(side=1,at=10^(-6:-2),labels=percent(10^(-6:-2)),lwd=0,lwd.ticks=1)
axis(side=2,at=(0:6)*50,labels=(0:6)*50,lwd=0,lwd.ticks=1,las=1)

# QC test 3. is the backward method sensitive to changes in AC at the lower end - singleton, doubleton etc.
find_af_filter(ac=1,an=2*60706,ci=.95)
find_af_filter(ac=2,an=2*60706,ci=.95)
find_af_filter(ac=3,an=2*60706,ci=.95)
find_af_filter(ac=4,an=2*60706,ci=.95)

# QC test 4. does the backward method compute in reasonable time?
# (note just using ac_popmax and an_popmax for this initial test, though we'll ultimately want to find the
# highest af_filter across all populations)
exac_test = head(exac,n=200000)
start_time = proc.time()
mapply(find_af_filter, exac_test$ac_popmax, exac_test$an_popmax, ci=.95)
proc.time() - start_time
# takes roughly in the range of ~200 seconds for 200K variants

# QC test 5. for an actual set of ExAC variants, does running backwards and then forwards reproduce the actual AC reaosnably well?
exac_test$af_filter = mapply(find_af_filter, ac=exac_test$ac_popmax, an=exac_test$an_popmax, ci=.95)
exac_test$ac_forward = mapply(find_max_ac, af=exac_test$af_filter, an=exac_test$an_popmax, ci=.95)
cor.test(exac_test$ac_popmax, exac_test$ac_forward)
# yes, perfectly

# QC test 6. does the AF filter compare to the actual AF popmax in a reasonable way
#png('display/af_popmax_vs_af_filter.png',height=600,width=600)
par(mar=c(5,5,2,2))
plot(exac_test$af_popmax, exac_test$af_filter, pch=20, xlab='', ylab='',
     log='xy', xlim=10^c(-6,-1), ylim=10^c(-6,-1), axes=FALSE)
abline(a=0,b=1,col='red')
axis(side=1,at=10^(-5:-1),labels=percent(10^(-5:-1)),lwd=0,lwd.ticks=1)
axis(side=2,at=10^(-5:-1),labels=percent(10^(-5:-1)),lwd=0,lwd.ticks=1,las=1)
mtext(side=1,line=3,text='AF popmax',font=2)
mtext(side=2,line=4,text='AF filter',font=2)
#dev.off()
# yes, all points are below the diagonal, and the correlation is fairly tight
# do a sanity check that much of the noisiness can be explained by differences in AN:
find_af_filter(ac=100,an=100000)
find_af_filter(ac=10,an=10000)
find_af_filter(ac=1,an=1000)
# ok, good...

# QC test 7. does the find highest AF filter function work correctly and give the same answer 
# when called two different ways?
# first, call with pmax and separate mapply calls
exac_test$af_filter_max = pmax(mapply(find_af_filter,ac=exac_test$ac_afr,an=exac_test$an_afr),
                               mapply(find_af_filter,ac=exac_test$ac_amr,an=exac_test$an_amr),
                               mapply(find_af_filter,ac=exac_test$ac_eas,an=exac_test$an_eas),
                               mapply(find_af_filter,ac=exac_test$ac_nfe,an=exac_test$an_nfe),
                               mapply(find_af_filter,ac=exac_test$ac_sas,an=exac_test$an_sas))

# now call with find_highest_af_filter which bundles all of that together
ptm = proc.time()
test_filtering_data = find_highest_af_filter(ac_list=as.list(exac_test[,c("ac_afr","ac_amr","ac_eas","ac_nfe","ac_sas")]),
                                          an_list=as.list(exac_test[,c("an_afr","an_amr","an_eas","an_nfe","an_sas")]))
exac_test$af_filter_max_2 = test_filtering_data[["af_filter"]]
exac_test$af_filter_pop = test_filtering_data[["af_filter_pop"]]
proc.time() - ptm # 310 seconds for 200K variants
# this predicts just over 4 hours to compute af_filter for all of exac

# check that the correlation between these two methods is perfect
cor.test(exac_test$af_filter_max, exac_test$af_filter_max_2)
# it is.

# and do the results still look reasonable w.r.t. af_popmax?
plot(exac_test$af_popmax, exac_test$af_filter_max_2,pch=20)
# yep

# ok, next actually compute af_filter for all of exac (N.B. this block took about 3.5 hours to complete)
filtering_data = find_highest_af_filter(ac_list=as.list(exac[,c("ac_afr","ac_amr","ac_eas","ac_nfe","ac_sas")]),
                                        an_list=as.list(exac[,c("an_afr","an_amr","an_eas","an_nfe","an_sas")]))
exac$af_filter = filtering_data[["af_filter"]]
exac$af_filter_pop = filtering_data[["af_filter_pop"]]

# write out results for public release
af_filter_data = exac[,c("chrom","pos","ref","alt","af_filter","af_filter_pop")]
write.table(af_filter_data,"af_filter_data.tsv",sep='\t',row.names=F,col.names=T,quote=F)


# a few QC checks.
# most variants should fall on the diagonal where popmax == af_filter_pop
table(exac[,c("popmax","af_filter_pop")])
# af_filter_pop
# popmax  afr     amr     eas     nfe    none     sas
# AFR  786392     202      56  101735  807486    6256
# AMR    1804  387449     181   83244  589595    7134
# EAS    1073     641  488158   56861  648524    7575
# NFE      24       1       9 1273470 2499020    1332
# SAS     139      50      29   61030 1147140  764055
# ok, that looks good.
# And the correlation should be tight:
cor.test(exac$af_filter, exac$af_popmax)
# which it is: cor = 0.9881759

### FIGURE

# variants remaining after filters
# things to loop over in generating summary stats
# allele frequency thresholds
af_limits = c(1e-6, 10^-5.5, 1e-5, 10^-4.5, 1e-4, 10^-3.5, 1e-3, 10^-2.5, 1e-2, 5e-2)
ancestries = c(names(pop_names)[c(1,2,3,5,7)]) # "afr" "amr" "eas" "nfe" "sas"
filters = c("af_filter")
# only doing dominant analysis so don't need this line:
# zygosities = c('any','hom')

# We don't have permission to release individual genotype-level data for the "leave-out 500" for the
# simulated Mendelian analysis, but such data are needed in order to compute the filtering statistics
# for Figure 3A. Thus, this part of the code can only be reproduced in-house. recalculate_stats is
# therefore set to FALSE by default for public release. We who are running the code locally can set to
# TRUE to recompute the filtering statistics. Note that this involves re-computing the af_filter
# on the 60,206 individuals who aren't left out, plus a bunch of costly joins between ExAC and the lv500
# dataset. Thus, this block of code takes several hours and should be run overnight.
recalculate_stats = FALSE
if ('lv500.table.gz' %in% list.files('../exac_papers/misc_data/') & recalculate_stats) {
  lv = read.table('../exac_papers/misc_data/lv500.table.gz',sep='\t',header=TRUE,comment.char='',quote='')
  colnames(lv) = tolower(colnames(lv))
  lv = lv[,-which(colnames(lv)=='ac_hemi')] 
  
  # join relevant cols to ExAC
  lv$pos_id = paste(lv$chrom, formatC(lv$pos,width=9,flag='0'), lv$ref, lv$alt, sep='_')
  exac$pos_id = paste(exac$chrom, formatC(exac$pos,width=9,flag='0'), exac$ref, exac$alt, sep='_')
  # figure out how to match lv500 to ExAC
  match_indices = match(lv$pos_id, exac$pos_id)
  
  # compute AC and AN for the 60,206 non-left-out individuals
  for (colname in colnames(lv)[7:26]) {
    # compute ACs in ExAC minus the leave-out 500
    lv[,paste("exac",colname,sep="_")] = exac[match_indices,colname] - lv[,colname]
  }
  
  # re-calculate af_filter without the leave-out individuals in there
  lv$af_filter = find_highest_af_filter(ac_list=as.list(lv[,c("exac_ac_afr","exac_ac_amr","exac_ac_eas","exac_ac_nfe","exac_ac_sas")]),
                                        an_list=as.list(lv[,c("exac_an_afr","exac_an_amr","exac_an_eas","exac_an_nfe","exac_an_sas")]))[["af_filter"]]
  
  # copy `use` and `category` from exac (and only use if lv500 also has ac_adj > 0)
  lv$use = exac$use[match_indices] & lv$ac_adj > 0
  lv$category = exac$category[match_indices]
  
  # some constants and stats
  n_indiv = {} # count the number of ExAC individuals in each continental ancestry
  for (ancestry in ancestries) {
    n_indiv[[ancestry]] = max(lv[,paste("an",ancestry,sep="_")])/2
  }
  sum(n_indiv) # this will be 500 
  
  # now calculate summary stats for Figure 3A
  stats = expand.grid(ancestries,filters,af_limits)
  colnames(stats) = c('ancestry','filter','af_limit')
  stats$ancestry = as.character(stats$ancestry)
  stats$total_vars = 0
  stats$filtered_vars = 0
  stats$proportion_removed = 0.0
  for (af_limit in af_limits) {
    for (filter in filters) {
      # only consider "use" variants that are missense, missense-equivalent, or protein-truncating
      include_in_total = lv$use & !is.na(lv$category) & lv$category %in% c(mis_like,lof_like)
      # note the < logic in next line. because af_filter is an af for which 95%CI AC < obs AC, 
      # users should remove variants whose af_filter is _greater than or equal to_ their maximum credible AF
      # conversely, the "filtered" set of variants includes only those with af_filter strictly _less than_
      # the user's maximum credible AF
      include_in_filtered = include_in_total & (lv[,filter] < af_limit) 
      for (ancestry in ancestries) {
        # find target row of stats table
        row = stats$ancestry == ancestry & stats$filter == filter & stats$af_limit == af_limit 
        exac_colname = paste("ac",ancestry,sep="_") # column of interest in exac
        # below, these n_indiv terms cancel in the final calculation, but are included as a reminder
        # that conceptually what we are computing is the average number of variants _in each person_
        # of a given ancestry that are removed by an allele frequency filter
        average_variants_without_af_filter = sum(lv[include_in_total,exac_colname]) / n_indiv[ancestry]
        average_variants_with_af_filter = sum(lv[include_in_filtered,exac_colname]) / n_indiv[ancestry] 
        proportion_removed = 1 - average_variants_with_af_filter / average_variants_without_af_filter
        stats$proportion_removed[row] = proportion_removed # store it in the table
        stats$total_vars[row] = average_variants_without_af_filter
        stats$filtered_vars[row] = average_variants_with_af_filter
      }
    }
  }
  # write out results so that others can at least reproduce Figure 3A, given this table of stats
  write.table(stats,'data/filtering_stats.tsv',sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
} else {
  # if not re-calculating stats, read them in and then you can reproduce Figure 3A
  stats = read.table('data/filtering_stats.tsv',sep='\t',quote='',header=TRUE)
}

# get colors and x axis values for plotting Figure 3A. if you end up with a Figure 3A that lacks
# the curves themselves and has just the axes, it's probably because you forgot to run this line
stats$acol = pop_colors[stats$ancestry]

# this saves a version of Figure 3A with gridlines. The actual version of Figure 3A for publication
# is produced in analysis/Figures.Rmd as a multi-panel plot with Figure 3B, from a modified version of this code.
png('figures/variants_after_filtering.png',width=800,height=500,res=100)
par(mar=c(4,5,1,1))
plot(NA, NA, xlim=c(-2,-6), ylim=c(0,550), ann=FALSE, axes=FALSE, yaxs='i')
abline(h=(0:5)*100, lwd=.5)
for (ancestry in ancestries) {
  rows = stats$ancestry==ancestry & stats$af_limit < .05
  points(x=log10(stats$af_limit[rows]), y=stats$filtered_vars[rows], type='l', lwd=5, col=stats$acol[rows])
}
axis(side=1, at=-2:-6, labels=c("1%","0.1%","0.01%","1e-5","1e-6"), lwd=0, lwd.ticks=1)
axis(side=2, at=(0:5)*100, labels=(0:5)*100, lwd=0, lwd.ticks=0, las=2)
mtext(side=1, line=2.5, text='your maximum credible population AF')
mtext(side=2, line=3.5, text='mean predicted protein-altering variants per exome')
legend('topright',pop_names[ancestries],col=pop_colors[ancestries],text.font=2,text.col=pop_colors[ancestries],lwd=5)
dev.off()

# statistics quoted in text -- averages across populations at different filtering thresholds
mean(stats$filtered_vars[stats$af_limit==.001])
mean(stats$filtered_vars[stats$af_limit==1e-6])


