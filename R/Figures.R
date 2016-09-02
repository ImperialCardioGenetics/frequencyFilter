#setwd("./R")
#source("../R/setup.R")
my.CI<-0.95

### Figure 2

variants<-read.table("../data-raw/LMM_OXF_path_likelyPath_VUS_variants_allExACpopulations_withLowestClass.txt", header=TRUE, sep="\t")
maxHCMPenTerm<-9
maxHCM<-5
variants$Class <- factor(variants$Class, c("Pathogenic", "Likely Pathogenic", "VUS"))

fullplot <- ggplot(subset(variants, (Disease=="HCM" & !(Gene %in% c('GLA','DMD','TAZ','EMD','LAMP2','VBP1','FHL1','KCNE1L','ZIC3','COL4A5','FLNA','MED12')))), aes(x=CaseCount, y=ExACCountALL, colour=Class, shape=Class)) +
  geom_point(size=2, alpha=0.6, position=position_jitter(width=0.5,height=0.5)) +
  geom_hline(yintercept=maxHCMPenTerm, linetype=2, colour="dodgerblue4") +
  geom_hline(yintercept=maxHCM, linetype=2, colour="lightskyblue3") +
  scale_color_manual(values=c("#ED1E24", "#FF9912", "gray38")) +
  scale_shape_manual(values=c(19,19,1)) +
  scale_y_continuous(limits = c(-2,178), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-2,108), expand = c(0, 0)) +
  theme(axis.text=element_text(size=12),axis.title=element_blank(),legend.position = c(.7, .8), 
        legend.text=element_text(size=14), legend.title = element_blank(), legend.key = element_blank(), 
        axis.line.x = element_line(color="grey"), axis.line.y = element_line(color="grey"), panel.background = element_blank(), 
        axis.ticks = element_line(color="grey"), plot.background = element_rect(linetype="solid", color="grey"))

zoomplot <- ggplot(subset(variants, (Disease=="HCM" & !(Gene %in% c('GLA','DMD','TAZ','EMD','LAMP2','VBP1','FHL1','KCNE1L','ZIC3','COL4A5','FLNA','MED12')))), aes(x=CaseCount, y=ExACCountALL, colour=Class, shape=Class)) +
  geom_point(size=2, alpha=0.6, position=position_jitter(width=0.5,height=0.5)) +
  geom_hline(yintercept=(maxHCMPenTerm+0.2), linetype=2, colour="dodgerblue4") +
  geom_hline(yintercept=(maxHCM+0.2), linetype=2, colour="lightskyblue3") +
  scale_color_manual(values=c("#ED1E24", "#FF9912", "gray38")) +
  scale_shape_manual(values=c(19,19,1)) +
  labs(x="Count cases",y="Count ExAC") +
  #ylim(0,30) +
  #xlim(0,70) +
  scale_y_continuous(limits = c(-0.5,30), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0,70), expand = c(0, 0)) +
  annotate("text",x=67,y=(maxHCMPenTerm+0.6), label="AC = 9",colour="dodgerblue4") +
  annotate("text",x=67,y=(maxHCM+0.6), label="AC = 5",colour="lightskyblue3") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        legend.position="none",  axis.line.x = element_line(color="grey"), axis.line.y = element_line(color="grey"), panel.background = element_blank(), 
        axis.ticks = element_line(color="grey"))

png("../figures/Figure2.png",width=1200,height=1000,res=150)
subvp<-viewport(width=.5,height=.5,x=.745,y=.745)

zoomplot
print(zoomplot)
print(fullplot,vp=subvp)
dev.off()

### Figure 3

# Details for Exome plot (figure 3a)
stats = read.table('../data/filtering_stats.tsv',sep='\t',quote='',header=TRUE)
pops <- c('afr', 'amr', 'eas', 'fin', 'nfe', 'sas')
color_amr = k_amr = '#ED1E24'
color_eur = k_eur = '#6AA5CD'
color_afr = k_afr = '#941494'
color_sas = k_sas = '#FF9912'
color_eas = k_eas = '#108C44'
color_oth = k_oth = '#ABB9B9'
color_mde = k_mde = '#000080'
color_nfe = k_nfe = color_eur
color_fin = k_fin = '#002F6C'

pop_colors = c('afr' = color_afr,'amr' = color_amr,'eas' = color_eas,'nfe' = color_nfe,'sas' = color_sas)
pop_names = c('afr' = 'African','amr' = 'Latino','eas' = 'East Asian','nfe' = 'European','sas' = 'South Asian')

ancestries = c(names(pop_names)[c(1,2,3,4,5)])
stats$acol = pop_colors[stats$ancestry]

# Compute Odds Ratios and Confidence Intervals from ExAC and Internal HCM data
#install.packages("tidyr")
options(stringsAsFactors = FALSE)
exac_data<-read.table("../data-raw/ExACVariants_0.001_sarcomeric.txt", header=TRUE, sep="\t") %>% 
  tbl_df %>%
  transmute(Gene,CodingVar,AC=ExACAlleleCount,exacAF=ExACFreq,group="control")
hcm_data<-read.table("../data-raw/HCM_BRU_Variants_and_freqs.txt", header=TRUE, sep="\t") %>% 
  tbl_df %>%
  transmute(Gene,CodingVar,AC=CaseCount,exacAF=ExACFreq,group="case")

data <- rbind(exac_data,hcm_data)
data$bin <- cut(data$exacAF,breaks=c(1,0.001,0.0005,0.0001,0.00005,0.00001,-1),labels=c("0.00001","0.00005","0.0001","0.0005","0.001","1"))

assignGeneGroup <- function(x){
  if(x=="MYH7"){return(x)} else
    if(x=="MYBPC3"){return(x)} else
    {return("OTHER")}
}
applyGeneGroup <- function(x){
  sapply(x,assignGeneGroup)
}

assignPopSize <- function(x,y){
  if(y=="case"){return(322)} else
    if(x=="MYH7"){return(60469)} else
      if(x=="MYBPC3"){return(45793.5)} else
      {return(57855)}
}
applyPopSize<-function(x,y){
  mapply(assignPopSize,x,y)
}

ef_data<-data %>%
  mutate(Gene=applyGeneGroup(Gene)) %>%
  group_by(Gene,bin,group) %>%
  summarise(totAC=sum(AC)) %>%
  mutate(size=applyPopSize(Gene,group)) %>%
  transmute(group,odds=totAC/(size-totAC),seprecurse=((1/totAC)+(1/(size-totAC)))) %>%
  unite(oddsSEpre,odds,seprecurse,sep="_") %>%
  spread(key=group,value=oddsSEpre,fill="0_0") %>%
  separate(case,c("CaseOdds","CaseSEPre"),sep="_") %>%
  separate(control,c("ControlOdds","ControlSEPre"),sep="_") %>%
  transmute(OR=(as.numeric(CaseOdds)/as.numeric(ControlOdds)),SE=(sqrt(as.numeric(CaseSEPre)+as.numeric(ControlSEPre)))) %>%
  transmute(OR,LowerCI=(exp(log(OR)-(1.96*SE))),UpperCI=(exp(log(OR)+(1.96*SE)))) %>%
  mutate(ExACfreq=as.numeric(bin)) 

gene_list <- c('MYBPC3','MYH7','OTHER')
gene_names = c('MYBPC3' = 'MYBPC3','MYH7' = 'MYH7','OTHER' = 'Other Sarcomeric')
gene_colors = c('MYH7' = color_afr,'MYBPC3' = color_nfe,'OTHER' = color_sas)
genes = c(names(gene_names)[c(1,2,3)])
ef_data$acol = gene_colors[ef_data$Gene]

# Plot multipanel plot
png('../figures/Figure3.png',width=3800,height=1500,res=150)
m<-rbind(c(1,2))
layout(m)

par(mar=c(5,6,3,2))
plot(NA, NA, xlim=c(-2,-6), ylim=c(0,550), ann=FALSE, axes=FALSE, yaxs='i')
#abline(h=(0:5)*100, lwd=.5)
abline(h=0)
abline(v=-1.84)
for (ancestry in ancestries) {
  rows = stats$ancestry==ancestry & stats$af_limit < .05
  points(x=log10(stats$af_limit[rows]), y=stats$filtered_vars[rows], type='l', lwd=5, col=stats$acol[rows])
}
axis(side=1, at=-2:-6, labels=c("1","0.1","0.01","0.001","0.0001"), lwd=0, lwd.ticks=1, cex.axis=2)
axis(side=2, at=(0:5)*100, labels=(0:5)*100, lwd=0, lwd.ticks=0, las=2, cex.axis=2)
mtext(side=1, line=3.5, text='Maximum credible population AF (%)', cex=2.5)
mtext(side=2, line=4.5, text='Mean protein-altering variants per exome', cex=2.5)
legend('topright',pop_names[ancestries],col=pop_colors[ancestries],text.font=2,text.col=pop_colors[ancestries],lwd=5, cex=2)
mtext("(a)", side=3, adj=-0.115, line=0.5, cex=2.8)

par(mar=c(5,6,3,2))
plot(NA, NA, xlim=c(5,1), ylim=c(0,62), ann=FALSE, axes=FALSE, yaxs='i', log='x')
#abline(h=(0:6)*10, lwd=.5)
abline(h=0)
abline(v=5.33)
for (gene in genes) {
  rows = ef_data$Gene==gene
  points(ef_data$ExACfreq[rows], ef_data$OR[rows], type='b', lwd=5, pch=19, col=ef_data$acol[rows])
  arrows(ef_data$ExACfreq[rows], ef_data$LowerCI[rows], ef_data$ExACfreq[rows], ef_data$UpperCI[rows], length=0.02, angle=90, code=3, col=ef_data$acol[rows])
}
axis(side=1, at=c(5,4,3,2,1), labels=c("0.1","0.05","0.01","0.005","0.001"), lwd=0, lwd.ticks=1, cex.axis=2)
axis(side=2, at=(0:6)*10, labels=(0:6)*10, lwd=0, lwd.ticks=0, las=2, cex.axis=2)
mtext(side=1, line=3.5, text='Maximum ExAC frequency (%)', cex=2.5)
mtext(side=2, line=4.0, text='Odds ratio', cex=2.5)
legend('topleft',gene_names[genes],col=gene_colors[genes],text.font=2,text.col=gene_colors[genes],lwd=5, cex=2)
mtext("(b)", side=3, adj=-0.13, line=0.5, cex=2.8)

dev.off()
