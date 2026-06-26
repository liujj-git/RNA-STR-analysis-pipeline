# =============================================================================
# concordance.R  (revision build)
# The validated Section 2.4 pipeline is UNCHANGED. This file adds, at the end,
# two analyses requested by the reviewers:
#   * concordance vs deviation from the reference copy number       [Reviewer 2]
#   * allele-specific expression: allelic dropout + allelic balance [Reviewer 1, pt 4]
# The only edits to the original body are: (i) extraction of the PDP/AB/DAB
# FORMAT fields near the top (kept in memory for the balance analysis), and
# (ii) one '<<-' -> '<-' typo fix. All previously validated numbers reproduce.
# =============================================================================
rm(list=ls())
library(readxl)
library(vcfR)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(EnvStats)
library(data.table)
library(ggsci)
library(ggvenn)
####input
##config
config<-read.table("../GRCh38.hipstr_reference.refine.bed")
##sample_info
repeated_samples<-as.data.frame(read_xlsx("../VBsampleinfo.xlsx",sheet = "repeatability"))
colnames(repeated_samples)<-repeated_samples[1,]
repeated_samples<-repeated_samples[-1,]
sample_info<-as.data.frame(read_xlsx("../VBsampleinfo.xlsx",sheet = "concordance"))
colnames(sample_info)<-sample_info[1,]
sample_info<-sample_info[-1,]
##rawvcf,GT,DP
raw_vcf<-read.vcfR("../pstr.recode.vcf.gz")
GT<-extract.gt(raw_vcf,element = "GT")
DP<-extract.gt(raw_vcf,element = "DP")

## --- ASE per-allele FORMAT fields (kept in memory for the allelic-balance section at the end) ---
.ase_have <- all(c("PDP","AB","DAB") %in% strsplit(raw_vcf@gt[1, "FORMAT"], ":")[[1]])
if (.ase_have) {
  PDP_mat <- extract.gt(raw_vcf, element = "PDP")
  AB_mat  <- extract.gt(raw_vcf, element = "AB")
  DAB_mat <- extract.gt(raw_vcf, element = "DAB")
} else warning("VCF lacks PDP/AB/DAB; the allelic-balance section will be skipped (point the input at the raw HipSTR VCF to enable it).")
DP<-as.data.frame(DP)
DP<-data.frame("locus"=rownames(DP),gather(DP,key = "sample",value = "DP"))
DP<-na.omit(DP)
DP$DP <- as.numeric(DP$DP)
##pSTR group (annotation)
pSTR_group<-read.table("../2_annotation/STR_annotation_group.txt",header = T)


####allele sequence
##processing
#allele sequence
allele_ladder<-data.frame("locus"=raw_vcf@fix[,3],
                          "Ref"=raw_vcf@fix[,4],
                          "Alt"=raw_vcf@fix[,5])
allele_ladder<-na.omit(gather(allele_ladder,key = "allele_group",value = "sequence",-locus)) #Ref & Alt
temp1<-allele_ladder[!grepl(",",allele_ladder$sequence),]
temp2<-allele_ladder[grepl(",",allele_ladder$sequence),] #多种Alt allele
#allele group
temp3<-do.call(rbind,apply(temp2, 1, function(x){
  alt_allele<-unlist(strsplit(x[3],","))
  data.frame("locus"=x[1],
             "allele_group"=1:length(alt_allele),
             "sequence"=alt_allele,
             row.names = NULL)
}))  #多种Alt allele拆分后给予次序label （0，1，2，3...）
temp3<-temp3[nchar(temp3$sequence)!=0,] #个别alt存在错误(手动删除Human_STR_904666)
temp1<-rbind(temp1,temp3)
temp1$allele_group[temp1$allele_group=="Ref"]<-0
temp1$allele_group[temp1$allele_group=="Alt"]<-1
#ncopy label
temp1$period_size<-as.numeric(config[match(temp1$locus,config[,6]),4]) #config中提取period_size
temp2<-nchar(temp1[,3])%/%temp1[,4]
temp3<-nchar(temp1[,3])%%temp1[,4]
temp3[temp3!=0]<-paste(".",temp3[temp3!=0],sep = "")
temp3[temp3==0]<-""
temp1$ncopy<-paste(temp2,temp3,sep = "") #ncopy
temp1$ncopy_label<-paste(temp1$locus,temp1$ncopy,sep = "@") #整合基因座及ncopy生成ncopy label，用于区分同长不同序列的allele
temp2<-as.data.frame(table(temp1$ncopy_label))
temp2[,1]<-as.character(temp2[,1])
temp2<-temp2[temp2$Freq!=1,1] #提取同长不同序列的allele的ncopy label
setDT(temp1) #data.table
temp1[!temp1$ncopy_label%in%temp2,ncopy_label := ncopy] #非目标ncopy label替换回ncopy
temp1[temp1$ncopy_label%in%temp2,ncopy_label :={
  paste(ncopy, seq_along(ncopy), sep = "_")
},by=ncopy_label] #目标ncopy label以ncopy及次序进行组合并替换
temp1$unique_label<-paste(temp1$locus,temp1$allele_group,sep = "@") #整合基因座及次序label生成unique label
allele_ladder<-temp1
#output:allele_sequence.concordance.txt
allele_ladder<-allele_ladder[,c("locus","period_size","allele_group","sequence","ncopy","ncopy_label")]
allele_ladder<-allele_ladder[order(allele_ladder$locus),]
write.table(allele_ladder,file="allele_sequence.concordance.txt",col.names = T,row.names = F,quote = F,sep = "\t")
allele_sequence<-allele_ladder
allele_sequence$ncopy<-as.numeric(allele_sequence$ncopy)
rm(raw_vcf,allele_ladder)

####RNA重复样本一致性情况
##分型情况
#每个样本的分型位点数目
sample_GTcounts<-data.frame()
for (i in seq_along(repeated_samples$individual)) {
  temp1<-data.frame("individual"=repeated_samples$individual[i],
                    "group"="duplication1",
                    "sample"=repeated_samples$R1[i],
                    "GTcounts"=length(na.omit(GT[,colnames(GT)==repeated_samples$R1[i]])))
  temp2<-data.frame("individual"=repeated_samples$individual[i],
                    "group"="duplication2",
                    "sample"=repeated_samples$R2[i],
                    "GTcounts"=length(na.omit(GT[,colnames(GT)==repeated_samples$R2[i]])))
  sample_GTcounts<-rbind(sample_GTcounts,temp1,temp2)
}
sample_GTcounts<-sample_GTcounts[order(sample_GTcounts$GTcounts,decreasing = T),]
sample_GTcounts$individual<-factor(sample_GTcounts$individual,levels = unique(sample_GTcounts$individual))
sample_GTcounts$sample<-factor(sample_GTcounts$sample,levels =sample_GTcounts$sample)
sample_GTcounts$group<-factor(sample_GTcounts$group,levels = c("duplication1","duplication2"))
repeatedsample_GTcounts_barplot<-ggplot(sample_GTcounts, aes(x = individual, y = GTcounts)) +
  geom_bar(aes(fill = group),
           stat = "identity", position = position_dodge(0.8),
           width = 0.7)+
  scale_color_rickandmorty()+
  scale_fill_rickandmorty()+
  theme_bw()+
  labs(title = "Genotype Counts of Repeated RNA Samples",x="Individuals",y="Genotype Counts")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10,angle = 90),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12))
repeatedsample_GTcounts_barplot
ggsave(repeatedsample_GTcounts_barplot,file="repeatedsample_GTcounts_barplot.pdf",height = 6,width = 12)
#共有位点占比
sharedlocus_counts_df<-data.frame()
for (i in repeated_samples$individual) {
  temp1<-as.character(repeated_samples$R1[repeated_samples$individual==i])
  temp2<-as.character(repeated_samples$R2[repeated_samples$individual==i])
  temp3<-GT[,match(c(temp2,temp1),colnames(GT))]
  temp4<-data.frame("individual"=i,
                    "duplication1"=temp1,
                    "duplication2"=temp2,
                    "duplication1_locus_counts"=length(na.omit(temp3[,1])),
                    "duplication2_locus_counts"=length(na.omit(temp3[,2])),
                    "shared_locus_counts"=nrow(na.omit(temp3)))
  sharedlocus_counts_df<-rbind(sharedlocus_counts_df,temp4)
}
##分型一致性情况
#concordance_repeated_df
concordance_repeated_df<-data.frame()
for (i in repeated_samples$individual) {
  temp1<-repeated_samples$R1[repeated_samples$individual==i]
  temp2<-repeated_samples$R2[repeated_samples$individual==i]
  print(paste(temp2,temp1,collapse = "--"))
  temp3<-as.data.frame(na.omit(GT[,match(c(temp2,temp1),colnames(GT))]))
  temp3<-gather(data.frame("locus"=rownames(temp3),temp3),key = "sample",value = "GT",-"locus")
  temp3<-separate(temp3,col = "GT",into=c("v1","v2"),sep = "\\|")
  temp3[as.numeric(temp3$v1)>as.numeric(temp3$v2),c("v1","v2")]<-
    temp3[as.numeric(temp3$v1)>as.numeric(temp3$v2),c("v2","v1")] #更改allele次序，否则误判"0,1"与"1,0"
  temp3<-unite(temp3,c("v1","v2"),col = "GT",sep = ",") #更改分型间隔
  temp3<-spread(temp3,key = "sample",value = "GT")
  rownames(temp3)<-temp3$locus
  temp3<-temp3[,-1]
  temp4<-temp3[temp3[,1]==temp3[,2],]
  df1<-data.frame("individual"=i,
                  "duplication1"=temp1,
                  "duplication2"=temp2,
                  "locus"=rownames(temp4),
                  "GT1"=temp4[,1],
                  "GT2"=temp4[,2],
                  "concordance"=2) #identical GT counts
  concordance_repeated_df<-rbind(concordance_repeated_df,df1)
  temp3<-temp3[temp3[,1]!=temp3[,2],]
  if(nrow(temp3)!=0){
    df2<-data.frame() #half_identical or completely_different
    for (p in 1:nrow(temp3)) {
      GT1<-sort(as.numeric(unlist(strsplit(temp3[p,1],","))))
      GT2<-sort(as.numeric(unlist(strsplit(temp3[p,2],","))))
      temp4<-ifelse(length(intersect(GT1,GT2))==1,1,#分型只有一半相同（共享1个等位基因）
                    0) #分型无交集（共享0个等位基因）
      df2<-rbind(df2,data.frame("individual"=i,
                                "duplication1"=temp1,
                                "duplication2"=temp2,
                                "locus"=rownames(temp3)[p],
                                "GT1"=paste(GT1,collapse = ","),
                                "GT2"=paste(GT2,collapse = ","),
                                "concordance"=temp4))
    }
    concordance_repeated_df<-rbind(concordance_repeated_df,df2)
  }
}
concordance_repeated_df$DP1<-as.numeric(DP$DP[match(paste(concordance_repeated_df$duplication1,
                                               concordance_repeated_df$locus,sep="@"),
                                         paste(DP$sample,DP$locus,sep="@"))])
concordance_repeated_df$DP2<-as.numeric(DP$DP[match(paste(concordance_repeated_df$duplication2,
                                               concordance_repeated_df$locus,sep="@"),
                                         paste(DP$sample,DP$locus,sep="@"))])
write.table(concordance_repeated_df,file="concordance_repeated_df.txt",
            col.names = T,row.names = F,quote = F,sep="\t")
#重复样本对之间的IBS评分组成（0,10,20,30 DP阈值）
concordance_repeated_sample_df<-data.frame()
concordance_repeated_sample_df2<-data.frame()
for (p in c(0,10,20,30)) { #设置DP阈值
  temp4<-concordance_repeated_df[concordance_repeated_df$DP1>=p&concordance_repeated_df$DP2>=p,]
  if(nrow(temp4)!=0){
    for (i in unique(temp4$individual)) {
      temp1<-as.data.frame(table(temp4[temp4$individual==i,"concordance"]))
      temp1$percentage<-temp1$Freq/sum(temp1$Freq)
      temp2<-data.frame("DPthreshold"=p,
                        "individual"=i,
                        "group"=as.character(temp1$Var1),
                        "percentage"=temp1$percentage,
                        "counts"=temp1$Freq)
      concordance_repeated_sample_df<-rbind(concordance_repeated_sample_df,temp2)
      temp3<-data.frame("DPthreshold"=p,
                        "individual"=i,
                        "counts"=sum(temp1$Freq))
      concordance_repeated_sample_df2<-rbind(concordance_repeated_sample_df2,temp3)
    }
  }
}
temp1<-concordance_repeated_sample_df2[concordance_repeated_sample_df2$DPthreshold==0,]
temp1<-temp1[order(temp1$counts,decreasing = F),]
concordance_repeated_sample_df$individual<-factor(concordance_repeated_sample_df$individual,
                                                  levels = temp1$individual)
concordance_repeated_sample_df$DPthreshold<-factor(concordance_repeated_sample_df$DPthreshold,
                                                  levels = c(0,10,20,30))
concordance_repeated_sample_df$group<-factor(concordance_repeated_sample_df$group,
                                    levels = c(0,1,2),
                                    labels = c("non_identical","partial_identical","identical"))
concordance_repeated_sample_plot<-ggplot(concordance_repeated_sample_df[concordance_repeated_sample_df$DPthreshold==0,],
                                         aes(x=individual,y=counts,fill=group))+
  geom_bar(stat = "identity",position = "stack")+
  # geom_text(data=concordance_repeated_sample_df2[concordance_repeated_sample_df2$DPthreshold==0,],
  #           aes(x=individual,y=counts,label=counts),
  #           size=3,fontface = "bold",inherit.aes = F,hjust=-0.2)+
  coord_flip()+
  scale_fill_startrek()+
  theme_bw()+
  labs(title = "The Concordance between Repeated RNA Sample Pairs",
       x="Individual",y="Loci Counts")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12))
concordance_repeated_sample_plot
ggsave(concordance_repeated_sample_plot,file="concordance_repeated_sample_plot.pdf",width = 12,height = 12)
concordance_repeated_sample_plot2<-ggplot(concordance_repeated_sample_df[concordance_repeated_sample_df$DPthreshold!=0,],
       aes(x=individual,y=counts,fill=group))+
  geom_bar(stat = "identity",position = "stack")+
  facet_wrap(~DPthreshold)+
  coord_flip()+
  scale_fill_startrek()+
  theme_bw()+
  labs(title = "The Concordance between Repeated RNA Sample Pairs under a Sequence of DP Thresholds",
       x="Individual",y="Loci Counts")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12))
concordance_repeated_sample_plot2
ggsave(concordance_repeated_sample_plot2,file="concordance_repeated_sample_plot2.pdf",width = 15,height = 6)
#IBS1位点的step差异
halfidentical_repeated_df<-concordance_repeated_df[concordance_repeated_df$concordance==1,
                                                   c("individual","duplication1","duplication2",
                                                     "locus","GT1","GT2","DP1","DP2")]
halfidentical_repeated_df<-separate(halfidentical_repeated_df,
                                    col="GT1",into = c("GT1_A1","GT1_A2"),sep = ",")
halfidentical_repeated_df<-separate(halfidentical_repeated_df,
                                    col="GT2",into = c("GT2_A1","GT2_A2"),sep = ",")
for (i in 1:nrow(halfidentical_repeated_df)) {
  print(paste(round(i/nrow(halfidentical_repeated_df),4)*100,"%",sep = ""))
  temp1<-halfidentical_repeated_df[i,c("GT1_A1","GT1_A2")]
  temp2<-halfidentical_repeated_df[i,c("GT2_A1","GT2_A2")]                                    
  temp3<-unique(intersect(temp1,temp2)) #共有allele
  halfidentical_repeated_df$dif_GT1[i]<-ifelse(length(temp1[temp1!=temp3])!=0,temp1[temp1!=temp3],temp3)
  halfidentical_repeated_df$dif_GT2[i]<-ifelse(length(temp2[temp2!=temp3])!=0,temp2[temp2!=temp3],temp3)
}
temp1<-match(paste(halfidentical_repeated_df$locus,halfidentical_repeated_df$dif_GT1,sep = "@"),
             paste(allele_sequence$locus,allele_sequence$allele_group,sep = "@"))
temp2<-match(paste(halfidentical_repeated_df$locus,halfidentical_repeated_df$dif_GT2,sep = "@"),
             paste(allele_sequence$locus,allele_sequence$allele_group,sep = "@"))
halfidentical_repeated_df$ncopy1<-as.numeric(allele_sequence$ncopy[temp1])
halfidentical_repeated_df$ncopy2<-as.numeric(allele_sequence$ncopy[temp2])
halfidentical_repeated_df$dif_step<-abs(halfidentical_repeated_df$ncopy1-halfidentical_repeated_df$ncopy2)
halfidentical_repeated_df$sequence1<-allele_sequence$sequence[temp1]
halfidentical_repeated_df$sequence2<-allele_sequence$sequence[temp2]
halfidentical_repeated_df<-na.omit(halfidentical_repeated_df) #过滤后的vcf中存在错误格式的alt allele，需要过滤
halfidentical_repeated_df2<-data.frame()
for (p in c(0,10,20,30)) {
  temp1<-halfidentical_repeated_df[halfidentical_repeated_df$DP1>=p&
                                      halfidentical_repeated_df$DP2>=p,]
  if(nrow(temp1)!=0){
    temp2<-as.data.frame(table(temp1$dif_step[temp1$dif_step%%1==0&
                                            temp1$dif_step<=10])) #统计10步以内的整数step差异的分布
    temp2<-temp2[order(temp2$Freq,decreasing = T),]
    halfidentical_repeated_df2<-rbind(halfidentical_repeated_df2,data.frame("DPthreshold"=p,temp2))
  }
}
colnames(halfidentical_repeated_df2)<-c("DPthreshold","dif_step","counts")
halfidentical_repeated_df2$dif_step<-as.numeric(as.character(halfidentical_repeated_df2$dif_step))
halfidentical_repeated_df2$DPthreshold<-factor(halfidentical_repeated_df2$DPthreshold,levels = c(0,10,20,30))
halfidentical_repeated_df2$group[halfidentical_repeated_df2$dif_step==0]<-"sequence_discordance"
halfidentical_repeated_df2$group[halfidentical_repeated_df2$dif_step>0]<-"positive_length_discordance"
halfidentical_repeated_df2$group[halfidentical_repeated_df2$dif_step<0]<-"negative_length_discordance"
stepdif_halfidentical_plot<-ggplot(halfidentical_repeated_df2[halfidentical_repeated_df2$DPthreshold==0,],
       aes(x=dif_step,y=counts,fill=group))+
  geom_bar(stat = 'identity') +
  geom_text(aes(label=counts),size=3.5,vjust=-0.2)+
  scale_fill_manual(values = c('black', '#DC1623','#2D6DB1'))+
  scale_x_continuous(limits = c(-0.5,10.5),breaks = seq(0,10,1))+
  theme_bw()+
  labs(title = "The Distribution of Step Difference among Partial-identical Genotype Pairs between RNA-RNA Sample Pairs",
       x="Step Difference",y="Genotype Pair Counts")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 12),
        axis.text.x = element_text(face = "bold", color = "black",size = 8),
        axis.title.x = element_text(face = "bold", color = "black",size = 10),
        axis.text.y = element_text(face = "bold", color = "black",size = 8),
        axis.title.y = element_text(face = "bold", color = "black",size = 10),
        legend.position = "none")
stepdif_halfidentical_plot
ggsave(stepdif_halfidentical_plot,file="stepdif_halfidentical_plot.pdf",width = 10,height = 6)
stepdif_halfidentical_plot2<-ggplot(halfidentical_repeated_df2[halfidentical_repeated_df2$DPthreshold!=0,],
                                   aes(x=dif_step,y=counts,fill=group))+
  geom_bar(stat = 'identity') +
  facet_wrap(~DPthreshold,nrow=3)+
  scale_fill_manual(values = c('black', '#DC1623','#2D6DB1'))+
  scale_x_continuous(limits = c(-0.5,10.5),breaks = seq(0,10,1))+
  theme_bw()+
  labs(title = "The Distribution of Step Difference among Partial-identical Genotype Pairs between RNA-RNA Sample Pairs",
       x="Step Difference",y="Genotype Pair Counts")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 12),
        axis.text.x = element_text(face = "bold", color = "black",size = 8),
        axis.title.x = element_text(face = "bold", color = "black",size = 10),
        axis.text.y = element_text(face = "bold", color = "black",size = 8),
        axis.title.y = element_text(face = "bold", color = "black",size = 10),
        legend.position = "none")
stepdif_halfidentical_plot2
ggsave(stepdif_halfidentical_plot2,file="stepdif_halfidentical_plot2.pdf",width = 10,height = 6)
#同长不同序列(isoallele variant,iav)在每个基因的half-identical分型中的占比
iav_locus<-as.data.frame(table(halfidentical_repeated_df$locus))
colnames(iav_locus)<-c("locus","halfidentical_genotypepair_counts")
temp1<-as.data.frame(table(halfidentical_repeated_df$locus[halfidentical_repeated_df$dif_step==0]))
iav_locus$iav_counts<-temp1$Freq[match(iav_locus$locus,temp1$Var1)]
iav_locus$iav_counts[is.na(iav_locus$iav_counts)]<-0
iav_locus$percentage<-iav_locus$iav_counts/iav_locus$halfidentical_genotypepair_counts
#iav的变异位置、变异碱基df，用于后续3个分析
iav_sequence<-halfidentical_repeated_df[halfidentical_repeated_df$dif_step==0,
                                        c("locus","ncopy1","ncopy2","sequence1","sequence2")]
df1<-data.frame()
for (i in 1:nrow(iav_sequence)) {
  print(paste(round(i/nrow(iav_sequence)*100,2),"%",sep=""))
  temp1<-data.frame("S1"=unlist(strsplit(iav_sequence$sequence1[i],"")),
                    "S2"=unlist(strsplit(iav_sequence$sequence2[i],"")))
  temp1$concordance[temp1$S1==temp1$S2]<-1
  temp1$concordance[temp1$S1!=temp1$S2]<-0
  temp2<-data.frame("locus"=iav_sequence$locus[i],
                    "allele_length"=nchar(iav_sequence$sequence1[i]),
                    "dif_location"=as.numeric(rownames(temp1[temp1$concordance==0,])),
                    "base1"=temp1$S1[temp1$concordance==0],
                    "base2"=temp1$S2[temp1$concordance==0])
  for (p in 1:nrow(temp2)) {
    temp2$changed_base_pattern[p]<-paste(sort(as.character(temp2[p,c("base1","base2")])),collapse = "_")
  }
  df1<-rbind(df1,temp2)
}
df1$dif_location_percentage<-df1$dif_location/df1$allele_length
iav_sequence<-df1
iav_sequence$changed_base_pattern<-gsub("_","|",iav_sequence$changed_base_pattern)
#(1)同一个iav中的Point Heterogeneity Polymorphism(PHP)个数（locus与allele长度相同的iav中有几种dif_location）
iav_difbase_counts<-data.frame("iav"=unique(paste(iav_sequence$locus,"_",iav_sequence$allele_length,
                               "@",iav_sequence$dif_location,sep="")))
iav_difbase_counts<-separate(iav_difbase_counts,col="iav",into = c("iav","dif_location"),sep = "@")
iav_difbase_counts<-as.data.frame(table(iav_difbase_counts$iav))
PHPcounts_histogram<-ggplot(iav_difbase_counts,aes(x=Freq))+
  geom_histogram(fill='#2D6DB1',bins = max(iav_difbase_counts$Freq),alpha=0.6,position = "identity")+
  scale_x_continuous(breaks = seq(1,max(iav_difbase_counts$Freq),1))+
  labs(y="Genotype Pair Counts",x="Base Variant Counts")+
  theme_bw()+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12))
PHPcounts_histogram
ggsave(PHPcounts_histogram,file="PHPcounts_histogram.pdf",width = 6,height = 3)
#(2)iav PHP碱基变化类型占比
iav_changed_base_pattern<-data.frame()
for (i in c(unique(pSTR_group$final_group),"Total")) {
  if(i!="Total"){
    temp1<-pSTR_group$locus[pSTR_group$final_group==i]
    temp2<-iav_sequence[iav_sequence$locus%in%temp1,]
    temp3<-as.data.frame(table(temp2$changed_base_pattern[
      paste(temp2$locus,temp2$allele_length,sep = "_")%in%as.character(iav_difbase_counts[iav_difbase_counts$Freq==1,1])
    ]))#基于上一步过滤非单PHP locus-allele
    colnames(temp3)<-c("pattern","counts")
    temp3$percentage<-temp3$counts/sum(temp3$counts)
  }else{
    temp2<-iav_sequence
    temp3<-as.data.frame(table(temp2$changed_base_pattern[
      paste(temp2$locus,temp2$allele_length,sep = "_")%in%as.character(iav_difbase_counts[iav_difbase_counts$Freq==1,1])
    ]))#基于上一步过滤非单PHP locus-allele
    colnames(temp3)<-c("pattern","counts")
    temp3$percentage<-temp3$counts/sum(temp3$counts)
  }
  iav_changed_base_pattern<-rbind(iav_changed_base_pattern,
                                  data.frame("region"=i,temp3))
}
iav_changed_base_pattern$label<-paste(iav_changed_base_pattern$counts,
                                      " (",round(iav_changed_base_pattern$percentage,4)*100,"%)",sep = "")
od<-iav_changed_base_pattern[iav_changed_base_pattern$region=="Total",]
od<-od[order(od$percentage,decreasing = T),]
iav_changed_base_pattern$pattern<-factor(iav_changed_base_pattern$pattern,
                                         levels = as.character(od$pattern))
iav_changed_base_pattern$region<-factor(iav_changed_base_pattern$region,
                                        levels = c("CDS","5_UTR","3_UTR","intron","intergenic","Total"))
iav_changed_base_pattern_pieplot<-ggplot(iav_changed_base_pattern, aes(x = "", y = percentage, fill = pattern)) +
  geom_col(color = "black") +
  geom_label(aes(label =label),size = 2, color = rep("black",nrow(iav_changed_base_pattern)),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  facet_wrap(~region,nrow=2)+
  coord_polar(theta = "y") +
  theme_void()+
  guides(fill = guide_legend(title = "Base Variant Patterns")) +
  scale_fill_rickandmorty()
iav_changed_base_pattern_pieplot
ggsave(iav_changed_base_pattern_pieplot,file="iav_changed_base_pattern_pieplot.pdf",width = 15,height = 12)
#(3)iav PHP基变化类型在不同位置下的变化，另附iav不同位置下PHP的占比
iav_changed_base_pattern2<-iav_sequence[paste(iav_sequence$locus,iav_sequence$allele_length,sep = "_")%in%
                                          as.character(iav_difbase_counts[iav_difbase_counts$Freq==1,1]),
                                        c("locus","changed_base_pattern","dif_location_percentage")] #同样过滤非单PHP locus-allele
temp1<-as.data.frame(table(iav_changed_base_pattern2$changed_base_pattern))
temp1<-temp1[order(temp1$Freq,decreasing = T),]
iav_changed_base_pattern2$changed_base_pattern<-factor(iav_changed_base_pattern2$changed_base_pattern,
                                                       levels = as.character(temp1$Var1))
colnames(iav_changed_base_pattern2)[2]<-"Base_Variant_patterns"
iav_changed_base_pattern_regions<-ggplot(iav_changed_base_pattern2,
                                         aes(dif_location_percentage,after_stat(count),
                                             color=Base_Variant_patterns))+ 
  geom_density(position='identity', linewidth=1)+
  scale_color_rickandmorty()+
  scale_x_continuous(limits =  c(0,1))+
  labs(title = "The Base Variant Pattern across Different Allele Areas of Sequence-discordant Genotype Pairs",
       y="Counts",x="Allele Areas")+
  theme_bw()+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12))
iav_changed_base_pattern_regions
ggsave(iav_changed_base_pattern_regions,file="iav_changed_base_pattern_regions.pdf",width = 8,height = 5)


####RNA-DNA一致性情况(DNA分型要用DP=3过滤)
##分型情况
#每个样本的分型位点数目
sample_GTcounts2<-data.frame()
for (i in seq_along(unique(c(sample_info$RNA,sample_info$DNA)))) {
  temp1<-unique(c(sample_info$RNA,sample_info$DNA))[i]
  print(temp1)
  if(temp1%in%sample_info$RNA){
    temp2<-data.frame("individual"=temp1,
                      "group"="RNA",
                      "sample"=temp1,
                      "GTcounts"=length(na.omit(GT[,colnames(GT)==temp1])))
  }else{
    temp2<-data.frame("individual"=unique(sample_info$RNA[sample_info$DNA==temp1]),
                      "group"="DNA",
                      "sample"=temp1,
                      "GTcounts"=length(na.omit(GT[rownames(GT)%in%DP$locus[DP$sample==temp1&DP$DP>=3] #DNA DP>=3
                        ,colnames(GT)==temp1])))
  }
  sample_GTcounts2<-rbind(sample_GTcounts2,temp2)
} #individual以RNA样本名命名
sample_GTcounts2<-sample_GTcounts2[order(sample_GTcounts2$GTcounts,decreasing = T),]
sample_GTcounts2$individual<-factor(sample_GTcounts2$individual,levels = unique(sample_GTcounts2$individual))
sample_GTcounts2$sample<-factor(sample_GTcounts2$sample,levels =unique(sample_GTcounts2$sample))
sample_GTcounts2$group<-factor(sample_GTcounts2$group,levels = c("RNA","DNA"))
sample_GTcounts2_pointplot<-ggplot(sample_GTcounts2, aes(x=sample, y=GTcounts,col=group)) + 
  geom_point(stat='identity', size=3)  +
  facet_wrap(~group,scales = "free_x",nrow = 1)+
  scale_color_rickandmorty()+
  theme_bw()+
  labs(title = "Genotype Counts of Compared RNA/DNA Samples",
       x="Samples",y="Genotype Counts")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10,angle=90),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12),
        legend.position = "none")
sample_GTcounts2_pointplot
ggsave(sample_GTcounts2_pointplot,file="sample_GTcounts2_pointplot.pdf",height = 6,width = 15)
#共有位点占比
sharedlocus_counts_df2<-data.frame()
for (i in 1:nrow(sample_info)) {
  temp1<-as.character(sample_info$RNA[i])
  temp2<-as.character(sample_info$DNA[i])
  print(c(temp1,temp2))
  temp3<-GT[rownames(GT)%in%DP$locus[DP$sample==temp1&DP$DP>=3],#DNA DP>=3
            match(c(temp1,temp2),colnames(GT))]
  temp4<-data.frame("individual"=temp1, #individual以RNA样本名命名
                    "RNA"=temp1,
                    "DNA"=temp2,
                    "RNA_locus_counts"=length(na.omit(temp3[,1])),
                    "DNA_locus_counts"=length(na.omit(temp3[,2])),
                    "shared_locus_counts"=nrow(na.omit(temp3)))
  sharedlocus_counts_df2<-rbind(sharedlocus_counts_df2,temp4)
}
##分型一致性情况
#concordance_compared_df
concordance_compared_df<-data.frame()
for (i in 1:nrow(sample_info)) {
  temp1<-as.character(sample_info$RNA[i])
  temp2<-as.character(sample_info$DNA[i])
  print(paste(temp2,temp1,collapse = "--"))
  temp3<-as.data.frame(na.omit(GT[rownames(GT)%in%DP$locus[DP$sample==temp1&DP$DP>=3], #DNA DP>=3
                                  match(c(temp1,temp2),colnames(GT))]))
  temp3<-gather(data.frame("locus"=rownames(temp3),temp3),key = "sample",value = "GT",-"locus")
  temp3<-separate(temp3,col = "GT",into=c("v1","v2"),sep = "\\|")
  temp3[as.numeric(temp3$v1)>as.numeric(temp3$v2),c("v1","v2")]<-
    temp3[as.numeric(temp3$v1)>as.numeric(temp3$v2),c("v2","v1")] #更改allele次序，否则误判"0,1"与"1,0"
  temp3<-unite(temp3,c("v1","v2"),col = "GT",sep = ",") #更改分型间隔
  temp3<-spread(temp3,key = "sample",value = "GT")
  rownames(temp3)<-temp3$locus
  temp3<-temp3[,-1]
  temp4<-temp3[temp3[,1]==temp3[,2],]
  df1<-data.frame("individual"=temp1, #individual以RNA样本名命名
                  "RNA"=temp1,
                  "DNA"=temp2,
                  "locus"=rownames(temp4),
                  "GT_RNA"=temp4[,colnames(temp4)==temp1],
                  "GT_DNA"=temp4[,colnames(temp4)==temp2],
                  "concordance"=2) #identical GT counts
  concordance_compared_df<-rbind(concordance_compared_df,df1)
  temp3<-temp3[temp3[,1]!=temp3[,2],]
  if(nrow(temp3)!=0){
    df2<-data.frame() #half_identical or completely_different
    for (p in 1:nrow(temp3)) {
      GT1<-sort(as.numeric(unlist(strsplit(temp3[p,colnames(temp3)==temp1],","))))
      GT2<-sort(as.numeric(unlist(strsplit(temp3[p,colnames(temp3)==temp2],","))))
      temp4<-ifelse(length(intersect(GT1,GT2))==1,1,#分型只有一半相同（共享1个等位基因）
                    0) #分型无交集（共享0个等位基因）
      df2<-rbind(df2,data.frame("individual"=temp1, #individual以RNA样本名命名
                                "RNA"=temp1,
                                "DNA"=temp2,
                                "locus"=rownames(temp3)[p],
                                "GT_RNA"=paste(GT1,collapse = ","),
                                "GT_DNA"=paste(GT2,collapse = ","),
                                "concordance"=temp4))
    }
    concordance_compared_df<-rbind(concordance_compared_df,df2)
  }
}
concordance_compared_df$DP_RNA<-as.numeric(DP$DP[match(paste(concordance_compared_df$RNA,
                                                             concordance_compared_df$locus,sep="@"),
                                                    paste(DP$sample,DP$locus,sep="@"))])
concordance_compared_df$DP_DNA<-as.numeric(DP$DP[match(paste(concordance_compared_df$DNA,
                                                             concordance_compared_df$locus,sep="@"),
                                                    paste(DP$sample,DP$locus,sep="@"))])
write.table(concordance_compared_df,file="concordance_compared_df.txt",
            col.names = T,row.names = F,quote = F,sep="\t")
#特殊的DP_RNA分布，不得不作图进行额外解释
concordance_compared_RNADP_distribution<-ggplot(concordance_compared_df[concordance_compared_df$DP_RNA<100,], aes(x=DP_RNA))+ 
  geom_density(position='identity', linewidth=1)+
  scale_x_continuous(limits =  c(0,100))+
  scale_y_continuous(limits =  c(0,0.1))+
  labs(title = "The DP distribution of RNA-pSTR genotypes within RNA-DNA",
       y="Percentage",x="DP")+
  theme_bw()+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12))
concordance_compared_RNADP_distribution
ggsave(concordance_compared_RNADP_distribution,file="concordance_compared_RNADP_distribution.pdf",width = 10,height = 6)
nrow(concordance_compared_df[concordance_compared_df$DP_RNA<=10,])/nrow(concordance_compared_df)
nrow(concordance_compared_df[concordance_compared_df$DP_RNA>=30,])/nrow(concordance_compared_df)
#RNA-DNA样本对之间的IBS评分组成（0,10,20,30 DP阈值）
concordance_compared_sample_df<-data.frame()
concordance_compared_sample_df2<-data.frame()
for (p in c(0,10,20,30)) { #设置DP阈值
  temp4<-concordance_compared_df[concordance_compared_df$DP_RNA>=p,]
  if(nrow(temp4)!=0){
    for (i in unique(temp4$individual)) {
      temp1<-as.data.frame(table(temp4[temp4$individual==i,"concordance"]))
      temp1$percentage<-temp1$Freq/sum(temp1$Freq)
      temp2<-data.frame("DPthreshold"=p,
                        "individual"=i,
                        "group"=as.character(temp1$Var1),
                        "percentage"=temp1$percentage,
                        "counts"=temp1$Freq)
      concordance_compared_sample_df<-rbind(concordance_compared_sample_df,temp2)
      temp3<-data.frame("DPthreshold"=p,
                        "individual"=i,
                        "counts"=sum(temp1$Freq))
      concordance_compared_sample_df2<-rbind(concordance_compared_sample_df2,temp3)
    }
  }
}
temp1<-concordance_compared_sample_df2[concordance_compared_sample_df2$DPthreshold==0,]
temp1<-temp1[order(temp1$counts,decreasing = F),]
concordance_compared_sample_df$individual<-factor(concordance_compared_sample_df$individual,
                                                  levels = temp1$individual)
concordance_compared_sample_df$DPthreshold<-factor(concordance_compared_sample_df$DPthreshold,
                                                   levels = c(0,10,20,30))
concordance_compared_sample_df$group<-factor(concordance_compared_sample_df$group,
                                             levels = c(0,1,2),
                                             labels = c("non_identical","partial_identical","identical"))
concordance_compared_sample_plot<-ggplot(concordance_compared_sample_df[concordance_compared_sample_df$DPthreshold==0,],
                                         aes(x=individual,y=counts,fill=group))+
  geom_bar(stat = "identity",position = "stack")+
  # geom_text(data=concordance_compared_sample_df2[concordance_compared_sample_df2$DPthreshold==0,],
  #           aes(x=individual,y=counts,label=counts),
  #           size=3,fontface = "bold",inherit.aes = F,hjust=-0.2)+
  coord_flip()+
  scale_fill_startrek()+
  theme_bw()+
  labs(title = "The Concordance bewteen RNA-DNA Sample Pairs",
       x="Individual",y="Loci Counts")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 6),
        axis.title.y = element_text(face = "bold", color = "black",size = 12))
concordance_compared_sample_plot
ggsave(concordance_compared_sample_plot,file="concordance_compared_sample_plot.pdf",width = 12,height = 20)
concordance_compared_sample_plot2<-ggplot(concordance_compared_sample_df[concordance_compared_sample_df$DPthreshold!=0,],
                                          aes(x=individual,y=counts,fill=group))+
  geom_bar(stat = "identity",position = "stack")+
  facet_wrap(~DPthreshold)+
  coord_flip()+
  scale_fill_startrek()+
  theme_bw()+
  labs(title = "The Concordance between RNA-DNA Sample Pairs under a Sequence of DP Thresholds",
       x="Individual",y="Loci Counts")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 5),
        axis.title.y = element_text(face = "bold", color = "black",size = 12))
concordance_compared_sample_plot2
ggsave(concordance_compared_sample_plot2,file="concordance_compared_sample_plot2.pdf",width = 15,height = 6)
#IBS1位点的step差异
halfidentical_compared_df<-concordance_compared_df[concordance_compared_df$concordance==1,
                                                   c("individual","RNA","DNA",
                                                     "locus","GT_RNA","GT_DNA","DP_RNA","DP_DNA")]
halfidentical_compared_df<-separate(halfidentical_compared_df,
                                    col="GT_RNA",into = c("GT_RNA_A1","GT_RNA_A2"),sep = ",")
halfidentical_compared_df<-separate(halfidentical_compared_df,
                                    col="GT_DNA",into = c("GT_DNA_A1","GT_DNA_A2"),sep = ",")
for (i in 1:nrow(halfidentical_compared_df)) {
  print(paste(round(i/nrow(halfidentical_compared_df),4)*100,"%",sep = ""))
  temp1<-halfidentical_compared_df[i,c("GT_RNA_A1","GT_RNA_A2")]
  temp2<-halfidentical_compared_df[i,c("GT_DNA_A1","GT_DNA_A2")]                                    
  temp3<-unique(intersect(temp1,temp2)) #共有allele
  halfidentical_compared_df$dif_GT_RNA[i]<-ifelse(length(temp1[temp1!=temp3])!=0,temp1[temp1!=temp3],temp3)
  halfidentical_compared_df$dif_GT_DNA[i]<-ifelse(length(temp2[temp2!=temp3])!=0,temp2[temp2!=temp3],temp3)
}
temp1<-match(paste(halfidentical_compared_df$locus,halfidentical_compared_df$dif_GT_RNA,sep = "@"),
             paste(allele_sequence$locus,allele_sequence$allele_group,sep = "@"))
temp2<-match(paste(halfidentical_compared_df$locus,halfidentical_compared_df$dif_GT_DNA,sep = "@"),
             paste(allele_sequence$locus,allele_sequence$allele_group,sep = "@"))
halfidentical_compared_df$ncopy_RNA<-as.numeric(allele_sequence$ncopy[temp1])
halfidentical_compared_df$ncopy_DNA<-as.numeric(allele_sequence$ncopy[temp2])
halfidentical_compared_df$dif_step<-halfidentical_compared_df$ncopy_RNA-halfidentical_compared_df$ncopy_DNA
halfidentical_compared_df$sequence_RNA<-allele_sequence$sequence[temp1]
halfidentical_compared_df$sequence_DNA<-allele_sequence$sequence[temp2]
halfidentical_compared_df<-na.omit(halfidentical_compared_df) #过滤后的vcf中存在错误格式的alt allele，需要过滤
halfidentical_compared_df2<-data.frame()
for (p in c(0,10,20,30)) {
  temp1<-halfidentical_compared_df[halfidentical_compared_df$DP_RNA>=p,]
  if(nrow(temp1)!=0){
    temp2<-as.data.frame(table(temp1$dif_step[temp1$dif_step%%1==0&
                                                temp1$dif_step<=10&
                                                temp1$dif_step>=-10])) #统计10步以内的整数step差异的分布
    temp2<-temp2[order(temp2$Freq,decreasing = T),]
    halfidentical_compared_df2<-rbind(halfidentical_compared_df2,data.frame("DPthreshold"=p,temp2))
  }
}
halfidentical_type_df<-data.frame("type"=c("sequence_discordance","positive_length_discordance","negative_length_discordance"),
                                  "repeated"=c(length(halfidentical_repeated_df$dif_step[halfidentical_repeated_df$dif_step==0]),
                                               length(halfidentical_repeated_df$dif_step[halfidentical_repeated_df$dif_step>0]),
                                               length(halfidentical_repeated_df$dif_step[halfidentical_repeated_df$dif_step<0])),
                                  "compared"=c(length(halfidentical_compared_df$dif_step[halfidentical_compared_df$dif_step==0]),
                                               length(halfidentical_compared_df$dif_step[halfidentical_compared_df$dif_step>0]),
                                               length(halfidentical_compared_df$dif_step[halfidentical_compared_df$dif_step<0]))) #dif_step type计数
halfidentical_type_df2<-data.frame("group"=c("repeated","compared"),
                                   "half_identical_counts"=c(sum(halfidentical_type_df$repeated),sum(halfidentical_type_df$compared)),
                                   "iav_counts"=c(halfidentical_type_df$repeated[halfidentical_type_df$type=="sequence_discordance"],
                                                  halfidentical_type_df$compared[halfidentical_type_df$type=="sequence_discordance"]),
                                  "iav_percentage"=c(halfidentical_type_df$repeated[halfidentical_type_df$type=="sequence_discordance"]/sum(halfidentical_type_df$repeated),
                                                     halfidentical_type_df$compared[halfidentical_type_df$type=="sequence_discordance"]/sum(halfidentical_type_df$compared)),
                                  "non_iav_counts"=c(sum(halfidentical_type_df$repeated[halfidentical_type_df$type!="sequence_discordance"]),
                                                     sum(halfidentical_type_df$compared[halfidentical_type_df$type!="sequence_discordance"])),
                                  "one_copy_dif_counts"=c(nrow(halfidentical_repeated_df[abs(halfidentical_repeated_df$dif_step)==1,]),
                                                          nrow(halfidentical_compared_df[abs(halfidentical_compared_df$dif_step)==1,])),
                                  "one_copy_dif_percentage"=c(nrow(halfidentical_repeated_df[abs(halfidentical_repeated_df$dif_step)==1,])/
                                                                sum(halfidentical_type_df$repeated[halfidentical_type_df$type!="sequence_discordance"]),
                                                              nrow(halfidentical_compared_df[abs(halfidentical_compared_df$dif_step)==1,])/
                                                                sum(halfidentical_type_df$compared[halfidentical_type_df$type!="sequence_discordance"]))) #dif_step type计数
colnames(halfidentical_compared_df2)<-c("DPthreshold","dif_step","counts")
halfidentical_compared_df2$dif_step<-as.numeric(as.character(halfidentical_compared_df2$dif_step))
halfidentical_compared_df2$DPthreshold<-factor(halfidentical_compared_df2$DPthreshold,levels = c(0,10,20,30))
halfidentical_compared_df2$group[halfidentical_compared_df2$dif_step==0]<-"sequence_discordance"
halfidentical_compared_df2$group[halfidentical_compared_df2$dif_step>0]<-"positive_length_discordance"
halfidentical_compared_df2$group[halfidentical_compared_df2$dif_step<0]<-"negative_length_discordance"
halfidentical_compared_df2$group<-factor(halfidentical_compared_df2$group,
                                         levels = c("sequence_discordance","negative_length_discordance","positive_length_discordance"))
stepdif_halfidentical_compared_plot<-ggplot(halfidentical_compared_df2[halfidentical_compared_df2$DPthreshold==0,],
                                   aes(x=dif_step,y=counts,fill=group))+
  geom_bar(stat = 'identity',color="black", alpha=0.3) +
  geom_point()+
  geom_line(linewidth = 0.8)+
  geom_text(aes(label=counts),size=3.5,vjust=-0.2)+
  scale_fill_manual(values = c('gray','gray', 'gray'))+
  scale_x_continuous(limits = c(-10.5,10.5),breaks = seq(-10,10,1))+
  theme_bw()+
  labs(title = "The Distribution of Step Difference among Partial-identical Genotype Pairs between RNA-DNA Sample Pairs",
       x="Step Difference",y="Genotype Pair Counts")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 12),
        axis.text.x = element_text(face = "bold", color = "black",size = 8),
        axis.title.x = element_text(face = "bold", color = "black",size = 10),
        axis.text.y = element_text(face = "bold", color = "black",size = 8),
        axis.title.y = element_text(face = "bold", color = "black",size = 10),
        legend.position = "none")
stepdif_halfidentical_compared_plot
ggsave(stepdif_halfidentical_compared_plot,file="stepdif_halfidentical_compared_plot.pdf",width = 14,height =12)
stepdif_halfidentical_compared_plot2<-ggplot(halfidentical_compared_df2[halfidentical_compared_df2$DPthreshold!=0,],
                                    aes(x=dif_step,y=counts,fill=group))+
  geom_bar(stat = 'identity',color="black", alpha=0.3) +
  geom_point()+
  geom_line(linewidth = 0.8)+
  facet_wrap(~DPthreshold,nrow=3)+
  scale_fill_manual(values = c('gray','gray', 'gray'))+
  scale_x_continuous(limits = c(-10.5,10.5),breaks = seq(-10,10,1))+
  theme_bw()+
  labs(title = "The Distribution of Step Difference among Partial-identical Genotype Pairs between RNA-DNA Sample Pairs",
       x="Step Difference",y="Genotype Pair Counts")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 12),
        axis.text.x = element_text(face = "bold", color = "black",size = 8),
        axis.title.x = element_text(face = "bold", color = "black",size = 10),
        axis.text.y = element_text(face = "bold", color = "black",size = 8),
        axis.title.y = element_text(face = "bold", color = "black",size = 10),
        legend.position = "none")
stepdif_halfidentical_compared_plot2
ggsave(stepdif_halfidentical_compared_plot2,file="stepdif_halfidentical_compared_plot2.pdf",width = 15,height = 6)
#同长不同序列的half-identical(isoallele variant,iav)在每个基因座分型中的占比
#要求至少有3个half-identical分型对，以及至少有1个iav
iav_locus2<-as.data.frame(table(halfidentical_compared_df$locus))
colnames(iav_locus2)<-c("locus","halfidentical_genotypepair_counts")
temp1<-as.data.frame(table(halfidentical_compared_df$locus[halfidentical_compared_df$dif_step==0]))
iav_locus2$iav_counts<-temp1$Freq[match(iav_locus2$locus,temp1$Var1)]
iav_locus2$iav_counts[is.na(iav_locus2$iav_counts)]<-0
iav_locus2$percentage<-iav_locus2$iav_counts/iav_locus2$halfidentical_genotypepair_counts
temp1<-rbind(data.frame("group"="RNA-RNA","iav_percentage"=iav_locus$percentage[
  iav_locus$halfidentical_genotypepair_counts>=3]), #至少3个half-identical
             data.frame("group"="RNA-DNA","iav_percentage"=iav_locus2$percentage[
               iav_locus2$halfidentical_genotypepair_counts>=3])) #至少3个half-identical
temp1$group<-factor(temp1$group,levels = c("RNA-RNA","RNA-DNA"))
my_comparison<-list(c("RNA-RNA","RNA-DNA"))
iav_percentage_locus_violin<-ggplot(temp1[temp1$iav_percentage!=0,], 
                                    aes(x = group, y = iav_percentage, fill=group)) +  #排除无iav基因座的干扰
  geom_boxplot(position=position_dodge(0.9))+
  geom_jitter(aes(fill=group),width =0.05,shape = 21,size=2)+
  stat_compare_means(comparisons = my_comparison)+
  scale_fill_lancet()+
  stat_n_text(vjust = 2)+
  theme_bw()+
  labs(x="group",y="The Sequence-discordance Percentage of Partial-identical Genotype Pairs")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 10),
        axis.text.x = element_text(face = "bold", color = "black",size = 8),
        axis.title.x = element_text(face = "bold", color = "black",size = 10),
        axis.text.y = element_text(face = "bold", color = "black",size = 8),
        axis.title.y = element_text(face = "bold", color = "black",size = 8),
        legend.position = "none") #RNA-RNA分型中iav的占比与DNA-RNA分型中iav的占比的比较
iav_percentage_locus_violin
ggsave(iav_percentage_locus_violin,file="iav_percentage_locus_violin.pdf",width = 3,height = 6)
#iav的变异位置、变异碱基df，用于后续3个分析
iav_sequence2<-halfidentical_compared_df[halfidentical_compared_df$dif_step==0,
                                        c("locus","ncopy_RNA","ncopy_DNA","sequence_RNA","sequence_DNA")] #要求过滤DP<10的GT后再统计
df1<-data.frame()
for (i in 1:nrow(iav_sequence2)) {
  print(paste(round(i/nrow(iav_sequence2)*100,2),"%",sep=""))
  temp1<-data.frame("S1"=unlist(strsplit(iav_sequence2$sequence_RNA[i],"")),
                    "S2"=unlist(strsplit(iav_sequence2$sequence_DNA[i],"")))
  temp1$concordance[temp1$S1==temp1$S2]<-1
  temp1$concordance[temp1$S1!=temp1$S2]<-0
  temp2<-data.frame("locus"=iav_sequence2$locus[i],
                    "allele_length"=nchar(iav_sequence2$sequence_RNA[i]),
                    "dif_location"=as.numeric(rownames(temp1[temp1$concordance==0,])),
                    "base1"=temp1$S1[temp1$concordance==0],
                    "base2"=temp1$S2[temp1$concordance==0])
  temp2$changed_base_pattern<-paste(temp2$base2,temp2$base1,sep = "_") #RNA-DNA PHP需指定DNA to RNA的变异方向
  df1<-rbind(df1,temp2)
}
df1$dif_location_percentage<-df1$dif_location/df1$allele_length
iav_sequence2<-df1
iav_sequence2$changed_base_pattern<-gsub("_","to",iav_sequence2$changed_base_pattern)
#(1)同一个iav中的Point Heterogeneity Polymorphism(PHP)个数（locus与allele长度相同的iav中有几种dif_location）
iav_difbase_counts2<-data.frame("iav"=unique(paste(iav_sequence2$locus,"_",iav_sequence2$allele_length,
                                                  "@",iav_sequence2$dif_location,sep="")))
iav_difbase_counts2<-separate(iav_difbase_counts2,col="iav",into = c("iav","dif_location"),sep = "@")
iav_difbase_counts2<-as.data.frame(table(iav_difbase_counts2$iav))
PHPcounts_histogram2<-ggplot(iav_difbase_counts2,aes(x=Freq))+
  geom_histogram(fill='#2D6DB1',bins = max(iav_difbase_counts2$Freq),alpha=0.6,position = "identity")+
  scale_x_continuous(breaks = seq(1,max(iav_difbase_counts2$Freq),1))+
  labs(y="Genotype Pair Counts",x="Base Variant Counts")+
  theme_bw()+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12))
PHPcounts_histogram2
ggsave(PHPcounts_histogram2,file="PHPcounts_histogram2.pdf",width = 6,height = 3)
#(2)iav PHP碱基变化类型占比
iav_changed_base_pattern_compared<-data.frame()
for (i in c(unique(pSTR_group$final_group),"Total")) {
  if(i!="Total"){
    temp1<-pSTR_group$locus[pSTR_group$final_group==i]
    temp2<-iav_sequence2[iav_sequence2$locus%in%temp1,]
    temp3<-as.data.frame(table(temp2$changed_base_pattern[
      paste(temp2$locus,temp2$allele_length,sep = "_")%in%as.character(iav_difbase_counts2[iav_difbase_counts2$Freq==1,1])
    ]))#基于上一步过滤非单PHP locus-allele
    colnames(temp3)<-c("pattern","counts")
    temp3$percentage<-temp3$counts/sum(temp3$counts)
  }else{
    temp2<-iav_sequence2
    temp3<-as.data.frame(table(temp2$changed_base_pattern[
      paste(temp2$locus,temp2$allele_length,sep = "_")%in%as.character(iav_difbase_counts2[iav_difbase_counts2$Freq==1,1])
    ]))#基于上一步过滤非单PHP locus-allele
    colnames(temp3)<-c("pattern","counts")
    temp3$percentage<-temp3$counts/sum(temp3$counts)
  }
  iav_changed_base_pattern_compared<-rbind(iav_changed_base_pattern_compared,
                                  data.frame("region"=i,temp3))
}
iav_changed_base_pattern_compared$label<-paste(iav_changed_base_pattern_compared$counts,
                                      " (",round(iav_changed_base_pattern_compared$percentage,4)*100,"%)",sep = "")
od<-iav_changed_base_pattern_compared[iav_changed_base_pattern_compared$region=="Total",]
od<-od[order(od$percentage,decreasing = T),]
iav_changed_base_pattern_compared$pattern<-factor(iav_changed_base_pattern_compared$pattern,
                                         levels = as.character(od$pattern))
iav_changed_base_pattern_compared$region<-factor(iav_changed_base_pattern_compared$region,
                                        levels = c("CDS","5_UTR","3_UTR","intron","intergenic","Total"))
iav_changed_base_pattern_pieplot2<-ggplot(iav_changed_base_pattern_compared,
                                                  aes(x = "", y = percentage, fill = pattern)) +
  geom_col(color = "black") +
  geom_label(aes(label =label),size = 2, color = rep("black",nrow(iav_changed_base_pattern_compared)),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  facet_wrap(~region,nrow=2)+
  coord_polar(theta = "y") +
  theme_void()+
  guides(fill = guide_legend(title = "Base Variant Patterns")) +
  scale_fill_rickandmorty()
iav_changed_base_pattern_pieplot2
ggsave(iav_changed_base_pattern_pieplot2,file="iav_changed_base_pattern_pieplot2.pdf",width = 18,height = 15)
#(3)iav PHP基变化类型在不同位置（百分位校正）下的变化，另附iav不同位置下PHP的占比
iav_changed_base_compared_pattern2<-iav_sequence2[paste(iav_sequence2$locus,iav_sequence2$allele_length,sep = "_")%in%
                                          as.character(iav_difbase_counts2[iav_difbase_counts2$Freq==1,1]),
                                        c("locus","changed_base_pattern","dif_location_percentage")] #同样过滤非单PHP locus-allele
temp1<-as.data.frame(table(iav_changed_base_compared_pattern2$changed_base_pattern))
temp1<-temp1[order(temp1$Freq,decreasing = T),]
iav_changed_base_compared_pattern2$changed_base_pattern<-factor(iav_changed_base_compared_pattern2$changed_base_pattern,
                                                       levels = as.character(temp1$Var1))
colnames(iav_changed_base_compared_pattern2)[2]<-"Base_Variant_patterns"
iav_changed_base_pattern_regions_compared<-ggplot(iav_changed_base_compared_pattern2,
                                         aes(dif_location_percentage,after_stat(count),
                                             color=Base_Variant_patterns))+ 
  geom_density(position='identity', linewidth=1)+
  scale_color_rickandmorty()+
  scale_x_continuous(limits =  c(0,1))+
  labs(title = "The Base_Variant Pattern across Different Allele Areas of Sequence-discordant Genotype Pairs",
       y="Counts",x="Allele Areas")+
  theme_bw()+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12))
iav_changed_base_pattern_regions_compared
ggsave(iav_changed_base_pattern_regions_compared,file="iav_changed_base_pattern_regions_compared.pdf",width = 6,height = 6)


####DP以外的不一致原因
##不同基因分区下的差异
#RNA-RNA
concordance_repeated_df$pSTRgroup<-pSTR_group[match(concordance_repeated_df$locus,pSTR_group$locus),
                                              "final_group"]
concordance_region_repeated_df<-data.frame()
concordance_region_repeated_df2<-data.frame()
for (i in unique(concordance_repeated_df$pSTRgroup)) {
  temp1<-as.data.frame(table(concordance_repeated_df[concordance_repeated_df$pSTRgroup==i,"concordance"]))
  temp1$percentage<-temp1$Freq/sum(temp1$Freq)
  temp2<-data.frame("region"=i,
                    "group"=as.character(temp1$Var1),
                    "percentage"=temp1$percentage,
                    "counts"=temp1$Freq)
  concordance_region_repeated_df<-rbind(concordance_region_repeated_df,temp2)
  temp3<-data.frame("region"=i,
                    "counts"=sum(temp1$Freq))
  concordance_region_repeated_df2<-rbind(concordance_region_repeated_df2,temp3)
}
concordance_region_repeated_df$region<-factor(concordance_region_repeated_df$region,
                                     levels = c("CDS","5_UTR","3_UTR","intron","intergenic"))
concordance_region_repeated_df$group<-factor(concordance_region_repeated_df$group,
                                    levels = c(0,1,2),
                                    labels = c("non_identical","partial_identical","identical"))
concordance_region_repeated_plot<-ggplot(concordance_region_repeated_df,aes(x=region,y=percentage,fill=group))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(data=concordance_region_repeated_df2,
            aes(x=region,y=1,label=counts),inherit.aes = F,vjust=-0.2)+
  scale_fill_startrek()+
  theme_bw()+
  labs(title = "The Concordance bewteen RNA-RNA Sample Pairs in Different Genomic Regions",
       x="Region",y="Percentage")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12))
concordance_region_repeated_plot
ggsave(concordance_region_repeated_plot,file="concordance_region_repeated_plot.pdf",width = 8,height = 5)
#RNA-DNA
concordance_compared_df$pSTRgroup<-pSTR_group[match(concordance_compared_df$locus,pSTR_group$locus),
                                              "final_group"]
concordance_region_compared_df<-data.frame()
concordance_region_compared_df2<-data.frame()
for (i in unique(concordance_compared_df$pSTRgroup)) {
  temp1<-as.data.frame(table(concordance_compared_df[concordance_compared_df$pSTRgroup==i,"concordance"]))
  temp1$percentage<-temp1$Freq/sum(temp1$Freq)
  temp2<-data.frame("region"=i,
                    "group"=as.character(temp1$Var1),
                    "percentage"=temp1$percentage,
                    "counts"=temp1$Freq)
  concordance_region_compared_df<-rbind(concordance_region_compared_df,temp2)
  temp3<-data.frame("region"=i,
                    "counts"=sum(temp1$Freq))
  concordance_region_compared_df2<-rbind(concordance_region_compared_df2,temp3)
}
concordance_region_compared_df$region<-factor(concordance_region_compared_df$region,
                                              levels = c("CDS","5_UTR","3_UTR","intron","intergenic"))
concordance_region_compared_df$group<-factor(concordance_region_compared_df$group,
                                             levels = c(0,1,2),
                                             labels = c("non_identical","partial_identical","identical"))
concordance_region_compared_plot<-ggplot(concordance_region_compared_df,aes(x=region,y=percentage,fill=group))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(data=concordance_region_compared_df2,aes(x=region,y=1,label=counts),inherit.aes = F,vjust=-0.2)+
  scale_fill_startrek()+
  theme_bw()+
  labs(title = "The Concordance bewteen RNA-DNA Sample Pairs in Different Genomic Regions",
       x="Region",y="Percentage")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12))
concordance_region_compared_plot
ggsave(concordance_region_compared_plot,file="concordance_region_compared_plot.pdf",width = 8,height = 5)
##等位基因平均长度的影响
#RNA-RNA
concordance_repeated_df$locus_length<-config[match(concordance_repeated_df$locus,config[,6]),3]-
  config[match(concordance_repeated_df$locus,config[,6]),2]+1
concordance_allele_length_repeated<-concordance_repeated_df[,c("locus","concordance","locus_length")]
concordance_allele_length_repeated$concordance<-factor(concordance_allele_length_repeated$concordance,
                                              levels = c(0,1,2),
                                              labels = c("non_identical","partial_identical","identical"))
my_comparison<-list(c("non_identical","partial_identical"),
                    c("partial_identical","identical"),
                    c("non_identical","identical"))
concordance_allele_length_repeated_plot<-ggplot(concordance_allele_length_repeated, 
                                                aes(x = concordance, y = locus_length, fill=concordance)) + 
  geom_violin(trim=FALSE,alpha=0.5) + 
  geom_boxplot(width=0.1,alpha=0.5,position=position_dodge(0.9))+
  stat_compare_means(comparisons = my_comparison)+
  scale_fill_rickandmorty()+
  stat_n_text(vjust = 2)+
  theme_bw()+
  labs(title = "The Reference Allele Length of RNA-RNA Genotype Pairs Grouped by Concordance",
       x="Concordance",y="Reference Allele Length")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12),
        legend.position = "none")
concordance_allele_length_repeated_plot
ggsave(concordance_allele_length_repeated_plot,file="concordance_allele_length_repeated_plot.pdf",
       width = 10,height = 6)
#RNA-DNA
concordance_compared_df$locus_length<-config[match(concordance_compared_df$locus,config[,6]),3]-
  config[match(concordance_compared_df$locus,config[,6]),2]+1
concordance_allele_length_compared<-concordance_compared_df[,c("locus","concordance","locus_length")]
concordance_allele_length_compared$concordance<-factor(concordance_allele_length_compared$concordance,
                                                       levels = c(0,1,2),
                                                       labels = c("non_identical","partial_identical","identical"))
concordance_allele_length_compared_plot<-ggplot(concordance_allele_length_compared, 
                                                aes(x = concordance, y = locus_length, fill=concordance)) + 
  geom_violin(trim=FALSE,alpha=0.5) + 
  geom_boxplot(width=0.1,alpha=0.5,position=position_dodge(0.9))+
  stat_compare_means(comparisons = my_comparison)+
  scale_fill_rickandmorty()+
  stat_n_text(vjust = 2)+
  theme_bw()+
  labs(title = "The Reference Allele Length of RNA-DNA Genotype Pairs Grouped by Concordance",
       x="Concordance",y="Reference Allele Length")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12),
        legend.position = "none")
concordance_allele_length_compared_plot
ggsave(concordance_allele_length_compared_plot,file="concordance_allele_length_compared_plot.pdf",
       width = 10,height = 6)


##补
#（1）发生RDD的loci在重复样本间的交集
repeated_samples_compared<-repeated_samples[repeated_samples$R1%in%sample_info$RNA&
                                              repeated_samples$R2%in%sample_info$RNA,]
repeated_samples_compared_df<-data.frame()
p=1
for (i in 1:nrow(repeated_samples_compared)) {
  temp1<-concordance_compared_df$locus[concordance_compared_df$RNA==repeated_samples_compared$R1[i]&
                                         concordance_compared_df$concordance==1]
  temp2<-concordance_compared_df$locus[concordance_compared_df$RNA==repeated_samples_compared$R2[i]&
                                         concordance_compared_df$concordance==1]
  temp3<-list(temp1,temp2)
  names(temp3)<-c(repeated_samples_compared$R1[i],repeated_samples_compared$R2[i])
  nm<-paste("p",p,sep = "")
  assign(nm,ggvenn(temp3,fill_color = c('#2D6DB1', '#DC1623'),
                   stroke_size = 0.5,text_size = 2.5,set_name_size = 3.5))
  p=p+1
  if (length(intersect(temp1,temp2))!=0){
    repeated_samples_compared_df<-rbind(repeated_samples_compared_df,
                                        data.frame("group"="partial_identical",
                                                   "individual"=repeated_samples_compared$individual[i],
                                                   "locus"=intersect(temp1,temp2)))
  }
}
repeated_compared_venn_halfidentical<-ggarrange(p1,p2,p3,p4,p5,
                                                p6,p7,p8,p9,p10,
                                                p11,p12,p13,p14,p15,ncol = 5,nrow=3)
repeated_compared_venn_halfidentical
ggsave(repeated_compared_venn_halfidentical,file="repeated_compared_venn_halfidentical.pdf",
       width = 15,height = 5)
p=1
for (i in 1:nrow(repeated_samples_compared)) {
  temp1<-concordance_compared_df$locus[concordance_compared_df$RNA==repeated_samples_compared$R1[i]&
                                         concordance_compared_df$concordance==0]
  temp2<-concordance_compared_df$locus[concordance_compared_df$RNA==repeated_samples_compared$R2[i]&
                                         concordance_compared_df$concordance==0]
  temp3<-list(temp1,temp2)
  names(temp3)<-c(repeated_samples_compared$R1[i],repeated_samples_compared$R2[i])
  nm<-paste("p",p,sep = "")
  assign(nm,ggvenn(temp3,fill_color = c('#2D6DB1', '#DC1623'),
                   stroke_size = 0.5,text_size = 2.5,set_name_size = 3.5))
  p=p+1
  if (length(intersect(temp1,temp2))!=0) {
    repeated_samples_compared_df<-rbind(repeated_samples_compared_df,
                                        data.frame("group"="non_identical",
                                                   "individual"=repeated_samples_compared$individual[i],
                                                   "locus"=intersect(temp1,temp2)))
  }
}
repeated_compared_venn_nonidentical<-ggarrange(p1,p2,p3,p4,p5,
                                                p6,p7,p8,p9,p10,
                                                p11,p12,p13,p14,p15,ncol = 5,nrow=3)
repeated_compared_venn_nonidentical
ggsave(repeated_compared_venn_nonidentical,file="repeated_compared_venn_nonidentical.pdf",
       width = 15,height = 5)
rm(p1,p2,p3,p4,p5,
   p6,p7,p8,p9,p10,
   p11,p12,p13,p14,p15)
#(2)交集loci的出现频次折线
temp1<-as.data.frame(table(repeated_samples_compared_df$locus[repeated_samples_compared_df$group=="partial_identical"]))
temp1<-as.data.frame(table(temp1$Freq))
temp2<-as.data.frame(table(repeated_samples_compared_df$locus[repeated_samples_compared_df$group=="non_identical"]))
temp2<-as.data.frame(table(temp2$Freq))
temp3<-data.frame("Wrong_Locus_Frequency"=1:nrow(repeated_samples_compared))
temp3$partial_identical<-temp1$Freq[match(temp3$Wrong_Locus_Frequency,temp1$Var1)]
temp3$non_identical<-temp2$Freq[match(temp3$Wrong_Locus_Frequency,temp2$Var1)]
temp3[is.na(temp3)]<-0
temp3<-gather(temp3,key = "Discordant_Type",value = "Loci_Counts",-Wrong_Locus_Frequency)
temp3$Discordant_Type<-factor(temp3$Discordant_Type,levels = c("partial_identical","non_identical"))
Discordance_freq_line<-ggplot(temp3,aes(x=Wrong_Locus_Frequency, y=Loci_Counts, group=Discordant_Type, col=Discordant_Type)) +
  geom_line(size=1)+
  scale_x_continuous(breaks =  seq(1,nrow(repeated_samples_compared),1))+
  scale_color_lancet()+
  theme_bw()+
  labs(x="Discordance Frequency",y="Loci Counts")+
  theme(axis.title = element_text(face = "bold", color = "black",size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 10),
        axis.title.x = element_text(face = "bold", color = "black",size = 12),
        axis.text.y = element_text(face = "bold", color = "black",size = 10),
        axis.title.y = element_text(face = "bold", color = "black",size = 12))
Discordance_freq_line
ggsave(Discordance_freq_line,file="Discordance_freq_line.pdf",
       width = 6,height = 4)


summary(sample_GTcounts2$GTcounts[sample_GTcounts2$group=="DNA"])
summary(sample_GTcounts2$GTcounts[sample_GTcounts2$group=="RNA"])
summary(sharedlocus_counts_df2$shared_locus_counts)
summary(concordance_repeated_sample_df$percentage[as.character(concordance_repeated_sample_df$DPthreshold)=="30"&
                                                    as.character(concordance_repeated_sample_df$group)=="identical"])







# =============================================================================
# NEW ANALYSES (revision) — appended to the validated concordance pipeline
#   (1) concordance vs deviation from the reference copy number    [Reviewer 2]
#   (2) allele-specific expression: dropout + allelic balance      [Reviewer 1, pt 4]
# Uses in-memory objects: concordance_compared_df, allele_sequence, config, and
# the PDP/AB/DAB matrices extracted near the top. No file/VCF re-reading.
# =============================================================================
ase_theme <- theme_bw() +
  theme(axis.title  = element_text(face = "bold", colour = "black", size = 14),
        axis.text.x = element_text(face = "bold", colour = "black", size = 10),
        axis.text.y = element_text(face = "bold", colour = "black", size = 10))
AB_SIG <- log10(0.05)   # ~ -1.30 ; AB below this = significant allele bias (p < 0.05)

ccd <- concordance_compared_df
al  <- allele_sequence

ncopy_lu <- setNames(as.numeric(al$ncopy), paste(al$locus, al$allele_group, sep = "@"))
ref_lu   <- setNames(as.numeric(al$ncopy[al$allele_group == 0]), al$locus[al$allele_group == 0])
ase_idx     <- function(gt, k) vapply(strsplit(gt, ","), function(x) x[k], character(1))
ase_gf      <- function(mat, locus, sample) mat[cbind(match(locus, rownames(mat)), match(sample, colnames(mat)))]
ase_dpbin   <- function(x) cut(x, c(-Inf, 10, 20, 30, Inf), c("<10","10-19","20-29",">=30"), right = FALSE)
ase_copybin <- function(x) cut(round(abs(x)), c(-Inf, 0.5, 1.5, 2.5, 3.5, 4.5, Inf), c("0","1","2","3","4",">=5"))

ccd$rna_i1 <- ase_idx(ccd$GT_RNA, 1); ccd$rna_i2 <- ase_idx(ccd$GT_RNA, 2)
ccd$dna_i1 <- ase_idx(ccd$GT_DNA, 1); ccd$dna_i2 <- ase_idx(ccd$GT_DNA, 2)
ccd$rna_c1 <- ncopy_lu[paste(ccd$locus, ccd$rna_i1, sep = "@")]
ccd$rna_c2 <- ncopy_lu[paste(ccd$locus, ccd$rna_i2, sep = "@")]
ccd$dna_c1 <- ncopy_lu[paste(ccd$locus, ccd$dna_i1, sep = "@")]
ccd$dna_c2 <- ncopy_lu[paste(ccd$locus, ccd$dna_i2, sep = "@")]
ccd$ref_c  <- ref_lu[ccd$locus]
ccd <- ccd[complete.cases(ccd[, c("rna_c1","rna_c2","dna_c1","dna_c2","ref_c")]), ]

## (1) concordance vs deviation from reference  [Reviewer 2] ------------------
ccd$ref_dev  <- pmax(abs(ccd$dna_c1 - ccd$ref_c), abs(ccd$dna_c2 - ccd$ref_c))
ccd$dev_bin  <- ase_copybin(ccd$ref_dev)
ccd$conc_lab <- factor(ccd$concordance, levels = c(2, 1, 0),
                       labels = c("identical","partial_identical","non_identical"))
dev_tab <- as.data.frame(table(dev_bin = ccd$dev_bin, concordance = ccd$conc_lab))
dev_tot <- tapply(dev_tab$Freq, dev_tab$dev_bin, sum)
dev_tab$pct   <- 100 * dev_tab$Freq / dev_tot[as.character(dev_tab$dev_bin)]
dev_tab$n_bin <- dev_tot[as.character(dev_tab$dev_bin)]
write.table(dev_tab, "concordance_vs_refdeviation.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
ggsave("concordance_vs_refdeviation.pdf",
       ggplot(dev_tab, aes(dev_bin, pct, fill = concordance)) + geom_col(colour = "black") +
         scale_fill_startrek() + ase_theme +
         labs(title = "RNA-DNA concordance vs deviation from reference",
              x = "Max allele deviation from reference (copies)", y = "Percentage"),
       width = 9, height = 5)

## (2a) allelic dropout  [Reviewer 1, pt 4] ----------------------------------
# scope: DNA length-heterozygotes (two alleles, different length); dropout =
# RNA called homozygous (one allele not seen). Computed from GT only.
lh <- ccd[ccd$dna_i1 != ccd$dna_i2 & ccd$dna_c1 != ccd$dna_c2, ]
lh$rna_hom    <- lh$rna_i1 == lh$rna_i2
lh$dCopy_bin  <- ase_copybin(abs(lh$dna_c1 - lh$dna_c2))
lh$RNA_DP_bin <- ase_dpbin(lh$DP_RNA)
dropout_overall <- mean(lh$rna_hom)
drop_by_copy <- aggregate(rna_hom ~ dCopy_bin, lh, function(x) c(n = length(x), rate = mean(x)))
drop_by_copy <- data.frame(dCopy_bin = drop_by_copy$dCopy_bin,
                           n = as.integer(drop_by_copy$rna_hom[, "n"]), dropout_rate = drop_by_copy$rna_hom[, "rate"])
drop_by_dp <- aggregate(rna_hom ~ RNA_DP_bin, lh, function(x) c(n = length(x), rate = mean(x)))
drop_by_dp <- data.frame(RNA_DP_bin = drop_by_dp$RNA_DP_bin,
                         n = as.integer(drop_by_dp$rna_hom[, "n"]), dropout_rate = drop_by_dp$rna_hom[, "rate"])
write.table(drop_by_copy, "ase_dropout_by_copydiff.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(drop_by_dp,   "ase_dropout_by_depth.txt",    col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
ggsave("ase_dropout_by_copydiff.pdf",
       ggplot(drop_by_copy, aes(dCopy_bin, 100 * dropout_rate)) +
         geom_col(fill = "#DC1623", colour = "black") +
         geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", 100 * dropout_rate, n)), vjust = -0.2, size = 3) +
         ase_theme + labs(title = "Allelic dropout vs allele-length difference (DNA length-hets)",
                          x = "Copy-number difference between the two DNA alleles", y = "RNA homozygous-call rate (%)"),
       width = 8, height = 5)
ggsave("ase_dropout_by_depth.pdf",
       ggplot(drop_by_dp, aes(RNA_DP_bin, 100 * dropout_rate)) +
         geom_col(fill = "#2D6DB1", colour = "black") +
         geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", 100 * dropout_rate, n)), vjust = -0.2, size = 3) +
         ase_theme + labs(title = "Allelic dropout vs RNA depth (DNA length-hets)",
                          x = "RNA DP", y = "RNA homozygous-call rate (%)"),
       width = 8, height = 5)

## (2b) allelic balance  [Reviewer 1, pt 4] ----------------------------------
# scope: RNA het, length-het calls; balance from PDP (fractional reads per
# haplotype); AB = HipSTR log10 allele-bias p-value.
if (.ase_have) {
  bal <- ccd[ccd$rna_i1 != ccd$rna_i2 & ccd$rna_c1 != ccd$rna_c2, ]
  bal$PDP <- ase_gf(PDP_mat, bal$locus, bal$RNA)
  bal$AB  <- suppressWarnings(as.numeric(ase_gf(AB_mat,  bal$locus, bal$RNA)))
  bal$DAB <- suppressWarnings(as.numeric(ase_gf(DAB_mat, bal$locus, bal$RNA)))
  ase_minor <- function(s) {
    v <- suppressWarnings(as.numeric(strsplit(s, "\\|")[[1]]))
    if (length(v) != 2 || any(is.na(v)) || sum(v) == 0) return(NA_real_)
    min(v) / sum(v)
  }
  bal$minor_frac <- vapply(bal$PDP, ase_minor, numeric(1))
  bal <- bal[!is.na(bal$minor_frac), ]
  bal$dCopy_bin  <- ase_copybin(abs(bal$rna_c1 - bal$rna_c2))
  bal$RNA_DP_bin <- ase_dpbin(bal$DP_RNA)
  bal$AB_sig <- !is.na(bal$AB) & bal$AB < AB_SIG
  write.table(bal[, c("individual","RNA","locus","GT_RNA","DP_RNA","PDP","minor_frac","AB","DAB","dCopy_bin","RNA_DP_bin","AB_sig")],
              "ase_balance.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  bal_summary <- rbind(
    data.frame(stratum = "all", level = "all", n = nrow(bal),
               median_minor_frac = median(bal$minor_frac), pct_AB_significant = 100 * mean(bal$AB_sig)),
    do.call(rbind, lapply(levels(bal$dCopy_bin), function(L) { s <- bal[bal$dCopy_bin == L, ]; if (!nrow(s)) return(NULL)
      data.frame(stratum = "copy_diff", level = L, n = nrow(s), median_minor_frac = median(s$minor_frac), pct_AB_significant = 100 * mean(s$AB_sig)) })),
    do.call(rbind, lapply(levels(bal$RNA_DP_bin), function(L) { s <- bal[bal$RNA_DP_bin == L, ]; if (!nrow(s)) return(NULL)
      data.frame(stratum = "RNA_DP", level = L, n = nrow(s), median_minor_frac = median(s$minor_frac), pct_AB_significant = 100 * mean(s$AB_sig)) })))
  write.table(bal_summary, "ase_balance_summary.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  ggsave("ase_balance_distribution.pdf",
         ggplot(bal, aes(minor_frac)) + geom_histogram(binwidth = 0.02, fill = "#2D6DB1", colour = "white") +
           geom_vline(xintercept = 0.5, linetype = "dashed") + ase_theme +
           labs(title = "Allelic balance at RNA heterozygous calls", x = "Minor-allele read fraction (PDP)", y = "Genotype count"),
         width = 8, height = 5)
  ggsave("ase_balance_by_copydiff.pdf",
         ggplot(bal, aes(dCopy_bin, minor_frac, fill = dCopy_bin)) + geom_boxplot(outlier.size = 0.4) +
           geom_hline(yintercept = 0.5, linetype = "dashed") + scale_fill_lancet() + ase_theme +
           theme(legend.position = "none") +
           labs(title = "Allelic balance vs allele-length difference",
                x = "Copy-number difference between the two RNA alleles", y = "Minor-allele read fraction (PDP)"),
         width = 8, height = 5)
}

## (3) master revision file + verification of the existing §2.4 concordance ---
per_ind <- do.call(rbind, lapply(unique(ccd$individual), function(i) {
  s <- ccd[ccd$individual == i, ]
  data.frame(individual = i, identical = mean(s$concordance == 2),
             partial = mean(s$concordance == 1), non_identical = mean(s$concordance == 0))
}))
conc_rv <- character(0); rvadd <- function(...) conc_rv[[length(conc_rv) + 1]] <<- sprintf(...)
rvadd("# concordance.R revision values (auto-generated %s)", as.character(Sys.Date()))
rvadd("# matched RNA-DNA genotype pairs: %d  |  individuals: %d", nrow(ccd), nrow(per_ind))
rvadd("")
rvadd("## [VERIFY existing RNA-DNA concordance] -- vs identical 96.12%%, partial 3.7%%, non-identical 0.2%%")
rvadd("identical\t%.2f%%", 100 * mean(per_ind$identical))
rvadd("partial_identical\t%.2f%%", 100 * mean(per_ind$partial))
rvadd("non_identical\t%.2f%%", 100 * mean(per_ind$non_identical))
rvadd("")
rvadd("## [R2: CONCORDANCE vs REFERENCE DEVIATION] -> concordance_vs_refdeviation.txt")
for (b in levels(ccd$dev_bin)) {
  s <- ccd[ccd$dev_bin == b, ]
  if (nrow(s)) rvadd("dev=%s\tn=%d\tidentical %.1f%%\tpartial %.1f%%\tnon %.1f%%",
                     b, nrow(s), 100 * mean(s$concordance == 2), 100 * mean(s$concordance == 1), 100 * mean(s$concordance == 0))
}
rvadd("")
rvadd("## [R1.4a: ALLELIC DROPOUT] (DNA length-hets; dropout = RNA homozygous)")
rvadd("length-heterozygous DNA genotypes\t%d", nrow(lh))
rvadd("overall dropout rate\t%.2f%%", 100 * dropout_overall)
for (i in seq_len(nrow(drop_by_copy)))
  rvadd("dCopy=%s\tn=%d\tdropout %.2f%%", as.character(drop_by_copy$dCopy_bin[i]), drop_by_copy$n[i], 100 * drop_by_copy$dropout_rate[i])
for (i in seq_len(nrow(drop_by_dp)))
  rvadd("RNA_DP=%s\tn=%d\tdropout %.2f%%", as.character(drop_by_dp$RNA_DP_bin[i]), drop_by_dp$n[i], 100 * drop_by_dp$dropout_rate[i])
rvadd("")
if (.ase_have) {
  rvadd("## [R1.4b: ALLELIC BALANCE] (RNA het length-het calls) -> ase_balance_summary.txt")
  rvadd("RNA het length-het calls\t%d", nrow(bal))
  rvadd("median minor-allele fraction\t%.3f", median(bal$minor_frac))
  rvadd("significant allele bias (AB<%.2f)\t%.2f%%", AB_SIG, 100 * mean(bal$AB_sig))
} else rvadd("## [R1.4b: ALLELIC BALANCE] skipped -- PDP/AB/DAB not present in the VCF")
## (4) RNA-RNA concordance vs deviation from reference  [Reviewer 2] ----------
# Parallels analysis (1) but for the 25 technical-replicate (RNA-RNA) pairs.
# Uses the validated per-pair table concordance_repeated_df (concordance 2/1/0,
# GT1/GT2 = comma-separated allele indices, DP1/DP2).
rrd <- concordance_repeated_df
rrd$r1_i1 <- ase_idx(rrd$GT1, 1); rrd$r1_i2 <- ase_idx(rrd$GT1, 2)
rrd$r2_i1 <- ase_idx(rrd$GT2, 1); rrd$r2_i2 <- ase_idx(rrd$GT2, 2)
rrd$ref_c <- ref_lu[rrd$locus]
rrd$c_r1a <- ncopy_lu[paste(rrd$locus, rrd$r1_i1, sep = "@")]
rrd$c_r1b <- ncopy_lu[paste(rrd$locus, rrd$r1_i2, sep = "@")]
rrd$c_r2a <- ncopy_lu[paste(rrd$locus, rrd$r2_i1, sep = "@")]
rrd$c_r2b <- ncopy_lu[paste(rrd$locus, rrd$r2_i2, sep = "@")]
rrd <- rrd[complete.cases(rrd[, c("ref_c","c_r1a","c_r1b","c_r2a","c_r2b")]), ]
# deviation = max distance of any allele (either replicate) from the reference
rrd$ref_dev  <- pmax(abs(rrd$c_r1a - rrd$ref_c), abs(rrd$c_r1b - rrd$ref_c),
                     abs(rrd$c_r2a - rrd$ref_c), abs(rrd$c_r2b - rrd$ref_c))
rrd$dev_bin  <- ase_copybin(rrd$ref_dev)
rrd$conc_lab <- factor(rrd$concordance, levels = c(2, 1, 0),
                       labels = c("identical","partial_identical","non_identical"))
rr_dev_tab <- as.data.frame(table(dev_bin = rrd$dev_bin, concordance = rrd$conc_lab))
rr_dev_tot <- tapply(rr_dev_tab$Freq, rr_dev_tab$dev_bin, sum)
rr_dev_tab$pct   <- 100 * rr_dev_tab$Freq / rr_dev_tot[as.character(rr_dev_tab$dev_bin)]
rr_dev_tab$n_bin <- rr_dev_tot[as.character(rr_dev_tab$dev_bin)]
write.table(rr_dev_tab, "concordance_vs_refdeviation_RNARNA.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
ggsave("concordance_vs_refdeviation_RNARNA.pdf",
       ggplot(rr_dev_tab, aes(dev_bin, pct, fill = concordance)) + geom_col(colour = "black") +
         scale_fill_startrek() + ase_theme +
         labs(title = "RNA-RNA concordance vs deviation from reference",
              x = "Max allele deviation from reference (copies)", y = "Percentage"),
       width = 9, height = 5)

## (5) discordant vs identical: read-depth comparison  [Reviewer 2] ----------
# Test whether non+partial calls are enriched at low DP relative to identical
# calls (Reviewer 2 asked for an explicit two-group comparison, not just the
# threshold sweep). Large N -> report median DP per group and an effect size,
# not only the p-value.
dp_test <- function(df, dp, label) {
  d <- data.frame(grp = ifelse(df$concordance == 2, "identical", "discordant"), dp = dp)
  d <- d[is.finite(d$dp), ]
  w  <- suppressWarnings(wilcox.test(dp ~ grp, data = d))     # 'discordant' is the first level
  n1 <- sum(d$grp == "discordant"); n2 <- sum(d$grp == "identical")
  rbis <- 1 - 2 * unname(w$statistic) / (as.numeric(n1) * as.numeric(n2))  # rank-biserial (as.numeric avoids integer overflow)
  data.frame(comparison = label,
             n_discordant = n1, median_DP_discordant = median(d$dp[d$grp == "discordant"]),
             n_identical  = n2, median_DP_identical  = median(d$dp[d$grp == "identical"]),
             p_value = w$p.value, rank_biserial = rbis)
}
dp_disc_tab <- rbind(
  dp_test(ccd, ccd$DP_RNA,            "RNA-DNA (RNA depth)"),
  dp_test(rrd, pmin(rrd$DP1, rrd$DP2), "RNA-RNA (min replicate depth)"))
write.table(dp_disc_tab, "dp_discordant_vs_identical.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
ggsave("dp_discordant_vs_identical.pdf",
       ggplot(rbind(
         data.frame(comparison = "RNA-DNA", grp = ifelse(ccd$concordance == 2, "identical", "discordant"), dp = ccd$DP_RNA),
         data.frame(comparison = "RNA-RNA", grp = ifelse(rrd$concordance == 2, "identical", "discordant"), dp = pmin(rrd$DP1, rrd$DP2))),
         aes(grp, log2(dp + 1), fill = grp)) + geom_boxplot(outlier.size = 0.3) +
         facet_wrap(~comparison) + scale_fill_lancet() + ase_theme + theme(legend.position = "none") +
         labs(title = "Read depth: discordant vs identical genotype calls", x = NULL, y = "log2(DP + 1)"),
       width = 8, height = 5)

rvadd("")
rvadd("## [R2: RNA-RNA CONCORDANCE vs REFERENCE DEVIATION] -> concordance_vs_refdeviation_RNARNA.txt")
for (b in levels(rrd$dev_bin)) {
  s <- rrd[rrd$dev_bin == b, ]
  if (nrow(s)) rvadd("dev=%s\tn=%d\tidentical %.1f%%\tpartial %.1f%%\tnon %.1f%%",
                     b, nrow(s), 100 * mean(s$concordance == 2), 100 * mean(s$concordance == 1), 100 * mean(s$concordance == 0))
}
rvadd("")
rvadd("## [R2: DP, DISCORDANT vs IDENTICAL] -> dp_discordant_vs_identical.txt")
for (i in seq_len(nrow(dp_disc_tab)))
  rvadd("%s\tmedian DP discordant=%.1f (n=%d) vs identical=%.1f (n=%d)\tWilcoxon p=%.2e\trank-biserial=%.3f",
        dp_disc_tab$comparison[i], dp_disc_tab$median_DP_discordant[i], dp_disc_tab$n_discordant[i],
        dp_disc_tab$median_DP_identical[i], dp_disc_tab$n_identical[i],
        dp_disc_tab$p_value[i], dp_disc_tab$rank_biserial[i])
writeLines(conc_rv, "revision_values_concordance.txt")

cat("concordance.R (ASE + reference-deviation + RNA-RNA deviation + DP-comparison additions) done.\n")



cat("\n===== §2.3 NUMBERS TO CONFIRM =====\n")
# 1. RNA-RNA / RNA-DNA headline
cat("\n-- RNA-DNA per-individual mean concordance --\n")
cat(sprintf("identical %.2f%% partial %.2f%% non %.2f%%\n",
            100*mean(per_ind$identical),100*mean(per_ind$partial),100*mean(per_ind$non_identical)))
# 2. RNA-DNA deviation table
cat("\n-- RNA-DNA concordance vs deviation --\n")
for (b in levels(ccd$dev_bin)) { s<-ccd[ccd$dev_bin==b,]; if(nrow(s))
  cat(sprintf("dev=%s n=%d id=%.1f%% part=%.1f%% non=%.1f%%\n",b,nrow(s),
              100*mean(s$concordance==2),100*mean(s$concordance==1),100*mean(s$concordance==0))) }
cat(sprintf("frac of all genotypes at reference (dev=0): %.1f%%\n",
            100*mean(ccd$dev_bin==levels(ccd$dev_bin)[1])))
# 3. RNA-RNA deviation table (confirm non-monotonic)
cat("\n-- RNA-RNA concordance vs deviation --\n")
for (b in levels(rrd$dev_bin)) { s<-rrd[rrd$dev_bin==b,]; if(nrow(s))
  cat(sprintf("dev=%s n=%d id=%.1f%% part=%.1f%% non=%.1f%%\n",b,nrow(s),
              100*mean(s$concordance==2),100*mean(s$concordance==1),100*mean(s$concordance==0))) }
# 4. DP discordant vs identical + effect size
cat("\n-- DP discordant vs identical --\n"); print(dp_disc_tab)
# 5. dropout: the x=0 cell mystery + by copy/dp
cat("\n-- dropout by copy-number difference (CHECK dCopy=0 exists?) --\n"); print(drop_by_copy)
cat("\n-- dropout by RNA depth --\n"); print(drop_by_dp)
cat(sprintf("overall dropout %.2f%% (n length-het = %d)\n",100*dropout_overall,nrow(lh)))
# 6. balance
if(.ase_have){cat("\n-- allelic balance --\n")
  cat(sprintf("n=%d median_minor=%.3f sig_bias=%.2f%%\n",nrow(bal),
              median(bal$minor_frac),100*mean(bal$AB_sig)))
  print(bal_summary[bal_summary$stratum=="RNA_DP",])}

source("fig3_revision_patch.R")
source("extra_values_III.R")
source("ase_patch.R")
source("fig4_v3.R")


test<-read.table("concordance_repeated_df.txt",header = T)
