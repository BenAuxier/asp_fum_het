library(ggplot2)
library(gggenomes)
library(cowplot)
library(RColorBrewer)
library(cowplot)
library(tidyverse)
library(rtracklayer)

setwd("/Users/ben/Desktop/asp_het/mapping")
#load data and divide----
snp_data <- read.delim("genotypes as frequencies rounded.June9.txt")

#now lets look at the header of the file
head(snp_data)

#make chromosome factor for later use with color
snp_data$CHROM <- as.factor(snp_data$CHROM)

#this file calls the 20 parent(NIR) as genotype 0, and 21 parent (CNX) as 1. The X20 and X21 are offspring

shannon_dist <- c()
for (i in 1:1000){
  temp_shannon <- data.frame(d <- rep(0,nrow(snp_data)))
  for (j in c(7,9,6,5,5,4,5,7,5,2,6,2,5,3,4,5)){
    group_size <- j
    group <- sample(colnames(snp_data[5:ncol(snp_data)]),group_size)
    groupr <- snp_data[,colnames(snp_data) %in% group]
    groupr_zero <- 1-rowSums(groupr)/ncol(groupr)
    groupr_one <- rowSums(groupr)/ncol(groupr)
    temp_shannon <- cbind(temp_shannon,groupr_zero*(log1p(groupr_zero)) + groupr_one*(log1p(groupr_one)))
  }
  temp_shannon[is.na(temp_shannon)] <- 0
  shannon_dist <- c(shannon_dist,max(rowSums(temp_shannon)/16))
}
sig_shannon <- quantile(shannon_dist,c(0.99))
#now lets make the compatiblity groups identified by Ke
#73/91/124/206 X 15/164/180
group1 <- snp_data[,colnames(snp_data) %in% c("X73","X91","X124","X206","X15","X164","X180")]
group1_one <- rowSums(group1)/ncol(group1)
group1_shannon <- (1-group1_one)*(log1p(1-group1_one)) + group1_one*(log1p(group1_one))

#40/48/123/137/214 X 19/142/208/249  and 142/166 X 20nia  
group2 <- snp_data[,colnames(snp_data) %in% c("X40","X48","X123","X137","X214","X19","X142","X208","X249","X166")]
#add parent 20 to this by adding a column of zeroes
group2$p20 <- 0
group2_one <- rowSums(group2)/ncol(group2)
group2_shannon <- (1-group2_one)*(log1p(1-group2_one)) + group2_one*(log1p(group2_one))

#69/165/215 X 27/140
group3 <- snp_data[,colnames(snp_data) %in% c("X69","X165","X215","X27","X140")]
#now drop 69 since contaminated
group3 <- snp_data[,colnames(snp_data) %in% c("X165","X215","X27","X140")]
group3_one <- rowSums(group3)/ncol(group3)
group3_shannon <- (1-group3_one)*(log1p(1-group3_one)) + group3_one*(log1p(group3_one))

#57/113/188/193 X 199 
group4 <- snp_data[,colnames(snp_data) %in% c("X57","X113","X188","X193","X199")]
group4_one <- rowSums(group4)/ncol(group4)
group4_shannon <- (1-group4_one)*(log1p(1-group4_one)) + group4_one*(log1p(group4_one))

#89 X 78/80/82
group5 <- snp_data[,colnames(snp_data) %in% c("X89","X78","X80","X82")]
group5_one <- rowSums(group5)/ncol(group5)
group5_shannon <- (1-group5_one)*(log1p(1-group5_one)) + group5_one*(log1p(group5_one))

#94/218/224 X 100/156   
group6 <- snp_data[,colnames(snp_data) %in% c("X94","X218","X224","X100","X156")]
group6_one <- rowSums(group6)/ncol(group6)
group6_shannon <- (1-group6_one)*(log1p(1-group6_one)) + group6_one*(log1p(group6_one))

#39/55/70/145/190 X 114/177  
group7 <- snp_data[,colnames(snp_data) %in% c("X39","X55","X70","X145","X190","X114","X177")]
group7_one <- rowSums(group7)/ncol(group7)
group7_shannon <- (1-group7_one)*(log1p(1-group7_one)) + group7_one*(log1p(group7_one))

#25/67 X 116/175/213    
group8 <- snp_data[,colnames(snp_data) %in% c("X25","X67","X116","X175","X213")]
group8_one <- rowSums(group8)/ncol(group8)
group8_shannon <- (1-group8_one)*(log1p(1-group8_one)) + group8_one*(log1p(group8_one))

#133 X 77   
group9 <- snp_data[,colnames(snp_data) %in% c("X133","X77")]
group9_one <- rowSums(group9)/ncol(group9)
group9_shannon <- (1-group9_one)*(log1p(1-group9_one)) + group9_one*(log1p(group9_one))

#165/215 X 27/140/144/158 
group10 <- snp_data[,colnames(snp_data) %in% c("X165","X215","X27","X140","X144","X158")]
group10_one <- rowSums(group10)/ncol(group10)
group10_shannon <- (1-group10_one)*(log1p(1-group10_one)) + group10_one*(log1p(group10_one))

#39/70 X 152/177/208/114 
group11 <- snp_data[,colnames(snp_data) %in% c("X39","X70","X152","X177","X208","X114")]
group11_one <- rowSums(group11)/ncol(group11)
group11_shannon <- (1-group11_one)*(log1p(1-group11_one)) + group11_one*(log1p(group11_one))

#69 X170/182
group12 <- snp_data[,colnames(snp_data) %in% c("X69","X170","X182")]
group12 <- snp_data[,colnames(snp_data) %in% c("X170","X182")]
group12_one <- rowSums(group12)/ncol(group12)
group12_shannon <- (1-group12_one)*(log1p(1-group12_one)) + group12_one*(log1p(group12_one))

#22/110/207 X 209 212
group13 <- snp_data[,colnames(snp_data) %in% c("X22","X110","X207","X209","X212")]
group13_one <- rowSums(group13)/ncol(group13)
group13_shannon <- (1-group13_one)*(log1p(1-group13_one)) + group13_one*(log1p(group13_one))

#69/224 X 156/237 
group14 <- snp_data[,colnames(snp_data) %in% c("X69","X224","X156","X237")]
#now drop 69 since it was contaminated with mixed strains
group14 <- snp_data[,colnames(snp_data) %in% c("X224","X156","X237")]
group14_one <- rowSums(group14)/ncol(group14)
group14_shannon <- (1-group14_one)*(log1p(1-group14_one)) + group14_one*(log1p(group14_one))

#206/226 X 180/253 
group15 <- snp_data[,colnames(snp_data) %in% c("X206","X226","X180","X253")]
group15_one <- rowSums(group15)/ncol(group15)
group15_shannon <- (1-group15_one)*(log1p(1-group15_one)) + group15_one*(log1p(group15_one))

#108/118/189 X 21nia/32   
group16 <- snp_data[,colnames(snp_data) %in% c("X108","X118","X189","X32")]
#now need to add the parent 21 which has 1 in all positions, since reference is p20
group16$p21 <- 1
group16_one <- rowSums(group16)/ncol(group16)
group16_shannon <- (1-group16_one)*(log1p(1-group16_one)) + group16_one*(log1p(group16_one))

#71/108/189 X 32/154   
group17 <- snp_data[,colnames(snp_data) %in% c("X71","X108","X189","X32","X154")]
group17_one <- rowSums(group17)/ncol(group17)
group17_shannon <- (1-group17_one)*(log1p(1-group17_one)) + group17_one*(log1p(group17_one))

#removed the group3_shannon, since they were contaminated with wildtypes, and didn't help
mean_shannon <- (group1_shannon + group2_shannon + group4_shannon +
  group5_shannon + group6_shannon + group7_shannon + group8_shannon +
  group9_shannon + group10_shannon + group11_shannon + group12_shannon +
  group13_shannon + group14_shannon + group15_shannon + group16_shannon +
  group17_shannon)/16


#plotting----
plotting <- data.frame(CHROM = snp_data$CHROM,
                       POS = snp_data$POS,
                       mean_shannon = mean_shannon)

#to make widths of chromosomes represent actual value, I got these from the Afum_p20.fna.fai file
fai <- c(4650176,4865329,4022249,3154315,3912066,3799661,1733642,1784527)

top_shannon <- sort(plotting$mean_shannon)[floor(length(mean_shannon))*0.95]

labels = data.frame(CHROM=c("chr2","chr5","chr6","chr8"),
                    POS  =c(3.9,1.0,2.9,1.07),
                    mean_shannon=c(0.66,0.68,0.67,0.67),
                    label=c("italic(het)*A","italic(het)*B","italic(het)*C","italic(het)*D"))

main <- ggplot(plotting) + geom_point(aes(x=POS/1000000,y=mean_shannon,col=CHROM),size=0.5) +
  geom_hline(aes(yintercept=sig_shannon),lty=3)+
  scale_x_continuous(breaks=c(0,1,2,3,4))+
  #labs(title="Ke+Joost with some uncertain")+
  scale_color_brewer(palette="Dark2") +
  facet_grid(cols=vars(CHROM),scales = "free_x",space="free_x") +
  theme_bw() + labs(x="Position (Mb)",y="Shannon\nDivergence") +
  geom_text(data=labels,aes(x=POS,y=mean_shannon,label=label),size=3,parse=T)+
  #geom_rect(data=joost_regions,aes(xmin=(start-50000)/1000000,xmax=(end+50000)/1000000,ymin=0.68,ymax=0.72),inherit.aes = F,size=3)+
  theme(panel.grid.minor.x = element_blank(),
        legend.position="none",
        panel.spacing= unit(0.2,"lines"),
        panel.grid = element_blank(),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.title.x = element_blank(),
        panel.background = element_rect(fill=NA),
        strip.background = element_rect(fill=NA))

main

loc1 <- ggplot(plotting[plotting$CHROM=="chr2",]) + theme_minimal()+
  coord_cartesian(xlim=c(4.675-0.03,4.675+0.03),ylim=c(0.5,0.7),expand=F)+
  scale_x_continuous(breaks=seq(from=4.65,to=4.70,by=0.01))+
  geom_point(aes(x=POS/1000000,y=mean_shannon),size=1.5,alpha=0.7,color="#D95F02")+
  geom_line(aes(x=POS/1000000,y=mean_shannon),lwd=0.5,color="grey40") +
  geom_hline(yintercept=sig_shannon,lty=3)+
  geom_vline(xintercept=c(4.666709,4.687709),lty=2)+ 
  theme_minimal() +
  theme(axis.title = element_blank(),
        panel.grid.minor=element_blank())

loc2 <- ggplot(plotting[plotting$CHROM=="chr5",]) + theme_minimal()+
  coord_cartesian(xlim=c(0.198-0.03,0.198+0.03),ylim=c(0.5,0.7),expand=F)+
  geom_point(aes(x=POS/1000000,y=mean_shannon),size=1.5,alpha=0.7,color="#66A61E")+
  geom_line(aes(x=POS/1000000,y=mean_shannon),lwd=0.5,color="grey40") +
  geom_hline(yintercept=sig_shannon,lty=3)+
  geom_vline(xintercept=c(0.183767,0.203767),lty=2)+
  theme(axis.title = element_blank(),
        axis.text.y=element_blank(),
        panel.grid.minor=element_blank())

loc3 <- ggplot(plotting[plotting$CHROM=="chr6",]) + theme_minimal()+
  coord_cartesian(xlim=c(3.62-0.03,3.62+0.03),ylim=c(0.5,0.7),expand=F)+
  geom_point(aes(x=POS/1000000,y=mean_shannon),size=1.5,alpha=0.7,color="#E6AB02")+
  geom_line(aes(x=POS/1000000,y=mean_shannon),lwd=0.5,color="grey40") +
  geom_hline(yintercept=sig_shannon,lty=3)+
  geom_vline(xintercept=c(3.607012,3.627012),lty=2)+
  theme(axis.title = element_blank(),
        #axis.text.y = element_blank(),
        panel.grid.minor=element_blank())

loc4 <- ggplot(plotting[plotting$CHROM=="chr8",]) + theme_minimal()+
  coord_cartesian(xlim=c(0.39-0.03,0.39+0.03),ylim=c(0.5,0.7),expand=F)+
  geom_point(aes(x=POS/1000000,y=mean_shannon),size=1.5,alpha=0.7,color="#666666")+
  geom_line(aes(x=POS/1000000,y=mean_shannon),lwd=0.5,color="grey40") +
  geom_hline(yintercept=sig_shannon,lty=3)+
  geom_vline(xintercept=c(0.386300,0.4063),lty=2)+
  theme(axis.title=element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank())

sub <- plot_grid(loc1,loc2,loc3,loc4,nrow=2,labels=c("B) hetC","C) hetB","D) hetC","E) hetD"))

#svg(paste0("asp_het_plotting",Sys.Date(),".svg"),width=8,height=6)
plot_grid(main,sub,nrow=2,rel_heights=c(0.35,0.65),labels=c("A)",""))
#dev.off()

#now we need to make use of gggenomes package. the first step is to do minimpa all vs. all alignment for each gene, using 
#https://thackl.github.io/gggenomes/articles/emales.html#compare-genome-synteny

setwd("/Users/ben/Desktop/asp_het/fine_mapping_genes/")
hetA_seqs <- read_seq_len("hetA.both.fasta") %>% arrange(desc(seq_id))
hetA_genes <- read_gff("hetA.both.gff3")
hetA_genes <- subset(hetA_genes,type == "gene")
hetA_genes$type <- "CDS"
#expect warnings here, no problem
hetA_links <- read_paf("hetA.paf")
#we want to subset to 10k up/downstream of highest value, which is 4677709
#hetA goes from 4644999, so the peak is at position 32710 in the seqence
#I want from position 4666709 to 4687709
hetA_gene <- gggenomes(seqs=hetA_seqs,hetA_genes,links=hetA_links) + 
  scale_x_continuous(breaks=c(3284,13284),labels=c(4.67,4.68))+
  geom_link(fill="grey70",col="grey70",offset=0) +geom_seq() + geom_gene(size=3,shape=5,fill="cornflowerblue")
hetA_comb <- plot_grid(hetA_gene,loc1,nrow=2,align="h")
hetA_comb

hetB_seqs <- read_seq_len("hetB.both.fasta") %>% arrange(desc(seq_id))
hetB_genes <- read_gff("hetB.both.gff3")
hetB_genes <- subset(hetB_genes,type == "gene")
hetB_genes$type <- "CDS"
#expect warnings here, no problem
hetB_links <- read_paf("hetB.paf")
hetB_gene <- gggenomes(seqs=hetB_seqs,hetB_genes,links=hetB_links) + 
  scale_x_continuous(breaks=c(6232,16232),labels=c("0.19","0.20"))+
  geom_link(fill="grey70",col="grey70",offset=0) +geom_seq() + geom_gene(size=3,shape=5,fill="cornflowerblue")
hetB_comb <- plot_grid(hetB_gene,loc2,nrow=2,align="h")
hetB_comb


hetC_seqs <- read_seq_len("hetC.both.fasta") %>% arrange(desc(seq_id))
hetC_genes <- read_gff("hetC.both.gff3")
hetC_genes <- subset(hetC_genes,type == "gene")
hetC_genes$type <- "CDS"
#expect warnings here, no problem
hetC_links <- read_paf("hetC.paf")
hetC_gene <- gggenomes(seqs=hetC_seqs,hetC_genes,links=hetC_links) + 
  scale_x_continuous(breaks=c(2986,12986),labels=c(3.61,3.62))+
  geom_link(fill="grey70",col="grey70",offset=0) +geom_seq() + geom_gene(size=3,shape=5,fill="cornflowerblue")
hetC_comb <- plot_grid(hetC_gene,loc3,nrow=2,align="h")
hetC_comb
hetD_seqs <- read_seq_len("hetD.both.fasta") %>% arrange(desc(seq_id))
hetD_genes <- read_gff("hetD.both.gff3")
hetD_genes <- subset(hetD_genes,type == "gene")
hetD_genes$type <- "CDS"
#expect warnings here, no problem
hetD_links <- read_paf("hetD.paf")
hetD_gene <- gggenomes(seqs=hetD_seqs,hetD_genes,links=hetD_links) + 
  scale_x_continuous(breaks=c(3692,13692),labels=c("0.39","0.40"))+
  geom_link(fill="grey70",col="grey70",offset=0) +geom_seq() + geom_gene(size=3,shape=5,fill="cornflowerblue")
hetD_comb <- plot_grid(hetD_gene,loc4,nrow=2,align="h")
hetD_comb
sub <- plot_grid(hetA_comb,hetB_comb,hetC_comb,hetD_comb,nrow=2,labels=c("B) hetA","C) hetB","D) hetC","E) hetD"))

plot_grid(main,sub,nrow=2,rel_heights=c(0.35,0.65),labels=c("A)",""))
svg(paste0("asp_het_plotting",Sys.Date(),".svg"),width=10,height=8)
plot_grid(main,sub,nrow=2,rel_heights=c(0.35,0.65),labels=c("A)",""))
dev.off()
#second analysis, for slow mutating locus----
#Original group 2 #40/48/123/137/214 X 19/142/208/249  and 142/166 X 20nia  
#1a 19/40/48/14/208/249
group1a <- snp_data[,colnames(snp_data) %in% c("X19","X40","X48","X142","X208","X249")]
group1a_one <- rowSums(group1a)/ncol(group1a)
group1a_zero <- 1-group1a_one
group1a_shannon <- group1a_zero*(log1p(group1a_zero)) + group1a_one*(log1p(group1a_one))
group1a_score <- rowSums(group1a)/ncol(group1a)

#1b 123/137/214
group1b <- snp_data[,colnames(snp_data) %in% c("X123","X137","X214")]
group1b_one <- rowSums(group1b)/ncol(group1b)
group1b_zero <- 1-group1b_one
group1b_shannon <- group1b_zero*(log1p(group1b_zero)) + group1b_one*(log1p(group1b_one))
group1b_score <- rowSums(group1b)/ncol(group1b)

#3a 188/199
group3a <- snp_data[,colnames(snp_data) %in% c("X188","X199")]
group3a_one <- rowSums(group3a)/ncol(group3a)
group3a_zero <- 1-group3a_one
group3a_shannon <- group3a_zero*(log1p(group3a_zero)) + group3a_one*(log1p(group3a_one))
group3a_score <- rowSums(group3a)/ncol(group3a)

#3b 57/133/193
group3b <- snp_data[,colnames(snp_data) %in% c("X57","X133","X193")]
group3b_one <- rowSums(group3b)/ncol(group3b)
group3b_zero <- 1-group3b_one
group3b_shannon <- group3b_zero*(log1p(group3b_zero)) + group3b_one*(log1p(group3b_one))
group3b_score <- rowSums(group3b)/ncol(group3b)

#7a 110/207/212
group7a <- snp_data[,colnames(snp_data) %in% c("X110","X207","X212")]
group7a_one <- rowSums(group7a)/ncol(group7a)
group7a_zero <- 1-group7a_one
group7a_shannon <- group7a_zero*(log1p(group7a_zero)) + group7a_one*(log1p(group7a_one))
group7a_score <- rowSums(group7a)/ncol(group7a)

#7b 22/209
group7b <- snp_data[,colnames(snp_data) %in% c("X22","X209")]
group7b_one <- rowSums(group7b)/ncol(group7b)
group7b_zero <- 1-group7b_one
group7b_shannon <- group7b_zero*(log1p(group7b_zero)) + group7b_one*(log1p(group7b_one))
group7b_score <- rowSums(group7b)/ncol(group7b)

#9a 25/116/175/213
group9a <- snp_data[,colnames(snp_data) %in% c("X25","X116","X175","X213")]
group9a_one <- rowSums(group9a)/ncol(group9a)
group9a_zero <- 1-group9a_one
group9a_shannon <- group9a_zero*(log1p(group9a_zero)) + group9a_one*(log1p(group9a_one))
group9a_score <- rowSums(group9a)/ncol(group9a)

#9b 67
#skip

#14a 15/73/180/206/251
group14a <- snp_data[,colnames(snp_data) %in% c("X15","X73","X180","X206","X251")]
group14a_one <- rowSums(group14a)/ncol(group14a)
group14a_zero <- 1-group14a_one
group14a_shannon <- group14a_zero*(log1p(group14a_zero)) + group14a_one*(log1p(group14a_one))
group14a_score <- rowSums(group14a)/ncol(group14a)

#14b 91/124/164
group14b <- snp_data[,colnames(snp_data) %in% c("X91","X124","X164")]
group14b_one <- rowSums(group14b)/ncol(group14b)
group14b_zero <- 1-group14b_one
group14b_shannon <- group14b_zero*(log1p(group14b_zero)) + group14b_one*(log1p(group14b_one))
group14b_score <- rowSums(group14b)/ncol(group14b)

#16a 32/71/154
group16a <- snp_data[,colnames(snp_data) %in% c("X32","X71","X154")]
group16a_one <- rowSums(group16a)/ncol(group16a)
group16a_zero <- 1-group16a_one
group16a_shannon <- group16a_zero*(log1p(group16a_zero)) + group16a_one*(log1p(group16a_one))
group16a_score <- rowSums(group16a)/ncol(group16a)

#16b 108/118/189
group16b <- snp_data[,colnames(snp_data) %in% c("X108","X118","X189")]
group16b_one <- rowSums(group16b)/ncol(group16b)
group16b_zero <- 1-group16b_one
group16b_shannon <- group16b_zero*(log1p(group16b_zero)) + group16b_one*(log1p(group16b_one))
group16b_score <- rowSums(group16b)/ncol(group16b)


shannon_dist5 <- c()
for (i in 1:1000){
  temp_shannon5 <- data.frame(d <- rep(0,nrow(snp_data)))
  for (j in c(6,3,2,3,3,2,4,5,3,3,3)){
    group_size <- j
    group <- sample(colnames(snp_data[5:ncol(snp_data)]),group_size)
    groupr <- snp_data[,colnames(snp_data) %in% group]
    groupr_zero <- 1-rowSums(groupr)/ncol(groupr)
    groupr_one <- rowSums(groupr)/ncol(groupr)
    temp_shannon5 <- cbind(temp_shannon5,groupr_zero*(log1p(groupr_zero)) + groupr_one*(log1p(groupr_one)))
  }
  temp_shannon5[is.na(temp_shannon5)] <- 0
  shannon_dist5 <- c(shannon_dist5,max(rowSums(temp_shannon5)/11))
}

sig_shannon5 <- quantile(shannon_dist5,0.95)

differences <- (abs(group1a_score - group1b_score)+
                  abs(group3a_score - group3b_score)+
                  abs(group7a_score - group7b_score)+
                  abs(group14a_score- group14b_score)+
                  abs(group16a_score- group16b_score))/5
shannon_loc5 <- (group1a_shannon + group1b_shannon + group2_shannon + group3a_shannon + group3b_shannon +
                   group7a_shannon + group7b_shannon + group9a_shannon +
                   group14a_shannon + group14b_shannon + group16a_shannon + group16b_shannon)/11

sub_plotting <- data.frame(
                           CHRS = snp_data$CHROM,
                           POS = snp_data$POS,
                           differences = differences)
label <- data.frame(CHRS="chr6",
                    POS=2.810000,
                    differences=0.89,
                    label="italic(het)*E")
loc_5 <- ggplot(sub_plotting) + theme_minimal() +
  geom_point(aes(x=POS/1000000,y=differences,col=CHRS),size=0.5) +
  scale_color_brewer(palette="Dark2") +
  scale_x_continuous(breaks=c(0,1,2,3,4))+
  scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8))+
  geom_hline(yintercept=sig_shannon5,lty=2)+
  geom_text(data=label,aes(x=POS,y=differences,label=label),parse=T)+
  facet_grid(cols=vars(CHRS),scales = "free_x",space="free_x") +
  theme_bw() + 
  labs(x="Position (Mb)",y="Between-group\nDifference") +
  theme(panel.grid = element_blank(),
        legend.position="none")
sub_loc5 <- ggplot(sub_plotting[sub_plotting$CHRS=="chr6",]) +
  scale_x_continuous(breaks=seq(from=1.50,to=1.70,by=0.02))+
  coord_cartesian(xlim=c(1.60-0.1,1.60+0.1),ylim=c(0.3,0.95))+
  geom_vline(aes(xintercept=1.56),lty=2)+
  geom_vline(aes(xintercept=1.65),lty=2)+
  geom_line(aes(x=POS/1000000,y=differences),size=0.5,color="grey40")+
  geom_point(aes(x=POS/1000000,y=differences),size=1.5,alpha=0.7,color="#E6AB02")+
  labs(x="Position (Mb)")+ theme_minimal() +
  theme(axis.title.y = element_blank(),
        panel.grid.minor=element_blank())
hetE_seqs <- read_seq_len("hetE.both.fasta") %>% arrange(desc(seq_id))
hetE_genes <- read_gff("hetE.both.gff3")
hetE_genes <- subset(hetE_genes,type == "gene")
hetE_genes$type <- "CDS"
#expect warnings here, no problem
hetE_links <- read_paf("hetE.paf")
hetE_gene <- gggenomes(seqs=hetE_seqs,hetE_genes,links=hetE_links) + 
  geom_link(fill="grey70",col="grey70",offset=0) +
  geom_seq() + geom_gene(size=3,shape=5,fill="cornflowerblue")+
  #starts at 1560000
  scale_x_continuous(breaks=c(0,20000,40000,60000,80000),labels=c("1.56","1.58","1.60","1.62","1.64"))
#hetE_gene
sub <- plot_grid(hetE_gene,sub_loc5,nrow=2)

svg(paste0("Supplemental.hetE.",Sys.Date(),".svg"),width=8,height=6)
plot_grid(loc_5,sub,nrow=2,labels=c("A)","B)"))
dev.off()
