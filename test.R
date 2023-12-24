library(data.table);library(tibble);library(dplyr)

# Rsq
info_all<-fread('C:/Users/minni/OneDrive - Harvard University/Research/VTE/Analysis/all_chrs.info')
info_all_Rsq<-info_all %>% select(SNP,Rsq)

# Sara's GW-sig 37 SNPs
sara37<-read.table("C:/Users/minni/OneDrive - Harvard University/Research/VTE/Sarah/VTE_37SNP_PRS_Lindstrom.txt")
colnames(sara37)<-c("chr","rsid","pos","a1","a0","weight_Lindstrom")
sara37$SNP <- paste(sara37$chr, sara37$pos, sep=":")
rs<-sara37$rsid;pos<-sara37$pos

sara37_info<-inner_join(sara37,info_all_Rsq,by="SNP")

write.csv(sara37_info,"C:/Users/minni/OneDrive - Harvard University/Research/VTE/Sarah/sara37_info.csv",row.names = F)

#LD matrix
ld_eur <- readRDS("C:/Users/minni/OneDrive - Harvard University/Research/VTE/Analysis/LD_matrix/ld.Feb_13_2022")
ld_afr <- readRDS("C:/Users/minni/OneDrive - Harvard University/Research/VTE/Analysis/LD_matrix/ld_AFR.Apr_02_2022")
#corr <- readRDS("C:/Users/minni/OneDrive - Harvard University/Research/VTE/Analysis/LD_matrix/corr.Feb_13_2022")
#df_beta <- readRDS("C:/Users/minni/OneDrive - Harvard University/Research/VTE/Analysis/LD_matrix/df_beta.Feb_13_2022") 
#h2_est <- readRDS("C:/Users/minni/OneDrive - Harvard University/Research/VTE/Analysis/LD_matrix/h2_est.Feb_13_2022") # updated - imputed BIM file

ld_eur<-as.data.frame(ld_eur)
rownames(ld_eur)<- info_snp_new_EUR$rsid
colnames(ld_eur)<-"ld_eur"
ld_eur <- tibble::rownames_to_column(ld_eur, "rsid")

ld_afr<-as.data.frame(ld_afr)
rownames(ld_afr)<- info_snp_new_AFR$rsid
colnames(ld_afr)<-"ld_afr"
ld_afr <- tibble::rownames_to_column(ld_afr, "rsid")

sara37_info_ld_eur<-left_join(sara37_info,ld_eur,by="rsid")
sara37_info_ld<-left_join(sara37_info_ld_eur,ld_afr,by="rsid")

############
## LDPRED2
############
setwd("C:/Users/minni/OneDrive - Harvard University/Research/VTE/GBMI/Attachments/LDPRED2")
info_snp_new_EUR<-readRDS("info_snp_EUR")
info_snp_new_AFR<-readRDS("info_snp_AFR")

beta_auto_EUR<-readRDS("beta_auto_EUR.Jun_10_2022") #604741     48
beta_auto_AFR<-readRDS("beta_auto_AFR.Apr_03_2022") #1184805      48

rownames(beta_auto_EUR)<- info_snp_new_EUR$rsid
rownames(beta_auto_AFR)<- info_snp_new_AFR$rsid

#check weight distributions
summary(rowMeans(beta_auto_EUR));quantile(rowMeans(beta_auto_EUR), probs = c(.9))
summary(rowMeans(beta_auto_AFR));quantile(rowMeans(beta_auto_AFR), probs = c(.9))

#beta_auto_EUR[rs,]
#beta_auto_EUR["rs6025",]
#beta_auto_EUR[which(rownames(beta_auto_EUR)=="rs3131972"),]

sara37_ldpred2_eur<-beta_auto_EUR[which(rownames(beta_auto_EUR) %in% rs),]
sara37_ldpred2_afr<-beta_auto_AFR[which(rownames(beta_auto_AFR) %in% rs),]

sara37_ldpred2_eur <- tibble::rownames_to_column(as.data.frame(sara37_ldpred2_eur), "rsid")
sara37_ldpred2_afr <- tibble::rownames_to_column(as.data.frame(sara37_ldpred2_afr), "rsid")

sara37_ldpred2_eur_merged <- left_join(sara37_info_ld,sara37_ldpred2_eur,by="rsid")
sara37_ldpred2_afr_merged <- left_join(sara37_info_ld,sara37_ldpred2_afr,by="rsid")

write.csv(sara37_ldpred2_eur_merged,"C:/Users/minni/OneDrive - Harvard University/Research/VTE/Sarah/sara37_ldpred2_eur_merged.csv",row.names = F)
write.csv(sara37_ldpred2_afr_merged,"C:/Users/minni/OneDrive - Harvard University/Research/VTE/Sarah/sara37_ldpred2_afr_merged.csv",row.names = F)
#info_snp_new_EUR[info_snp_new_EUR$rsid=="rs6025",]
#info_snp_new_AFR[info_snp_new_AFR$rsid=="rs6025",]

sara37_ldpred2_eur_rs<-info_snp_new_EUR[info_snp_new_EUR$rsid %in% c(rs),] #find by rsid #10
sara37_ldpred2_eur_pos<-info_snp_new_EUR[info_snp_new_EUR$pos %in% c(pos),] #find by position

sara37_ldpred2_afr_rs<-info_snp_new_AFR[info_snp_new_AFR$rsid %in% c(rs),] #find by rsid #15
sara37_ldpred2_afr_pos<-info_snp_new_AFR[info_snp_new_AFR$pos %in% c(pos),] #find by position

############
## PRSCSX
############
setwd("C:/Users/minni/OneDrive - Harvard University/Research/VTE/GBMI/Attachments/PRSCSX")
prscsx_chr1_afr<-read.table("weight_AFR_pst_eff_a1_b0.5_phiauto_chr1.txt",header = F)
prscsx_chr1_eur<-read.table("weight_EUR_pst_eff_a1_b0.5_phiauto_chr1.txt",header = F)

df_eur<-data.frame()
for(i in 1:22) {
  data <- read.table(paste0("weight_EUR_pst_eff_a1_b0.5_phiauto_chr",i,".txt"),header = F)
  df_eur <- rbind(df_eur,data)
}
colnames(df_eur)<-c("chr","rsid","pos","a0","a1","weight_prscsx_eur")
sara37_prscsx_eur <- df_eur %>% select(rsid,weight_prscsx_eur)

df_afr<-data.frame()
for(i in 1:22) {
  data <- read.table(paste0("weight_AFR_pst_eff_a1_b0.5_phiauto_chr",i,".txt"),header = F)
  df_afr <- rbind(df_afr,data)
}
colnames(df_afr)<-c("chr","rsid","pos","a0","a1","weight_prscsx_afr")
sara37_prscsx_afr <- df_afr %>% select(rsid,weight_prscsx_afr)

#check weight distributions
summary(sara37_prscsx_eur$weight_prscsx_eur);quantile(sara37_prscsx_eur$weight_prscsx_eur, probs = c(.9))
summary(sara37_prscsx_afr$weight_prscsx_afr);quantile(sara37_prscsx_afr$weight_prscsx_afr, probs = c(.9))

#sara37_prscsx_eur_merged <- left_join(sara37_info,sara37_prscsx_eur,by="rsid")
#sara37_prscsx_afr_merged <- left_join(sara37_info,sara37_prscsx_afr,by="rsid")

sara37_prscsx_eur_merged <- left_join(sara37_info_ld,sara37_prscsx_eur,by="rsid")
sara37_prscsx_merged <- left_join(sara37_prscsx_eur_merged,sara37_prscsx_afr,by="rsid")

#write.csv(sara37_prscsx_eur_merged,"C:/Users/minni/OneDrive - Harvard University/Research/VTE/Sarah/sara37_prscsx_eur_merged.csv",row.names = F)
#write.csv(sara37_prscsx_afr_merged,"C:/Users/minni/OneDrive - Harvard University/Research/VTE/Sarah/sara37_prscsx_afr_merged.csv",row.names = F)
write.csv(sara37_prscsx_merged,"C:/Users/minni/OneDrive - Harvard University/Research/VTE/Sarah/sara37_prscsx_merged.csv",row.names = F)

#prscsx_chr1_afr[prscsx_chr1_afr$rsid=="rs6025",]
#prscsx_chr1_eur[prscsx_chr1_eur$rsid=="rs6025",]

sara37_prscsx_eur_rs<-df_eur[df_eur$rsid %in% c(rs),] #find by rsid #10
sara37_prscsx_eur_pos<-df_eur[df_eur$pos %in% c(pos),] #find by position

sara37_prscsx_afr_rs<-df_afr[df_afr$rsid %in% c(rs),] #find by rsid #10
sara37_prscsx_afr_pos<-df_afr[df_afr$pos %in% c(pos),] #find by position

############################################################################################
## PK's code - CompareLindstromJee.Rd
# PRS odds ratio per s.d.  
OR <- 1.50

# Simulate PRS values for cases (subjects 1: 50000) and controls (50001:100000)
g                <- rep(-9,100000)                                        # PRS w/o F5
g[1:50000]       <- log(OR) + rnorm(50000)
g[50001:100000]  <- rnorm(50000)
g1               <- rep(-9,100000)                                        # PRS w/ F5
g1[1:50000]      <- log(OR) + rnorm(50000)
g1[1:50000]      <- g1[1:50000] + rbinom(50000,2,2.4*0.03/(2.4*.03+.97)) 
g1[50001:100000] <- rnorm(50000)
g1[50001:100000] <- g1[50001:100000] + rbinom(50000,2,.03)
d                <- c(rep(1,50000),rep(0,50000))

# Create genotype categories as in Sara's paper and run logistic regression
gcat  <- cut(g ,quantile(g ,prob=c(0,.05,.15,.25,.75,.85,.95,1)))
g1cat <- cut(g1,quantile(g1,prob=c(0,.05,.15,.25,.75,.85,.95,1)))
out   <- summary(glm(d~factor(gcat) ,family=binomial()))$coefficients
out1  <- summary(glm(d~factor(g1cat),family=binomial()))$coefficients

# Odds ratio comparing top 5% to 20-75th %iles
exp(out[7,1]  - out[4,1])
exp(out1[7,1] - out1[4,1])

# PRS odds ratio per s.d.  
OR <- 1.50

# Simulate PRS values for cases (subjects 1: 50000) and controls (50001:100000)
g <- rep(-9,100000)
g[1:50000] <- log(OR) + rnorm(50000)
g[50001:100000] <- rnorm(50000)
d <- c(rep(1,50000),rep(0,50000))

# Create genotype categories as in Sara's paper and run logistic regression
gcat <- cut(g,quantile(g,prob=c(0,.05,.15,.25,.75,.85,.95,1)))
out  <- summary(glm(d~factor(gcat),family=binomial()))$coefficients

# Odds ratio comparing top 5% to 20-75th %iles
exp(out[7,1] - out[4,1])
