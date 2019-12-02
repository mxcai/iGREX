library(reshape2)
library(ggplot2)
library(gridExtra)
setwd("./results")

############################################### Fig 1e: 012 mat: auto-correlation ###############################################
out1 <- read.table("geno_rho0.1_5prop23_SNRy0.3_n1000_4000.txt",header = T)
out2 <- read.table("geno_rho0.3_5prop23_SNRy0.3_n1000_4000.txt",header = T)
out3 <- read.table("geno_rho0.5_5prop23_SNRy0.3_n1000_4000.txt",header = T)
out4 <- read.table("geno_rho0.8_5prop23_SNRy0.3_n1000_4000.txt",header = T)

# rho <- c(replicate(3,c(0.1,0.3,0.5,0.8)))
rho <- c(0.1,0.3,0.5,0.8)
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:4){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:10],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",REML0="REML0",MoM="MoM",MoM0="MoM0",ss="IGREX-s"),
               component=ifelse(substr(out[,1],4,4)=="g","GREX","Alternative"),rho=rho[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVEt","method","component","rho")
out_all$method <- factor(out_all$method,levels = c("REML0","MoM0","REML","MoM","IGREX-s"))


Pd  <- ggplot(data=out_all,aes(x=component,y=PVEt,fill=method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(.~rho, labeller = label_bquote(cols=rho*"="*.(rho))) +
  scale_fill_brewer(palette="RdBu",labels=c(REML="REML",REML0=bquote(REML[0]),Mom="MoM",MoM0=bquote(MoM[0]),SummaryStat="IGREX-s")) +
  coord_cartesian(ylim = c(0.1,0.45)) +
  labs(title="d",y=expression(PVE))+
  theme(text = element_text(size=20))
pdf("Fig1_geno_autocorr.pdf",width=26,height=10)
Pd
dev.off()


################################################ Supp Fig 1: SNRy ~ n1 ###############################################
out1 <- read.table("geno_5prop23_SNRy0.3_n800_4000.txt",header = T)
out2 <- read.table("geno_5prop23_SNRy0.3_n1000_4000.txt",header = T)
out3 <- read.table("geno_5prop23_SNRy0.3_n2000_4000.txt",header = T)
out4 <- read.table("geno_5prop23_SNRy0.2_n800_4000.txt",header = T)
out5 <- read.table("geno_5prop23_SNRy0.2_n1000_4000.txt",header = T)
out6 <- read.table("geno_5prop23_SNRy0.2_n2000_4000.txt",header = T)
out7 <- read.table("geno_5prop23_SNRy0.1_n800_4000.txt",header = T)
out8 <- read.table("geno_5prop23_SNRy0.1_n1000_4000.txt",header = T)
out9 <- read.table("geno_5prop23_SNRy0.1_n2000_4000.txt",header = T)

n1 <- c(replicate(3,c(800,1000,2000)))
SNRy <- c(sapply(c(0.3,0.2,0.1),rep,3))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:9){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:10],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",REML0="REML0",MoM="MoM",MoM0="MoM0",ss="IGREX-s"),
               component=ifelse(substr(out[,1],4,4)=="g","GREX","Alternative"),n1=n1[i],SNRy=SNRy[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVEt","method","component","n1","PVEy")
out_all$method <- factor(out_all$method,levels = c("REML0","MoM0","REML","MoM","IGREX-s"))


Pa <- ggplot(data=out_all,aes(x=component,y=PVEt,fill=method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(PVEy~n1, labeller = label_bquote(cols="n"[r]*"="*.(n1),rows="PVE"[y]*"="*.(PVEy))) +
  scale_fill_brewer(palette="RdBu",labels=c(REML="REML",REML0=bquote(REML[0]),Mom="MoM",MoM0=bquote(MoM[0]))) +
  labs(title="a",y=expression(PVE))+
  theme(text = element_text(size=20))
pdf("Fig1_geno_PVEy_n1.pdf",width=13,height=10)
Pa
dev.off()


################################################ Supp Fig 2: PVEg ~ n1 ###############################################
out1 <- read.table("5prop14_SNRy0.3_n800_4000.txt",header = T)
out2 <- read.table("5prop14_SNRy0.3_n1000_4000.txt",header = T)
out3 <- read.table("5prop14_SNRy0.3_n2000_4000.txt",header = T)
out4 <- read.table("5prop23_SNRy0.3_n800_4000.txt",header = T)
out5 <- read.table("5prop23_SNRy0.3_n1000_4000.txt",header = T)
out6 <- read.table("5prop23_SNRy0.3_n2000_4000.txt",header = T)
out7 <- read.table("5prop32_SNRy0.3_n800_4000.txt",header = T)
out8 <- read.table("5prop32_SNRy0.3_n1000_4000.txt",header = T)
out9 <- read.table("5prop32_SNRy0.3_n2000_4000.txt",header = T)
out10 <- read.table("5prop41_SNRy0.3_n800_4000.txt",header = T)
out11 <- read.table("5prop41_SNRy0.3_n1000_4000.txt",header = T)
out12 <- read.table("5prop41_SNRy0.3_n2000_4000.txt",header = T)

n1 <- c(replicate(4,c(800,1000,2000)))
PVEg <- c(sapply(c(0.1,0.2,0.3,0.4),rep,3))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:12){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:10],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",REML0="REML0",MoM="MoM",MoM0="MoM0",ss="IGREX-s"),
               component=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),n1=n1[i],PVEg=PVEg[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","component","n1","PVEg")
out_all$method <- factor(out_all$method,levels = c("REML0","MoM0","REML","MoM","IGREX-s"))
levels(out_all$component) <- c("GREX","Alternative")

pdf("Supp_PVEg_n1.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=component,y=PVE,fill=method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = PVEg),linetype="dashed",color="blue") +
  geom_hline(aes(yintercept = 0.5-PVEg),linetype="dashed",color="red") +
  facet_grid(PVEg~n1, labeller = label_bquote(cols="n"[r]*"="*.(n1),rows="PVE"[GREX]*"="*.(PVEg))) +
  scale_fill_brewer(palette="RdBu",labels=c(REML="REML",REML0=bquote(REML[0]),Mom="MoM",MoM0=bquote(MoM[0]))) +
  labs(y=expression(PVE))+
  theme(text = element_text(size=20),legend.position = "bottom")
P
dev.off()


############################################### Supp Fig 3: piG ~ n1 ###############################################
out1 <- read.table("sparsityB0.2_sparsityG0.2_SNRy0.3_n800_4000.txt",header = T)
out2 <- read.table("sparsityB0.2_sparsityG0.5_SNRy0.3_n800_4000.txt",header = T)
out3 <- read.table("sparsityB0.2_sparsityG0.8_SNRy0.3_n800_4000.txt",header = T)
out4 <- read.table("sparsityB0.2_sparsityG0.2_SNRy0.3_n1000_4000.txt",header = T)
out5 <- read.table("sparsityB0.2_sparsityG0.5_SNRy0.3_n1000_4000.txt",header = T)
out6 <- read.table("sparsityB0.2_sparsityG0.8_SNRy0.3_n1000_4000.txt",header = T)
out7 <- read.table("sparsityB0.2_sparsityG0.2_SNRy0.3_n2000_4000.txt",header = T)
out8 <- read.table("sparsityB0.2_sparsityG0.5_SNRy0.3_n2000_4000.txt",header = T)
out9 <- read.table("sparsityB0.2_sparsityG0.8_SNRy0.3_n2000_4000.txt",header = T)

sparsityG <- c(replicate(4,c(0.2,0.5,0.8)))
n1 <- c(sapply(c(800,1000,2000),rep,3))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:9){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:6],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",MoM="MoM",ss="IGREX-s"),
               component=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),sparsityG=sparsityG[i],n1=n1[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","component","sparsity","n1")
out_all$method <- factor(out_all$method,levels = c("REML","MoM","IGREX-s"))
levels(out_all$component) <- c("GREX","Alternative")
pdf("Supp_piG_n1.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=component,y=PVE,fill=method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(sparsity~n1, labeller = label_bquote(cols="n"[r]*"="*.(n1),rows=pi[alpha]*"="*.(sparsity))) +
  scale_fill_brewer(palette="RdBu",labels=c(REML="REML",REML0=bquote(REML[0]),Mom="MoM",MoM0=bquote(MoM[0]))) +
  labs(y=expression(PVE))+
  theme(text = element_text(size=20),legend.position = "bottom")
P
dev.off()



################################################ Supp Fig 4: piB ~ SNRy ###############################################
out1 <- read.table("sparsityB0.2_SNRy0.1_n800_4000.txt",header = T)
out2 <- read.table("sparsityB0.5_SNRy0.1_n800_4000.txt",header = T)
out3 <- read.table("sparsityB0.8_SNRy0.1_n800_4000.txt",header = T)
out4 <- read.table("sparsityB0.2_SNRy0.2_n800_4000.txt",header = T)
out5 <- read.table("sparsityB0.5_SNRy0.2_n800_4000.txt",header = T)
out6 <- read.table("sparsityB0.8_SNRy0.2_n800_4000.txt",header = T)
out7 <- read.table("sparsityB0.2_SNRy0.3_n800_4000.txt",header = T)
out8 <- read.table("sparsityB0.5_SNRy0.3_n800_4000.txt",header = T)
out9 <- read.table("sparsityB0.8_SNRy0.3_n800_4000.txt",header = T)
out10 <- read.table("sparsityB0.2_SNRy0.5_n800_4000.txt",header = T)
out11 <- read.table("sparsityB0.5_SNRy0.5_n800_4000.txt",header = T)
out12 <- read.table("sparsityB0.8_SNRy0.5_n800_4000.txt",header = T)

sparsity <- c(replicate(4,c(0.2,0.5,0.8)))
SNRy <- c(sapply(c(0.1,0.2,0.3,0.5),rep,3))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:12){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:6],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",MoM="MoM",ss="IGREX-s"),
               component=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),sparsity=sparsity[i],SNRy=SNRy[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","component","sparsity","SNRy")
out_all$method <- factor(out_all$method,levels = c("REML","MoM","IGREX-s"))
levels(out_all$component) <- c("GREX","Alternative")
pdf("Supp_piB_PVEy.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=component,y=PVE,fill=method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(sparsity~SNRy, labeller = label_bquote(cols="PVE"[y]*"="*.(SNRy),rows=pi[beta ]*"="*.(sparsity))) +
  scale_fill_brewer(palette="RdBu",labels=c(REML="REML",REML0=bquote(REML[0]),Mom="MoM",MoM0=bquote(MoM[0]))) +
  labs(y=expression(PVE))+
  theme(text = element_text(size=20),legend.position = "bottom")
P
dev.off()

################################################ Supp Fig 5: 012 mat: m ~ n1 ###############################################
out1 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample100_prop23_SNRy0.3_n800_4000.txt",header = T)
out2 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample300_prop23_SNRy0.3_n800_4000.txt",header = T)
out3 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample500_prop23_SNRy0.3_n800_4000.txt",header = T)
out4 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample1000_prop23_SNRy0.3_n800_4000.txt",header = T)
out5 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample100_prop23_SNRy0.3_n1000_4000.txt",header = T)
out6 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample300_prop23_SNRy0.3_n1000_4000.txt",header = T)
out7 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample500_prop23_SNRy0.3_n1000_4000.txt",header = T)
out8 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample1000_prop23_SNRy0.3_n1000_4000.txt",header = T)
out9 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample100_prop23_SNRy0.3_n2000_4000.txt",header = T)
out10 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample300_prop23_SNRy0.3_n2000_4000.txt",header = T)
out11 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample500_prop23_SNRy0.3_n2000_4000.txt",header = T)
out12 <- read.table("/Users/cmx/Desktop/Predixcan/mHeritability/simulation_results/MoM_subsample1000_prop23_SNRy0.3_n2000_4000.txt",header = T)

subsample <- c(replicate(3,c(100,300,500,1000)))
n1 <- c(sapply(c(800,1000,2000),rep,4))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:12){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,1:2],id.vars =NULL)
  out <- cbind(out,method=sapply(substring(out[,1],6),switch,REML="REML",REML0="REML0",MoM="MoM",MoM0="MoM0",ss="SummaryStat"),
               component=ifelse(substr(out[,1],4,4)=="g","gene","SNP"),subsample=subsample[i],n1=n1[i])
  out <- out[,-1]
  out_all <- rbind(out_all,out)
}

names(out_all) <- c("PVE","method","component","subsample","n1")
out_all$method <- factor(out_all$method,levels = c("REML0","MoM0","REML","MoM","SummaryStat"))
levels(out_all$component) <- c("GREX","Alternative")

pdf("Supp_geno_m_n1.pdf",width=13,height=10)
P  <- ggplot(data=out_all,aes(x=component,y=PVE)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.2,linetype="dashed",color="blue") +
  geom_hline(yintercept = 0.3,linetype="dashed",color="red") +
  facet_grid(n1~subsample, labeller = label_bquote(cols="m"*"="*.(subsample),rows="n"[r]*"="*.(n1))) +
  scale_fill_brewer(palette="RdBu") +
  labs(y=expression(PVE))+
  theme(text = element_text(size=20),legend.position = "bottom")
P
dev.off()

################################################ Supp Fig 6: RHOGE: n1 ###############################################
library(ggplot2)
library(reshape2)
out1 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.1_n800_4000.txt",header = T)
out2 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.3_n800_4000.txt",header = T)
out3 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.5_n800_4000.txt",header = T)
out4 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.7_n800_4000.txt",header = T)
out5 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.9_n800_4000.txt",header = T)
out6 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.1_n1000_4000.txt",header = T)
out7 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.3_n1000_4000.txt",header = T)
out8 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.5_n1000_4000.txt",header = T)
out9 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.7_n1000_4000.txt",header = T)
out10 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.9_n1000_4000.txt",header = T)
out11 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.1_n2000_4000.txt",header = T)
out12 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.3_n2000_4000.txt",header = T)
out13 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.5_n2000_4000.txt",header = T)
out14 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.7_n2000_4000.txt",header = T)
out15 <- read.table("iGREXvsRHOGE_SNRz0.2_SNRy0.9_n2000_4000.txt",header = T)

SNRy <- c(replicate(6,c(0.1,0.3,0.5,0.7,0.9)))
n1 <- c(replicate(2,c(sapply(c(800,1000,2000),rep,5))))
n_rep <- nrow(out1)

out_all <- data.frame()
for(i in 1:15){
  out <- get(paste("out",i,sep=""))
  out <- melt(out[,c(1,2,3,5)],id.vars =NULL)
  out <- cbind(out,SNRy=SNRy[i],n1=n1[i])
  out_all <- rbind(out_all,out)
}
names(out_all) <- c("method","PVEy","PVE","n1")
out_all$PVE <- as.factor(out_all$PVE)
levels(out_all$method) <- c("REML","MoM","IGREX-s","RhoGE")


Pe  <- ggplot(data=out_all,aes(x=PVE,y=PVEy,fill=method)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = 0.2),linetype="dashed",color="blue") +
  facet_grid(.~n1,scales="free", labeller = label_bquote(cols="n"[r]*"="*.(n1))) +
  scale_fill_brewer(palette="RdBu") +
  labs(title="e",x=expression(PVE[y]), y=expression(PVE[GREX]))+
  theme(text = element_text(size=20))
pdf("/Users/cmx/Desktop/Predixcan/iGREX_paper/Fig1_RHOGE_n1_PVEt.pdf",width=26,height=10)
Pe
dev.off()

