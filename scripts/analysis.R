library(ggplot2)
library(gridExtra)
library(reshape2)

no_cpaths = read.csv("/home/kdasilva/master_project/sevens/results/new_paths_mapping_mix_o104_iai39_bl21.csv")

# 1. Tableau des gènes

data = read.csv("/home/kdasilva/master_project/sevens/results/table_genes_mapping_mix_o104_iai39_bl21.csv",sep = ",")
rownames(data) = data[,1]
data = data[,-1]
data[is.na(data)] = 0
data = data[,c(1,8:10,12:14,2:5,7,6,11)]
cnames = c("K12","IAI39","O104","Sakai","SE15","Santai","RM8426")
names(cnames) = c("NC_000913.3","NC_011750.1","NC_018658.1","NC_002695.2","NC_013654.1","NZ_CP007592.1","NZ_CP028116.1")
colnames(data)[1:7] = cnames[colnames(data)[1:7]]

## histogramme abondance des gènes
ggplot(data[-which(data$mean_abund==0),], aes(x=mean_abund)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white",binwidth=1)+
  geom_density(alpha=.2, fill="#FF6666") 

## sans gènes dupliqués
to_del = c()
for (i in 1:nrow(data)) {
  if (sum(data[i,1:7]>1)>0) {
    to_del = c(to_del,i)
  }
}
data2 = data[-to_del,]
ggplot(data2[-which(data2$mean_abund==0),], aes(x=mean_abund)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white",binwidth=1)+
  geom_density(alpha=.2, fill="#FF6666") 

# 1.1. tous les gènes considérés
#data

# 1.2. sans les gènes duppliqués
data = data2

data_strains = data.frame(matrix(0,7,length(seq(0,1,0.05))))
colnames(data_strains) = seq(0,1,0.05)
data_strains$strain = colnames(data)[1:7]
data_strains = melt(data_strains)
colnames(data_strains)[2] = "thr_cov"
data_strains$thr_cov = as.numeric(as.character(data_strains$thr_cov))
data_strains = data_strains[,-3]

## proportion de gènes détectés

data_strains$hl_uniq_genes_pc = 0
data_strains$hl_total_genes_pc = 0
for (i in 1:nrow(data_strains)) {
  s = data_strains$strain[i]
  thr = data_strains$thr_cov[i]
  data_strains$hl_uniq_genes_pc[i] = sum(data[which(data$coverage >= thr & data$mean_abund > 0),s] >0)/sum(data[,s] >0)
  data_strains$hl_total_genes_pc[i] = sum(data[which(data$coverage >= thr & data$mean_abund > 0),s])/sum(data[,s])
}

ggplot(data=data_strains, aes(x=thr_cov, y=hl_total_genes_pc, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Proportion de gènes détectés")

## abondance moyenne des gènes

### tous les gènes considérés même ceux à zéro
data_strains$mean_abund = 0
for (i in 1:nrow(data_strains)) {
  s = data_strains$strain[i]
  thr = data_strains$thr_cov[i]
  data_strains$mean_abund[i] = sum( data$mean_abund[which(data$coverage >= thr & data[,s] > 0)]/data[which(data$coverage >= thr & data[,s] > 0),s] )/sum(data[,s])
}

ggplot(data=data_strains, aes(x=thr_cov, y=mean_abund, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Abondance moyenne des gènes")

for (i in levels(factor(data_strains$thr_cov))) {
  print(data_strains$mean_abund[which(data_strains$strain=="O104" & data_strains$thr_cov==i)]/data_strains$mean_abund[which(data_strains$strain=="IAI39" & data_strains$thr_cov==i)])
}

### tous les gènes considérés sans ceux à zéro
data_strains$mean_abund = 0
for (i in 1:nrow(data_strains)) {
  s = data_strains$strain[i]
  thr = data_strains$thr_cov[i]
  data_strains$mean_abund[i] = mean( data$mean_abund[which(data$coverage >= thr & data[,s] > 0)]/data[which(data$coverage >= thr & data[,s] > 0),s] )
}

ggplot(data=data_strains, aes(x=thr_cov, y=mean_abund, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Abondance moyenne des gènes")

### moyenne uniquement sur les gènes détectés
data_strains$mean_abund = 0
for (i in 1:nrow(data_strains)) {
  s = data_strains$strain[i]
  thr = data_strains$thr_cov[i]
  data_strains$mean_abund[i] = mean( data$mean_abund[which(data$coverage >= thr & data[,s] > 0)]/data[which(data$coverage >= thr & data[,s] > 0),s] )
}

ggplot(data=data_strains, aes(x=thr_cov, y=mean_abund, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Abondance moyenne des gènes")

for (i in levels(factor(data_strains$thr_cov))) {
  print(data_strains$mean_abund[which(data_strains$strain=="O104" & data_strains$thr_cov==i)]/( data_strains$mean_abund[which(data_strains$strain=="IAI39" & data_strains$thr_cov==i)] + data_strains$mean_abund[which(data_strains$strain=="O104" & data_strains$thr_cov==i)] + data_strains$mean_abund[which(data_strains$strain=="K12" & data_strains$thr_cov==i)]) ) 
}
for (i in levels(factor(data_strains$thr_cov))) {
  print(data_strains$mean_abund[which(data_strains$strain=="IAI39" & data_strains$thr_cov==i)]/( data_strains$mean_abund[which(data_strains$strain=="IAI39" & data_strains$thr_cov==i)] + data_strains$mean_abund[which(data_strains$strain=="O104" & data_strains$thr_cov==i)] + data_strains$mean_abund[which(data_strains$strain=="K12" & data_strains$thr_cov==i)]) ) 
}
for (i in levels(factor(data_strains$thr_cov))) {
  print(data_strains$mean_abund[which(data_strains$strain=="K12" & data_strains$thr_cov==i)]/( data_strains$mean_abund[which(data_strains$strain=="IAI39" & data_strains$thr_cov==i)] + data_strains$mean_abund[which(data_strains$strain=="O104" & data_strains$thr_cov==i)] + data_strains$mean_abund[which(data_strains$strain=="K12" & data_strains$thr_cov==i)]) ) 
}

## modèle
library(glmnet)
mod <- cv.glmnet(as.matrix(data[,1:7]),as.numeric(data$uniq_count>0),intercept=F,alpha=1)
coef(mod,s="lambda.min") 

lm = lm(I(uniq_count + 0)~ 0 + O104 + IAI39 + K12,data)
coef(lm)
coef(lm)[1]/sum(coef(lm))
coef(lm)[2]/sum(coef(lm))
coef(lm)[3]/sum(coef(lm))

ggplot(data=data, aes(x=thr_cov, y=p_subst, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Proportion de substitutions (comptage unique)")









data = read.csv("/home/kdasilva/master_project/sevens/results/res_strainslevel_mapping_mix_o104_iai39_bl21.csv",sep = ",")
data = data[,-1]
cnames = c("K12","IAI39","O104","Sakai","SE15","Santai","RM8426")
names(cnames) = c("NC_000913.3","NC_011750.1","NC_018658.1","NC_002695.2","NC_013654.1","NZ_CP007592.1","NZ_CP028116.1")
data$strain = cnames[as.character(data$strain)]

nb_uniq_genes = c(5446, 4666, 5227, 4825, 4307, 4286, 4961)
names(nb_uniq_genes) = c('NZ_CP028116.1','NZ_CP007592.1','NC_002695.2','NC_011750.1','NC_000913.3','NC_013654.1','NC_018658.1')
names(nb_uniq_genes) = cnames[as.character(names(nb_uniq_genes))]

nb_total_genes = c(5546, 4693, 5299, 4863, 4306, 4985,4319)
names(nb_total_genes) = c('NZ_CP028116.1', 'NZ_CP007592.1', 'NC_002695.2', 'NC_011750.1', 'NC_013654.1', 'NC_018658.1', 'NC_000913.3')
names(nb_total_genes) = cnames[as.character(names(nb_total_genes))]


data$X0_norm = 0
for (i in levels(factor(data$thr_cov))) {
  data$X0_norm[which(data$thr_cov==i)] = data$X0[which(data$thr_cov==i)]/sum(data$X0[which(data$thr_cov==i)])*100
}

data$X1_norm = 0
for (i in levels(factor(data$thr_cov))) {
  data$X1_norm[which(data$thr_cov==i)] = data$X1[which(data$thr_cov==i)]/sum(data$X1[which(data$thr_cov==i)])*100
}

# uniq_genes %
data$hl_genes_uniq_pc = data$hl_genes_uniq
for (i in 1:nrow(data)) {
  data$hl_genes_uniq_pc[i] = data$hl_genes_uniq_pc[i]/nb_uniq_genes[data$strain[i]]
}

# total_genes %
data$hl_genes_total_pc = data$hl_genes_total
for (i in 1:nrow(data)) {
  data$hl_genes_total_pc[i] = data$hl_genes_total_pc[i]/nb_total_genes[data$strain[i]]
}

pdf(file = "/home/kdasilva/master_project/sevens/results/res_mapping_mix_o104_iai39_bl21.pdf", width = 9, height = 9, pointsize = 10)

ggplot(data=data, aes(x=thr_cov, y=hl_genes_uniq, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Nombre de gènes (uniques) détectés")

ggplot(data=data, aes(x=thr_cov, y=hl_genes_uniq_pc, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Proportion de gènes (uniques) détectés")

ggplot(data=data, aes(x=thr_cov, y=hl_genes_total, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Nombre de gènes (totaux) détectés")

ggplot(data=data, aes(x=thr_cov, y=hl_genes_total_pc, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Proportion de gènes (totaux) détectés")

ggplot(data=data, aes(x=thr_cov, y=mean_abund, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Abondance moyenne des gènes")

ggplot(data=data, aes(x=thr_cov, y=p_subst, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Proportion de substitutions (comptage unique)")

ggplot(data=data, aes(x=thr_cov, y=X0, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Nombre de reads sans mutations (comptage unique)")

ggplot(data=data, aes(x=thr_cov, y=X0_norm, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Nombre de reads sans mutations (comptage unique) en %")

ggplot(data=data, aes(x=thr_cov, y=X0+X1, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Nombre de reads présentant jusqu'à 1 substitution (comptage unique)")

ggplot(data=data, aes(x=thr_cov, y=X0_norm+X1_norm, color=strain)) +
  geom_point() + geom_line() + xlab("Seuil de couverture des gènes") + ylab("Nombre de reads présentant jusqu'à 1 substitution (comptage unique) en %")

dev.off()
