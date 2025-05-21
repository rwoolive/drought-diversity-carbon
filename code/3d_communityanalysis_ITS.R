

# Separate microbial composition analyses across timepoints prior to analysis.
# 1. Statistically estimate treatment effects on soil microbial communities using PERMANOVA. 
# 2. Use distance-based redundancy analyses to visualize microbial community compositions 
# across treatments and to quantify variation in microbial community composition explained 
# by soil chemical and functional properties. 
# 3. Linear discriminant analysis effect size (LEfSe) method was used to identify microbial 
# taxonomic groups whose abundances were associated with specific residue type and/or 
# drought treatments (Segata et al. 2011).
# 4. Use Pearson correlations to determine relationships between soil properties and soil 
# taxa identified by the LEfSe analysis. 


library(sjmisc)
library(vegan)
library(QsRutils)
library(tidyr)
library(dplyr)
library(ecodist)
library(ape)
library(lme4)
library(ggplot2)
library(emmeans)
library(nlme)
library(plyr)
library(multcomp)
library(ggpubr)
#install.packages("corrplot") 
library(corrplot)
#BiocManager::install("lefser")
library(lefser)




covercols <- c("chocolate", "goldenrod", "forestgreen", "darkcyan", "darkorchid")
mcols <- c("blue", "red")
hcols <- c("green", "hotpink")


#### NMDS & PERMANOVA: treatment effects on soil microbial communities


taxabund <- read.csv("raw-data/sequence/ITS/taxabund_rarefied.csv")

# asv abundances
abundances <- dplyr::select(taxabund, contains("Wooliver_ITS"))
#abundances <- readRDS(paste0("raw-data/sequence/ITS/seqtab/5seqtab.rds"))
abundances <- t(as.data.frame(abundances))
# taxonomic data
taxa <- (taxabund[,3:13])
#taxa <- read.csv(paste0("raw-data/sequence/ITS/taxonomy/5_taxonomy-assignment-with-ft.csv"))
taxa$asv_number <- as.numeric(taxa$X.1)

colnames(abundances) <- taxa$asv_number
divmat_pa <- proportions(as.matrix(abundances),1)
divmat_pa <- as.data.frame(divmat_pa)
rownames(abundances)[1:3]

# sequence metadata
mdatseq <- read.csv("raw-data/master-dataset-13c_ITS_rarefied.csv")
# set groupings
mdatseq$Residue.type <- as.factor(mdatseq$Residue.type)
mdatseq$Residue.type <- factor(mdatseq$Residue.type, levels(mdatseq$Residue.type)[rev(c(2,5,1,4,3))])
#levels(mdatseq$Residue.type) <- c("No residue", "Wheat", "Clover", "Wheat-clover mix", "Five-species mix")
mdatseq$Moisture.treatment <- factor(mdatseq$Moisture.treatment)
mdatseq$Harvest.day <- as.factor(mdatseq$Harvest.day)
mdatseq$Harvest.day2 <- mdatseq$Harvest.day
levels(mdatseq$Harvest.day2) <- c("24 hr", "6 mo")
mdatseq$Harvest.day2[which(is.na(mdatseq$Harvest.day2))] <- NA
mdatseq$file.names[1:3]
# duplicate 13C measurements from g-50 samples to paired g-00 samples
mdatseq <- mdatseq[order(mdatseq$Glucose.addition.treatment),]
mdatseq[which(mdatseq$Glucose.addition.treatment=="G-00"),c("X13soc", "X13poc", "X13maoc", "X13doc", "X13mbc", "X13co2", "CUE", "co2.native")] <- mdatseq[which(mdatseq$Glucose.addition.treatment=="G-50"),c("X13soc", "X13poc", "X13maoc", "X13doc", "X13mbc", "X13co2", "CUE", "co2.native")]
# order by sequence file names
mdatseq <- mdatseq[order(mdatseq$file.names),]
# remove sample that doesn't have diversity data
mdatseq <- mdatseq[-which(is.na(mdatseq$richness)==TRUE),] # sample number 93, ambient moisture, wheat-clover residue, glucose addition, 24 hr harvest
# mdatseq_treat <- mdatseq[,c("Harvest.day2", "Moisture.treatment", "Residue.type", "Lab.rep")]
# rownames(mdatseq_treat) <- mdatseq$file.names
# mdatseq_treat[1:3,]





### 1. permanova: subset to each season
mod.effects <- data.frame(effect=c("Moisture.treatment",
                                   "Residue.type",
                                   "Moisture.treatment:Residue.type",
                                   "Residual","Total")) 

# 24 hr
#see if glucose has an effect
divmat_pa_24hr <- divmat_pa[which(mdatseq$Harvest.day2=="24 hr"),]
mdatseq_24hr <- mdatseq[which(mdatseq$Harvest.day2=="24 hr"),]
permmod <- adonis2(divmat_pa_24hr ~ Moisture.treatment*Residue.type*Glucose.addition.treatment, 
                   data = mdatseq_24hr, 
                   method="bray", permutations = 199, 
                   strata = mdatseq_24hr$Lab.rep, 
                   by = "terms" )
permmod2 <- round(as.data.frame(permmod), 3)
write.csv(permmod2, paste0("Model-output/permanova/with glucose_ITS by timepoint_24 hr_rarefied.csv"))

# remove glucose from model
divmat_pa_24hr <- divmat_pa[which(mdatseq$Harvest.day2=="24 hr"),]
mdatseq_24hr <- mdatseq[which(mdatseq$Harvest.day2=="24 hr"),]
permmod <- adonis2(divmat_pa_24hr ~ Moisture.treatment*Residue.type, 
                   data = mdatseq_24hr, 
                   method="bray", permutations = 199, 
                   strata = mdatseq_24hr$Lab.rep, 
                   by = "terms" )
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", 
              permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'24hr' <- ssp



# 6 mo
#see if glucose has an effect
divmat_pa_6mo <- divmat_pa[which(mdatseq$Harvest.day2=="6 mo"),]
mdatseq_6mo <- mdatseq[which(mdatseq$Harvest.day2=="6 mo"),]
permmod <- adonis2(divmat_pa_6mo ~ Moisture.treatment*Residue.type*Glucose.addition.treatment, 
                   data = mdatseq_6mo, 
                   method="bray", permutations = 199, 
                   strata = mdatseq_6mo$Lab.rep, 
                   by = "terms" )
permmod2 <- round(as.data.frame(permmod), 3)
write.csv(permmod2, paste0("Model-output/permanova/with glucose_ITS by timepoint_6 mo_rarefied.csv"))

# remove glucose from model
divmat_pa_6mo <- divmat_pa[which(mdatseq$Harvest.day2=="6 mo"),]
mdatseq_6mo <- mdatseq[which(mdatseq$Harvest.day2=="6 mo"),]
permmod <- adonis2(divmat_pa_6mo ~ Moisture.treatment*Residue.type, 
                   data = mdatseq_6mo, 
                   method="bray", permutations = 199, 
                   strata = mdatseq_6mo$Lab.rep, 
                   by = "terms" )
permmod2 <- round(as.data.frame(permmod), 3)
# append to a mod.effects dataframe
ssp <- paste0("F=", permmod2$F, ", Df=", permmod2$Df, ", p=", 
              permmod2$`Pr(>F)`, ", R2=", permmod2$R2)
mod.effects$'6mo' <- ssp


write.csv(mod.effects, paste0("Model-output/permanova/ITS by timepoint_rarefied.csv"))





#### 2. dbRDA: microbial community composition and soil properties
mdatseq$X13residue <- mdatseq$X13soc - mdatseq$X13mbc
mdatseq$CSE <- (mdatseq$X13residue + mdatseq$X13mbc)/(mdatseq$X13residue + mdatseq$X13mbc + mdatseq$X13co2)
mdatseq$CSE[which(mdatseq$CSE<0)] <- NA
mdatseq$X13residue[which(mdatseq$X13residue<0)] <- NA


# 24hr rda

# # first determine which properties are co-linear so that we
# can exclude them from the rda if needed
divmat_pa_24hr <- divmat_pa[which(mdatseq$Harvest.day2=="24 hr"),]
mdatseq_24hr <- mdatseq[which(mdatseq$Harvest.day2=="24 hr"),]

# res <- cor(mdatseq_24hr[,c("CSE", "X13maoc", "X13poc", "X13co2", "X13mbc",
#                            "X13doc", "X13soc", "co2.native",
#                            "pom.c.g.kg",  "pom.cn", "pom.n.g.kg",
#                            "maom.c.g.kg", "maom.cn", "maom.n.g.kg",
#                            "bulk.soil.C.g.kg", "bulk.soil.cn",  "bulk.soil.N.g.kg",
#                            "doc", "co2", "mbc", "dna")], #
#            method="pearson", use="pairwise.complete.obs")
# res2 <- as.data.frame(round(res, 2))
# pdf("figures/4_correlations/soil-property-correlations-for-dbrda_24hr.pdf", width=12, height=12)
# corrplot(res, type = "upper", is.corr=T, method = "number", order = "hclust",
#          tl.col = "black", tl.srt = 45)
# dev.off()
# #inspect plot

myrda_24hr0 <- capscale(formula = divmat_pa_24hr ~CSE+ X13residue +X13maoc+ X13poc+ X13co2+ #X13mbc+ 
                          X13doc+ X13soc+ co2.native+
                          pom.c.g.kg+  pom.cn+ #pom.n.g.kg+
                          maom.c.g.kg+ maom.cn+ #maom.n.g.kg+ 
                          bulk.soil.C.g.kg+ #bulk.soil.cn+  #bulk.soil.N.g.kg+ 
                          doc+ co2, #+dna + mbc
                        distance = "bray", 
                        data = mdatseq_24hr, na.action = "na.omit")
#myrda_24hr <- ordistep(myrda_24hr0, perm.max = 50) 
myrda_24hr <- capscale(formula = divmat_pa_24hr ~ CSE + X13soc + pom.cn + maom.cn + bulk.soil.C.g.kg + co2, #X13co2 +  X13co2 + X13soc + co2 + maom.cn +  bulk.soil.C.g.kg, # +pom.cn+ CSE+
                       distance = "bray", 
                        data = mdatseq_24hr, na.action = "na.omit")
vif.cca(myrda_24hr); hist(vif.cca(myrda_24hr)) # Check that there are no redundant predictors in the model. A common rule is that values over 10 indicate redundant constraints.
mod <- round(as.data.frame(anova(myrda_24hr)), 3)
ssp <- paste0("F=", mod$F, ", Df=", mod$Df, ", p=", mod$`Pr(>F)`, ", R2=", round(mod$SumOfSqs/sum(mod$SumOfSqs),3))

##### append to a mod.effects dataframe
mod.effects <- data.frame(effect=c("model", "residual")) 

mod.effects$myrda_24hr <- ssp

## rda all
# factors associated with each axis
write.csv(as.data.frame(summary(myrda_24hr)$biplot), paste0("Model-output/db-rda/ITS_24hr_loadings_rarefied.csv"))
# site coordinates
write.csv(as.data.frame(summary(myrda_24hr)$site), paste0("Model-output/db-rda/ITS_24hr_axes_rarefied.csv"))
# variation explained by first two axes
write.csv(summary(myrda_24hr)$cont$importance, paste0("Model-output/db-rda/ITS_24hr_variance-explained_rarefied.csv"))
# variation explained by each soil property
mdatseq_24hr_nonas0 <- na.omit(mdatseq_24hr[,c("Residue.type", "Moisture.treatment", "Glucose.addition.treatment", "file.names",
                                               rownames(summary(myrda_24hr)$biplot))])
mdatseq_24hr_nonas <- na.omit(mdatseq_24hr[,c(rownames(summary(myrda_24hr)$biplot))]) 
ef <- envfit(myrda_24hr, mdatseq_24hr_nonas,
             choices = c(1,2), permutations = 0, na.rm=TRUE)
write.csv(round(cbind(ef$vectors$arrows, ef$vectors$r),2), paste0("Model-output/db-rda/ITS_24hr_variance-by-predictors_rarefied.csv"))



#### visualize: 
# with soil property eigenvectors
rdadat <- as.data.frame(summary(myrda_24hr)$sites) # sites
rdadat$Residue.type<- mdatseq_24hr_nonas0$Residue.type
rdadat$Moisture.treatment <- mdatseq_24hr_nonas0$Moisture.treatment
rdadat$Glucose.addition.treatment <- mdatseq_24hr_nonas0$Glucose.addition.treatment
rdadat$file.names <- mdatseq_24hr_nonas0$file.names
df2  <- data.frame(summary(myrda_24hr)$biplot) # loadings: only keep top rated for first two axes
df3 <- df2
labs2 <- c(expression(paste("CSE")),
              #expression(paste("MAO"^"13","C")),
              #expression(paste("PO"^"13","C")),
  #expression(paste("MB"^"13","C")),
  #expression(paste(" "^"13","CO"[2])),
  expression(paste("SO"^"13","C")),
  expression(paste("POM C:N")), 
  expression(paste("MAOM C:N")), 
  expression(paste("SOC")),  
  expression(paste("Cumulative CO"[2]))) 
  #expression(paste("DO"^"13","C")),
              #expression(paste("Priming")),
              #expression(paste("POM-C")), 
              #expression(paste("MAOM-C")), 
              #expression(paste("DOC")))
rdadat$Glucose.addition.treatment <- as.factor(rdadat$Glucose.addition.treatment)
levels(rdadat$Glucose.addition.treatment) <- c("No glucose", "Glucose")
summary(myrda_24hr)$cont$importance[2,1] # 0.4935619
summary(myrda_24hr)$cont$importance[2,2] # 0.06176933
db24hr <- ggplot(rdadat, aes(CAP1, CAP2)) +
  geom_point(aes(shape=Moisture.treatment, color=Residue.type, alpha=Glucose.addition.treatment)) +
  theme_bw() +
  scale_color_manual(values = covercols) +
  scale_alpha_manual(values = c(0.5,1)) +
  labs(x=("RDA 1 (49.4%)"), y="RDA 2 (6.1%)") + # amend according to variance explained
  lims(x=c(-1.2,1.5)) + 
  guides(shape=guide_legend(title="Moisture treatment"), color=guide_legend(title="Residue type"), alpha=guide_legend(title="Glucose addition")) +
  geom_segment(data=df3, aes(x=0, xend=CAP1, y=0, yend=CAP2), color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(inherit.aes=FALSE, data=df3, x=df3$CAP1,y=df3$CAP2,label=labs2, 
            hjust=0.5*(1-sign(df3$CAP1)),vjust=0.5*(1-sign(df3$CAP2)), 
            color="blue", size=2) +
  theme()
db24hr
pdf(paste0("figures/_communities/db-rda/ITS_24hr1_rarefied.pdf"), height=3.75, width=5.5)
db24hr
dev.off()






# 6mo rda


# # first determine which properties are co-linear so that we
# can exclude them from the rda if needed
divmat_pa_6mo <- divmat_pa[which(mdatseq$Harvest.day2=="6 mo"),]
mdatseq_6mo <- mdatseq[which(mdatseq$Harvest.day2=="6 mo"),]

# find NAs in data
drop <- mdatseq_6mo$file.names[which(is.na(mdatseq_6mo$CSE))]

# drop from both dfs
divmat_pa_6mo <- divmat_pa_6mo[-which(rownames(divmat_pa_6mo) %in% drop),]
mdatseq_6mo <- mdatseq_6mo[-which(mdatseq_6mo$file.names %in% drop),]

# res <- cor(mdatseq_6mo[,c("CSE", "X13residue", "X13maoc", "X13poc", "X13co2", "X13mbc",
#                           "X13doc", "X13soc", "co2.native",
#                           "pom.c.g.kg",  "pom.cn", "pom.n.g.kg",
#                           "maom.c.g.kg", "maom.cn", "maom.n.g.kg",
#                           "bulk.soil.C.g.kg", "bulk.soil.cn",  "bulk.soil.N.g.kg",
#                           "doc", "co2", "mbc", "dna")],
#            method="pearson", use="pairwise.complete.obs")
# res2 <- as.data.frame(round(res, 2))
# pdf("figures/4_correlations/soil-property-correlations-for-dbrda_6mo.pdf", width=12, height=12)
# corrplot(res, type = "upper", is.corr=T, method = "number", order = "hclust",
#          tl.col = "black", tl.srt = 45)
# dev.off()
# #inspect plot

myrda_6mo0 <- capscale(formula = divmat_pa_6mo ~CSE+ X13maoc+ X13poc+ X13co2+ #X13mbc+ X13residue +
                         X13doc+ co2.native+ #X13soc+ 
                         pom.c.g.kg+  pom.cn+ #pom.n.g.kg+
                         maom.c.g.kg+ maom.cn+ #maom.n.g.kg+ 
                         bulk.soil.C.g.kg+ #bulk.soil.cn+  #bulk.soil.N.g.kg+ 
                         doc+ co2 + mbc, #+dna+
                       distance = "bray", 
                       data = mdatseq_6mo, na.action = "na.omit")
#myrda_6mo <- ordistep(myrda_6mo0, perm.max = 50) 
myrda_6mo <- capscale(formula = divmat_pa_6mo ~ X13maoc + maom.cn + bulk.soil.C.g.kg + co2, #doc + X13co2 + 
                      distance = "bray", 
                       data = mdatseq_6mo, na.action = "na.omit")
vif.cca(myrda_6mo); hist(vif.cca(myrda_6mo)) # Check that there are no redundant predictors in the model. A common rule is that values over 10 indicate redundant constraints.
mod <- round(as.data.frame(anova(myrda_6mo)), 3)
ssp <- paste0("F=", mod$F, ", Df=", mod$Df, ", p=", mod$`Pr(>F)`, ", R2=", round(mod$SumOfSqs/sum(mod$SumOfSqs),3))

##### append to a mod.effects dataframe
#mod.effects <- data.frame(effect=c("model", "residual")) 

mod.effects$myrda_6mo <- ssp
# export mod effects
write.csv(t(mod.effects), paste0("Model-output/db-rda/ITS_dbrda_rarefied.csv"))

## rda all
# factors associated with each axis
write.csv(as.data.frame(summary(myrda_6mo)$biplot), paste0("Model-output/db-rda/ITS_6mo_loadings_rarefied.csv"))
# site coordinates
write.csv(as.data.frame(summary(myrda_6mo)$site), paste0("Model-output/db-rda/ITS_6mo_axes_rarefied.csv"))
# variation explained by first two axes
write.csv(summary(myrda_6mo)$cont$importance, paste0("Model-output/db-rda/ITS_6mo_variance-explained_rarefied.csv"))
# variation explained by each soil property
mdatseq_6mo_nonas0 <- na.omit(mdatseq_6mo[,c("Residue.type", "Moisture.treatment", "Glucose.addition.treatment", "X13maoc","doc", "X13co2", "maom.cn", "bulk.soil.C.g.kg", "co2", #"X13mbc", 
                                             rownames(summary(myrda_6mo)$biplot) )]) 
mdatseq_6mo_nonas <- na.omit(mdatseq_6mo[,c(rownames(summary(myrda_6mo)$biplot) )]) 
ef <- envfit(myrda_6mo, mdatseq_6mo_nonas,
             choices = c(1,2), permutations = 0, na.rm=TRUE)
write.csv(round(cbind(ef$vectors$arrows, ef$vectors$r),2), paste0("Model-output/db-rda/ITS_6mo_variance-by-predictors_rarefied.csv"))



#### visualize: 
# with soil property eigenvectors
rdadat <- as.data.frame(summary(myrda_6mo)$sites) # sites
rdadat$Residue.type<- mdatseq_6mo_nonas0$Residue.type
rdadat$Moisture.treatment <- mdatseq_6mo_nonas0$Moisture.treatment
rdadat$Glucose.addition.treatment <- mdatseq_6mo_nonas0$Glucose.addition.treatment
df2  <- data.frame(summary(myrda_6mo)$biplot) # loadings: only keep top rated for first two axes
df3 <- df2
labs2 <- c(#expression(paste("CSE")), 
  expression(paste("MAO"^"13","C")), 
  #expression(paste("DOC")),
  #expression(paste(" "^"13","CO"[2])), 
  expression(paste("MAOM C:N")),  
  expression(paste("SOC")),
  expression(paste("Cumulative CO"[2]))) 
           #expression(paste("PO"^"13","C")),
           #expression(paste("DO"^"13","C")), 
           #expression(paste("SO"^"13","C")), 
           #expression(paste("Priming")),
           #expression(paste("POM-C")), 
           #expression(paste("POM C:N")), 
           #expression(paste("MAOM-C")), 
           #expression(paste("DOC")))
rdadat$Glucose.addition.treatment <- as.factor(rdadat$Glucose.addition.treatment)
levels(rdadat$Glucose.addition.treatment) <- c("No glucose", "Glucose")
summary(myrda_6mo)$cont$importance[2,1] #0.5249469
summary(myrda_6mo)$cont$importance[2,2] #0.05351638
db6mo <- ggplot(rdadat, aes(CAP1, CAP2)) +
  geom_point(aes(shape=Moisture.treatment, color=Residue.type, alpha=Glucose.addition.treatment)) +
  theme_bw() +
  scale_color_manual(values = covercols) +
  scale_alpha_manual(values = c(0.5,1)) +
  labs(x=("RDA 1 (52.5%)"), y="RDA 2 (5.5%)") + # amend according to variance explained
  lims(x=c(-1.3,1.5)) + guides(shape=guide_legend(title="Moisture treatment"), color=guide_legend(title="Residue type"), alpha=guide_legend(title="Glucose addition")) +
  geom_segment(data=df3, aes(x=0, xend=CAP1, y=0, yend=CAP2), color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(inherit.aes=FALSE, data=df3, x=df3$CAP1,y=df3$CAP2,label=labs2, 
            hjust=0.5*(1-sign(df3$CAP1)),vjust=0.5*(1-sign(df3$CAP2)), 
            color="blue", size=2)

db6mo
pdf(paste0("figures/_communities/db-rda/ITS_6mo1_rarefied.pdf"), height=3.75, width=5.5)
db6mo
dev.off()


ann1 <- ggplot() + 
  geom_text(aes(x=0, y=0, label ="A) Fungal community at 24 hr"), size = 4, hjust = 0.5) +
  theme_void()

ann2 <- ggplot() + 
  geom_text(aes(x=0, y=0, label ="B) Fungal community at 6 mo"), size = 4, hjust = 0.54) +
  theme_void()


# Fig: Ordination of microbial communities across different drought and residue type treatments 
# based on redundancy analysis. Arrow length and direction represent the strength of correlation 
# between soil properties and microbial communities. Results are shown for (A-B) fungal and 
# (C-D) bacterial communities of soils destructively harvested (A,C) 24 hours or (B,D) six 
# months after addition of glucose. 
db24hr <- db24hr +
  guides(color = guide_legend(override.aes = list(size=3)),
         alpha = guide_legend(override.aes = list(size=3)),
         shape = guide_legend(override.aes = list(size=3))) +
  labs(color="Residue type", shape="Moisture treatment", alpha="Glucose treatment")
db6mo <- db6mo  +
  guides(color = guide_legend(override.aes = list(size=3)),
         alpha = guide_legend(override.aes = list(size=3)),
         shape = guide_legend(override.aes = list(size=3)))+
  labs(color="Residue type", shape="Moisture treatment", alpha="Glucose treatment")

figure <- ggpubr::ggarrange(ann1, ann2, db24hr, db6mo, common.legend = T,legend="right",
                            ncol = 2, nrow = 2, heights=c(0.1,0.9)) 

figure
png("figures/_communities/db-rda/ITS_rarefied.png", 
    width = 2500, height =1100, res=300, units="px")
figure
dev.off()









#### 3. LEfSe: identified taxonomic groups whose abundances were 
# associated with specific treatments ** 24 hr **

# data("zeller14")
# ex_dat <- assays(zeller14)$exprs
# ex_names <- rownames(ex_dat)
# # rownames of my count data should look like this, with taxon names to the highest available specificity


# make sure count data and metadata are organized in the same order
data.frame(rownames(abundances), mdatseq$file.names)
# samples to remove
removes <- c(c("baseline1","baseline2","baseline3","baseline4"), 
             mdatseq$Sample.number[which(mdatseq$Harvest.day2=="6 mo")])
mdatseq_l <- mdatseq[-which(mdatseq$Sample.number %in% removes),]
abundances_t <- as.data.frame(t(abundances[-which(mdatseq$Sample.number %in% removes),]))

#### list row names of count data as taxa, with each row as a unique taxon
taxa2 <- taxa
taxa2$Genus <- paste0("g__", taxa2$Genus)
taxa2$Genus <- gsub(taxa2$Genus, pattern="g__NA", replacement="NA")
abundances_t$name <- paste(taxa2$Kingdom, taxa2$Phylum, taxa2$Class, taxa2$Order, taxa2$Family, 
                           taxa2$Genus, taxa2$Species, sep = "/") # , rownames(abundances_t)
abundances_t$name <- paste0(abundances_t$name, "/t__", taxa2$Primary_lifestyle) # , rownames(abundances_t)
abundances_t$name <- gsub(abundances_t$name, pattern="/NA", replacement="")
abundances_t$name <- gsub(abundances_t$name, pattern="/t_NA", replacement="")
dataname <- as.data.frame(abundances_t$name[1:50]) # check



# ######################################
# # different lefse approach: 
# library(microbiomeMarker)
# names(covercols) <- c("No residue", "Wheat", "Clover", "Wheat-clover mix", "Five-species mix")
# 
# # make sure count data and metadata are organized in the same order
# data.frame(rownames(abundances), mdatseq$file.names)
# # samples to remove
# removes <- c(c("baseline1","baseline2","baseline3","baseline4"))
# mdatseq_l <- mdatseq[-which(mdatseq$Sample.number %in% removes),]
# abundances_t <- as.data.frame(t(abundances[-which(mdatseq$Sample.number %in% removes),]))
# 
# #### list row names of count data as taxa, with each row as a unique taxon
# taxa2 <- taxa
# taxa2$Kingdom <- gsub("k__", "", taxa2$Kingdom)
# taxa2$Phylum <- gsub("p__", "", taxa2$Phylum)
# taxa2$Class <- gsub("c__", "", taxa2$Class)
# taxa2$Order <- gsub("o__", "", taxa2$Order)
# taxa2$Family <- gsub("f__", "", taxa2$Family)
# 
# 
# ### make a phyloseq object with:
# # asv_table: samples as columns, asvs as rows
# asv_table <- apply(abundances_t, 2, as.numeric)
# row.names(asv_table) <- taxa2$asv_number
# asv_table <- phyloseq::otu_table(as.matrix(asv_table), taxa_are_rows = T)
# # tax: clades as columns, asvs as rows
# tax <- taxa2[, 3:8]
# row.names(tax) <- taxa2$asv_number
# tax <- phyloseq::tax_table(as.matrix(tax))
# # sam: metadata as columns, samples as rows
# sam <- na.omit(mdatseq[c("file.names", "Moisture.treatment", "Residue.type", "Harvest.day")])
# row.names(sam) <- sam$file.names
# sam <- phyloseq::sample_data(sam)
# # phy_tree, and refseq: NA
# divexp_phy <- phyloseq::phyloseq(asv_table, tax, sam)
# 
# 
# ndis <- n_distinct(tax[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")])
# reso <- 400
# w_cutoff <- 0.05/ndis
# kw_cutoff <- 0.05/ndis
# l_cutoff <- 4.5
# 
# 
# 
# ### LEfSe for 24 hr data
# 
# # subset to timepoint
# divexp_phy_small <- phyloseq::subset_samples(
#   divexp_phy,
#   Harvest.day %in% c("T1")
# )
# # lefse analysis
# mm_lefse_0 <- run_lefse(divexp_phy_small, group = "Residue.type", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)
# 
# # which treatments are represented?
# gps <- unique(mm_lefse_0@marker_table$enrich_group)
# 
# # plot with cladogram
# p_overall <- plot_cladogram(mm_lefse_0, color = covercols[gps], 
#                             alpha = 0.5, clade_label_level = 5) +
#   theme(plot.margin = margin(0, 0, 0, 0)) #+   guides(fill = "none")  
# 
# ggpubr::ggexport(p_overall + labs(tag = "A) Fungi at 24 hr"), width = 3500, height = 2500, res = reso, 
#                  filename = paste0("figures/_lefse/phy_fungi_res_24hr.png"))
# 
# ### LEfSe for 6 mo data
# 
# # subset to timepoint
# divexp_phy_small <- phyloseq::subset_samples(
#   divexp_phy,
#   Harvest.day %in% c("T90")
# )
# # lefse analysis
# mm_lefse_0 <- run_lefse(divexp_phy_small, group = "Residue.type", wilcoxon_cutoff = w_cutoff, kw_cutoff = kw_cutoff, lda_cutoff = l_cutoff)
# 
# # which treatments are represented?
# gps <- unique(mm_lefse_0@marker_table$enrich_group)
# 
# # plot with cladogram
# p_overall <- plot_cladogram(mm_lefse_0, color = covercols[gps], 
#                             alpha = 0.5, clade_label_level = 5) +
#   theme(plot.margin = margin(0, 0, 0, 0)) +   guides(fill = "none")  
# 
# ggpubr::ggexport(p_overall + labs(tag = "B) Fungi at 6 mo"), width = 3500, height = 2500, res = reso, 
#                  filename = paste0("figures/_lefse/phy_fungi_res_6mo.png"))
# 
# 
# 
# 
# 
# ######################################



# sum abundances by taxon
abundances_t2 <- as.data.frame(abundances_t)
abundances_t3 <- abundances_t2 %>% 
  dplyr::group_by(name) %>%
  dplyr::summarise(across(where(is.numeric), sum, na.rm = TRUE))
abundances_t3 <- as.data.frame(abundances_t3)
abundances_t3$name <- gsub(abundances_t3$name, pattern="/", replacement="|")
rownames(abundances_t3) <- abundances_t3$name
abundances_t4 <- abundances_t3[,2:(dim(abundances_t3)[2])]
abundances_t <- abundances_t4
# export abundances by taxon
write.csv(abundances_t, "raw-data/sequence/ITS/Abundances by taxon_24hr.csv")


# dummy code residue treatment
mdatseq_l$Residue.type.nores <- as.factor(ifelse(mdatseq_l$Residue.type == 'No residue', 1, 0))
mdatseq_l$Residue.type.wheat <-  as.factor(ifelse(mdatseq_l$Residue.type == 'Wheat', 1, 0))
mdatseq_l$Residue.type.clove <-  as.factor(ifelse(mdatseq_l$Residue.type == 'Clover', 1, 0))
mdatseq_l$Residue.type.whecl <-  as.factor(ifelse(mdatseq_l$Residue.type == 'Wheat-clover mix', 1, 0))
mdatseq_l$Residue.type.fsm <-  as.factor(ifelse(mdatseq_l$Residue.type == 'Five-species mix', 1, 0))


#### lefse models by residue type for 24 hr 
sub <- which(mdatseq_l$Harvest.day2=="24 hr")
lef_24hr <- SummarizedExperiment(list(counts=abundances_t[,sub]), colData = mdatseq_l[sub,])

lef_24hr_nores2 <- lefser(lef_24hr, groupCol = "Residue.type.nores")
lefserPlot(lef_24hr_nores2, trim.names = FALSE ,colors = c(covercols[1], covercols[1]))
lef_24hr_nores2$Residue.type <- rep("No residue", dim(lef_24hr_nores2)[1])

lef_24hr_wheat2 <- lefser(lef_24hr, groupCol = "Residue.type.wheat")
lefserPlot(lef_24hr_wheat2, trim.names = FALSE ,colors = c(covercols[2], covercols[2]))
lef_24hr_wheat2$Residue.type <- rep("Wheat", dim(lef_24hr_wheat2)[1])

lef_24hr_clover2 <- lefser(lef_24hr, groupCol = "Residue.type.clove", )
lefserPlot(lef_24hr_clover2, trim.names = FALSE ,colors = c(covercols[3], covercols[3]))
lef_24hr_clover2$Residue.type <- rep("Clover", dim(lef_24hr_clover2)[1])

lef_24hr_whecl2 <- lefser(lef_24hr, groupCol = "Residue.type.whecl")
lefserPlot(lef_24hr_whecl2, trim.names = FALSE ,colors = c(covercols[4], covercols[4]))
lef_24hr_whecl2$Residue.type <- rep("Wheat-clover mix", dim(lef_24hr_whecl2)[1])

lef_24hr_fsm2 <- lefser(lef_24hr, groupCol = "Residue.type.fsm")
lefserPlot(lef_24hr_fsm2, trim.names = FALSE ,colors = c(covercols[5], covercols[5]))
lef_24hr_fsm2$Residue.type <- rep("Five-speices mix", dim(lef_24hr_fsm2)[1])

# merge all lefse results together
aa <- rbind(lef_24hr_nores2,lef_24hr_wheat2)
bb <- rbind(lef_24hr_clover2,lef_24hr_whecl2)
cc <- rbind(aa,bb)
dd <- rbind(cc,lef_24hr_fsm2)

# create columns for taxonomic group
dd$kingdom <- rep(NA, dim(dd)[1])
dd$phylum <- rep(NA, dim(dd)[1])
dd$class <- rep(NA, dim(dd)[1])
dd$order <- rep(NA, dim(dd)[1])
dd$family <- rep(NA, dim(dd)[1])
dd$genus <- rep(NA, dim(dd)[1])
dd$species <- rep(NA, dim(dd)[1])
dd$primarylifestyle <- rep(NA, dim(dd)[1])

for(i in 1:dim(dd)[1]){
  split_name <- strsplit(dd$features[i], "[|]")
  dd$kingdom[i] <- split_name[[1]][which(startsWith(split_name[[1]], "k__"))]
  dd$phylum[i] <- ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "p__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "p__"))])
  dd$class[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "c__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "c__"))])
  dd$order[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "o__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "o__"))])
  dd$family[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "f__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "f__"))])
  dd$genus[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "g__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "g__"))])
  dd$species[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "s__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "s__"))])
  dd$primarylifestyle[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "t_"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "t_"))])
}

dd$genus <- gsub(dd$genus, pattern = "g__", replacement="")
dd$species <- gsub(dd$species, pattern = "s__", replacement="")
dd$taxon <- paste(dd$genus, dd$species, sep=" ")
dd$taxon <- gsub(dd$taxon, pattern = " NA", replacement="")
dd$taxon[which(dd$taxon %in% c(" ", "NA"))] <- dd$family[which(dd$taxon %in% c(" ", "NA"))]
dd$taxon[which(dd$taxon %in% c(" ", NA))] <- dd$order[which(dd$taxon %in% c(" ", NA))]
dd$taxon[which(dd$taxon %in% c(" ", NA))] <- dd$class[which(dd$taxon %in% c(" ", NA))]
dd$taxon[which(dd$taxon %in% c(" ", NA))] <- dd$phylum[which(dd$taxon %in% c(" ", NA))]
dd$taxon[which(dd$taxon %in% c(" ", NA))] <- dd$kingdom[which(dd$taxon %in% c(" ", NA))]
dd$primarylifestyle <- gsub(dd$primarylifestyle, pattern="t__", replacement="")

dd$Residue.type <- as.factor(dd$Residue.type)
dd$Residue.type <- factor(dd$Residue.type, levels(dd$Residue.type)[c(3,4,1,5,2)])
pl <- unique(dd$primarylifestyle)

dd <- dd[order(dd$primarylifestyle, dd$taxon),]

# export lefse results 
write.csv(dd, "raw-data/sequence/ITS/Lefse_24hr_rarefied.csv")

# export  abundances for lefse taxa 
abundances_t_lefse <- abundances_t[which(rownames(abundances_t) %in% dd$features),]
write.csv(abundances_t_lefse, "raw-data/sequence/ITS/Abundances by taxon_24hr_Lefse_rarefied.csv")



# create a grid for heatmap with lefse taxa, effect sizes, and residue types
ddgrid <- expand.grid(Residue.type=levels(dd$Residue.type), taxon=dd$taxon)
ddgrid$score <- rep(NA, dim(ddgrid)[1])
ddgrid$primarylifestyle <- rep(NA, dim(ddgrid)[1])
for(i in 1:dim(ddgrid)[1]){
  score <-  dd$scores[which(dd$Residue.type==ddgrid$Residue.type[i] & dd$taxon==ddgrid$taxon[i])]
  ddgrid$score[i] <- ifelse(length(score)==0, NA, score)
  lifestyle <-  dd$primarylifestyle[which(dd$taxon==ddgrid$taxon[i])]
  ddgrid$primarylifestyle[i] <- ifelse(length(lifestyle)==0, NA, lifestyle)
}

ddgrid$primarylifestyle[which(ddgrid$primarylifestyle=="NA")] <- "unspecified"
ddgrid$primarylifestyle <- as.factor(ddgrid$primarylifestyle)
ddgrid$primarylifestyle <- factor(ddgrid$primarylifestyle, levels(ddgrid$primarylifestyle)[rev(c(2,4,6,1,3,5))])
ddgrid1 <- ddgrid[order(ddgrid$primarylifestyle, ddgrid$taxon),]
ddgrid1$taxon <- as.factor(ddgrid1$taxon)
ddgrid1$taxon <- factor(ddgrid1$taxon, unique(ddgrid1$taxon))

result <- ddgrid1 %>% 
  dplyr::group_by(primarylifestyle) %>% 
  dplyr::summarise(count = length(unique(taxon)), .groups = 'drop')
result$hal <- result$count/2
result <- result %>% mutate(cumulative = cumsum(count))
result$cumulativehal <- result$cumulative-result$hal

lab <- paste0(gsub(result$primarylifestyle, pattern="_", replacement=" "))
segs <- data.frame(x=5.5, xend=0.5,
                   y=result$cumulative+0.5, yend=result$cumulative+0.5,
                   label=lab, #gsub(lab, pattern=" ", replacement="\n"),
                   ylab=result$cumulativehal)
segs$ybar1 <- c(0.2, segs$yend[c(1:dim(segs)[1]-1)]+0.2)
segs$ybar2 <- c(segs$y-0.2)

ann1 <- ggplot() + 
  geom_text(aes(x=0.1, y=0, label ="A) Fungal taxa at 24 hr"), size = 4, hjust = 0.5) +
  theme_void()

ddgrid1$font <- rep(3,dim(ddgrid1)[1])
for(i in 1:dim(ddgrid1)[1]){
  if(str_contains(ddgrid1$taxon[i], "__")==TRUE) {ddgrid1$font[i] <- 1}
}


lefsep <-ggplot(ddgrid1, aes(fill= score)) + 
  geom_tile(color="black", aes(y=reorder(ddgrid1$taxon, -ddgrid1$primarylifestyle), x=Residue.type)) + 
  labs(y="")+
  guides(fill= guide_colorbar(title.position="top", title.hjust = 0.5, title="LEfSE score")) +
  scale_fill_gradient2(low="orchid2", mid = "white", high="seagreen3", guide="colorbar",na.value="white") +
  #scale_fill_gradient2(guide="colorbar",na.value="white", limits = c(-1, 1), low="orchid2", high="seagreen3", mid="white") +
  theme(axis.text.x = element_text(angle = 270, vjust = 0, hjust=0), 
        axis.text.y = element_text(size=7, face="italic"), 
        plot.margin = unit(c(0.1, 0.1, 0.1, 2.05), "cm"),
        axis.title.x=element_blank(), legend.position="top") +
  geom_segment(inherit.aes=FALSE, size=1, data=segs, aes(x=x, y=y, yend=yend, xend=xend)) +
  coord_cartesian(xlim=c(1,5), clip="off")+
  geom_segment(inherit.aes=FALSE, size=1, data=segs, aes(x=-5.5, xend=-5.5, y=ybar1, yend=ybar2)) +
  annotate("text", x = -5.75, y = segs$ylab+0.5, label = segs$label, size=3, hjust=1, vjust=0.5) #+
  #theme_minimal()
  
  #geom_text(inherit.aes=FALSE, data=segs,size=3, fontface=3, aes(x=x, y=y-0.5, label=segs$label), hjust = 1) +
  
lefsep

  figure <- ggpubr::ggarrange(ann1, lefsep, legend="top",
                              ncol = 1, nrow = 2, heights=c(0.1,0.9)) 
  figure
  png("figures/_lefse/1_ITS_24 hr_rarefied.png", 
      width = 1000, height =1500, res=300, units="px")
  figure
  dev.off()
  
  

#### 4. Pearson correlations: relationships between soil properties 
  #### from dbRDA and soil taxa from LEfSe **24 hr**
  # Table 3: Relationships between abundances of treatment-responsive microbial taxa identified 
  # with LEfSe algorithm and soil C-cycling properties identified with redundancy analysis. Results 
  # are shown for two timepoints: 24 hours and six months after glucose addition. Values represent 
  # Pearson correlation coefficients, and those with p-values < 0.05 are colored blue (positive) or 
  # red (negative). 
  
  # import  abundances for lefse taxa 
  abundances_t_lefse <- read.csv("raw-data/sequence/ITS/Abundances by taxon_24hr_Lefse_rarefied.csv")
  # order by primary lifestyle
  # create columns for taxonomic group
  abundances_t_lefse$kingdom <- rep(NA, dim(abundances_t_lefse)[1])
  abundances_t_lefse$phylum <- rep(NA, dim(abundances_t_lefse)[1])
  abundances_t_lefse$class <- rep(NA, dim(abundances_t_lefse)[1])
  abundances_t_lefse$order <- rep(NA, dim(abundances_t_lefse)[1])
  abundances_t_lefse$family <- rep(NA, dim(abundances_t_lefse)[1])
  abundances_t_lefse$genus <- rep(NA, dim(abundances_t_lefse)[1])
  abundances_t_lefse$species <- rep(NA, dim(abundances_t_lefse)[1])
  abundances_t_lefse$primarylifestyle <- rep(NA, dim(abundances_t_lefse)[1])
  for(i in 1:dim(abundances_t_lefse)[1]){
    split_name <- strsplit(abundances_t_lefse$X[i], "[|]")
    abundances_t_lefse$kingdom[i] <- split_name[[1]][which(startsWith(split_name[[1]], "k__"))]
    abundances_t_lefse$phylum[i] <- ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "p__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "p__"))])
    abundances_t_lefse$class[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "c__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "c__"))])
    abundances_t_lefse$order[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "o__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "o__"))])
    abundances_t_lefse$family[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "f__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "f__"))])
    abundances_t_lefse$genus[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "g__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "g__"))])
    abundances_t_lefse$species[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "s__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "s__"))])
    abundances_t_lefse$primarylifestyle[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "t_"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "t_"))])
  }
  abundances_t_lefse$genus <- gsub(abundances_t_lefse$genus, pattern = "g__", replacement="")
  abundances_t_lefse$species <- gsub(abundances_t_lefse$species, pattern = "s__", replacement="")
  abundances_t_lefse$taxon <- paste(abundances_t_lefse$genus, abundances_t_lefse$species, sep=" ")
  abundances_t_lefse$taxon <- gsub(abundances_t_lefse$taxon, pattern = " NA", replacement="")
  abundances_t_lefse$taxon[which(abundances_t_lefse$taxon %in% c(" ", "NA"))] <- abundances_t_lefse$family[which(abundances_t_lefse$taxon %in% c(" ", "NA"))]
  abundances_t_lefse$taxon[which(abundances_t_lefse$taxon %in% c(" ", NA))] <- abundances_t_lefse$order[which(abundances_t_lefse$taxon %in% c(" ", NA))]
  abundances_t_lefse$taxon[which(abundances_t_lefse$taxon %in% c(" ", NA))] <- abundances_t_lefse$class[which(abundances_t_lefse$taxon %in% c(" ", NA))]
  abundances_t_lefse$taxon[which(abundances_t_lefse$taxon %in% c(" ", NA))] <- abundances_t_lefse$phylum[which(abundances_t_lefse$taxon %in% c(" ", NA))]
  abundances_t_lefse$taxon[which(abundances_t_lefse$taxon %in% c(" ", NA))] <- abundances_t_lefse$kingdom[which(abundances_t_lefse$taxon %in% c(" ", NA))]
  abundances_t_lefse$primarylifestyle <- gsub(abundances_t_lefse$primarylifestyle, pattern="t__", replacement="")
  pl <- unique(abundances_t_lefse$primarylifestyle)
  abundances_t_lefse <- abundances_t_lefse[order(abundances_t_lefse$primarylifestyle, abundances_t_lefse$taxon),]
  
  
  
  # import  soil property data
  mdatseq <- read.csv("raw-data/master-dataset-13c_ITS_rarefied.csv")
  mdatseq$Residue.type <- as.factor(mdatseq$Residue.type)
  mdatseq$Residue.type <- factor(mdatseq$Residue.type, levels(mdatseq$Residue.type)[rev(c(2,5,1,4,3))])
  mdatseq$Moisture.treatment <- factor(mdatseq$Moisture.treatment)
  mdatseq$Harvest.day <- as.factor(mdatseq$Harvest.day)
  mdatseq$Harvest.day2 <- mdatseq$Harvest.day
  levels(mdatseq$Harvest.day2) <- c("24 hr", "6 mo")
  mdatseq$Harvest.day2[which(is.na(mdatseq$Harvest.day2))] <- NA
  # duplicate 13C measurements from g-50 samples to paired g-00 samples
  mdatseq <- mdatseq[order(mdatseq$Glucose.addition.treatment),]
  mdatseq[which(mdatseq$Glucose.addition.treatment=="G-00"),c("X13soc", "X13poc", "X13maoc", "X13doc", "X13mbc", "X13co2", "CUE", "co2.native")] <- mdatseq[which(mdatseq$Glucose.addition.treatment=="G-50"),c("X13soc", "X13poc", "X13maoc", "X13doc", "X13mbc", "X13co2", "CUE", "co2.native")]
  # order by sequence file names
  mdatseq <- mdatseq[order(mdatseq$file.names),]
  mdatseq <- mdatseq[-which(is.na(mdatseq$richness)==TRUE),] # sample number 93, ambient moisture, wheat-clover residue, glucose addition, 24 hr harvest
  mdatseq$X13residue <- mdatseq$X13soc - mdatseq$X13mbc
  mdatseq$CSE <- (mdatseq$X13residue + mdatseq$X13mbc)/(mdatseq$X13residue + mdatseq$X13mbc + mdatseq$X13co2)
  mdatseq$CSE[which(mdatseq$CSE<0)] <- NA
  mdatseq$X13residue[which(mdatseq$X13residue<0)] <- NA
  # subset to timepoint
  mdatseq_24hr <- mdatseq[which(mdatseq$Harvest.day2=="24 hr"),]
  
  abundances_t_lefse_t <- as.data.frame(t(abundances_t_lefse))
  colnames(abundances_t_lefse_t) <- abundances_t_lefse_t[dim(abundances_t_lefse_t)[1],]
  abundances_t_lefse_t2 <- abundances_t_lefse_t[rownames(abundances_t_lefse_t) %in% mdatseq_24hr$file.names,]
  abundances_t_lefse_t2$file.names <- rownames(abundances_t_lefse_t2)
  unique(mdatseq_24hr$file.names==abundances_t_lefse_t2$file.names) #check
  
# combine abundances and soil property data
newdat <- merge(mdatseq_24hr, abundances_t_lefse_t2, by="file.names",
                all.x = TRUE, all.y = TRUE)

terms.soil <- c(rownames(summary(myrda_24hr)$biplot)) 
terms.taxa <- colnames(abundances_t_lefse_t2)
terms.taxa <- terms.taxa[-which(terms.taxa %in% c("file.names"))]  



#newdat <- newdat[order(newdat$primarylifestyle, newdat$taxon),]


# create a grid for heatmap with lefse taxa, effect sizes, and residue types
ddgrid2 <- expand.grid(soil=colnames(newdat)[colnames(newdat) %in% terms.soil], taxon=colnames(newdat)[colnames(newdat) %in% terms.taxa])
ddgrid2$cor <- rep(NA, dim(ddgrid2)[1])
ddgrid2$pval <- rep(NA, dim(ddgrid2)[1])
ddgrid2$pvalindic <- rep(NA, dim(ddgrid2)[1])
ddgrid2$primarylifestyle <- rep(NA, dim(ddgrid2)[1])
for(i in 1:dim(ddgrid2)[1]){
  mod <- cor.test(newdat[,which(colnames(newdat)==ddgrid2$soil[i])], as.numeric(newdat[,which(colnames(newdat)==ddgrid2$taxon[i])]))
  ddgrid2$cor[i] <-  mod$estimate
  ddgrid2$pval[i] <- mod$p.value
  if(mod$p.value>0.05/length(ddgrid2$cor)){ddgrid2$cor[i] <- NA}
  if(mod$p.value<0.05/length(ddgrid2$cor)){ddgrid2$pvalindic[i] <- "*"}
  ddgrid2$significance.level[i] <- 0.05/length(ddgrid2$cor)  # this is the bon ferroni adjusted pvalue cutoff
  lifestyle <-  dd$primarylifestyle[which(dd$taxon==ddgrid2$taxon[i])]
  ddgrid2$primarylifestyle[i] <- ifelse(length(lifestyle)==0, NA, lifestyle)
}


ddgrid2$primarylifestyle[which(ddgrid2$primarylifestyle=="NA")] <- "unspecified"
ddgrid2$primarylifestyle <- as.factor(ddgrid2$primarylifestyle)
ddgrid2$primarylifestyle <- factor(ddgrid2$primarylifestyle, levels(ddgrid2$primarylifestyle)[rev(c(2,4,6,1,3,5))])#[rev(c(2,4,6,8,9,1,3,5,7))])
ddgrid3 <- ddgrid2[order(ddgrid2$primarylifestyle, ddgrid2$taxon),]
ddgrid3$taxon <- as.factor(ddgrid3$taxon)
ddgrid3$taxon <- factor(ddgrid3$taxon, unique(ddgrid3$taxon))

result <- ddgrid3 %>% 
  dplyr::group_by(primarylifestyle) %>% 
  dplyr::summarise(count = length(unique(taxon)), .groups = 'drop')
result <- result %>% mutate(cumulative = cumsum(count))

segs2 <- data.frame(x=length(unique(ddgrid3$soil))+0.5, xend=0,
                   y=result$cumulative+0.5, yend=result$cumulative+0.5,
                   label=paste0(gsub(result$primarylifestyle, pattern="_", replacement=" ")))

ann2 <- ggplot() + 
  geom_text(aes(x=0.1, y=0, label =""), size = 4, hjust = 0.5) +
  theme_void()

labs2 <- c(#expression(paste("CSE")),
  #expression(paste("MAO"^"13","C")), 
  #expression(paste("PO"^"13","C")),
  expression(paste("POM C:N")), 
  expression(paste("MAOM C:N")), 
  expression(paste("SOC")),  
  expression(paste("SO"^"13","C")), 
  #expression(paste(" "^"13","CO"[2])), 
  expression(paste("Cumulative CO"[2])),
  expression(paste("CSE")))
  #expression(paste("POM C:N")), 
  #expression(paste("DO"^"13","C")), 
  #expression(paste("Priming")),
  #expression(paste("POM-C")), 
  #expression(paste("MAOM-C")), 
  #expression(paste("DOC")))

p<-ggplot(ddgrid3, aes(fill= cor)) + 
  geom_tile(color="black", aes(y=reorder(ddgrid3$taxon, -ddgrid3$primarylifestyle), x=soil)) + 
  labs(y="")+
  scale_fill_gradient2(guide="colorbar",na.value="white", limits = c(-1, 1), low="orchid2", high="seagreen3", mid="white") +
  scale_x_discrete(labels=labs2) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0, hjust=0), 
        plot.margin = unit(c(0.1, 0.1, 0.3, 0), "cm"),
        axis.title.x=element_blank(), 
        #axis.text.y=element_blank(), 
        legend.position="top") +
  guides(fill=guide_colorbar(title="correlation", title.position="top", title.hjust=0.5)) +
  #geom_text(aes(x=soil, y=reorder(ddgrid3$taxon, -ddgrid3$primarylifestyle), label=pvalindic)) +
  geom_segment(inherit.aes=FALSE, size=1, data=segs2, aes(x=x, y=y, yend=yend, xend=xend)) #+
  #geom_text(inherit.aes=FALSE, data=segs2,size=3, fontface=3, aes(x=x, y=y-0.5, label=segs2$label), hjust = 1) 

p
figure <- ggpubr::ggarrange(ann2,  p, 
                            ncol = 1, nrow = 2, heights=c(0.1,0.9)) 
figure
png("figures/_lefse/1_ITS_24 hr_correlations_rarefied.png", 
    width = 1000, height =1500, res=300, units="px")
figure
dev.off()












#### 3. LEfSe: identified taxonomic groups whose abundances were 
# associated with specific treatments ** 6 mo **

# data("zeller14")
# ex_dat <- assays(zeller14)$exprs
# ex_names <- rownames(ex_dat)
# # rownames of my count data should look like this, with taxon names to the highest available specificity


# make sure count data and metadata are organized in the same order
data.frame(rownames(abundances), mdatseq$file.names)
# samples to remove
removes <- c(c("baseline1","baseline2","baseline3","baseline4"), 
             mdatseq$Sample.number[which(mdatseq$Harvest.day2=="24 hr")])
mdatseq_l <- mdatseq[-which(mdatseq$Sample.number %in% removes),]
abundances_t <- as.data.frame(t(abundances[-which(mdatseq$Sample.number %in% removes),]))

#### list row names of count data as taxa, with each row as a unique taxon
taxa2 <- taxa
taxa2$Genus <- paste0("g__", taxa2$Genus)
taxa2$Genus <- gsub(taxa2$Genus, pattern="g__NA", replacement="NA")
abundances_t$name <- paste(taxa2$Kingdom, taxa2$Phylum, taxa2$Class, taxa2$Order, taxa2$Family, 
                           taxa2$Genus, taxa2$Species, sep = "/") # , rownames(abundances_t)
abundances_t$name <- paste0(abundances_t$name, "/t__", taxa2$Primary_lifestyle) # , rownames(abundances_t)
abundances_t$name <- gsub(abundances_t$name, pattern="/NA", replacement="")
abundances_t$name <- gsub(abundances_t$name, pattern="/t_NA", replacement="")
dataname <- as.data.frame(abundances_t$name[1:50]) # check



# sum abundances by taxon
abundances_t2 <- as.data.frame(abundances_t)
abundances_t3 <- abundances_t2 %>% 
  dplyr::group_by(name) %>%
  dplyr::summarise(across(where(is.numeric), sum, na.rm = TRUE))
abundances_t3 <- as.data.frame(abundances_t3)
abundances_t3$name <- gsub(abundances_t3$name, pattern="/", replacement="|")
rownames(abundances_t3) <- abundances_t3$name
abundances_t4 <- abundances_t3[,2:(dim(abundances_t3)[2])]
abundances_t <- abundances_t4
# export abundances by taxon
write.csv(abundances_t, "raw-data/sequence/ITS/Abundances by taxon_6mo_rarefied.csv")


# dummy code residue treatment
mdatseq_l$Residue.type.nores <- as.factor(ifelse(mdatseq_l$Residue.type == 'No residue', 1, 0))
mdatseq_l$Residue.type.wheat <-  as.factor(ifelse(mdatseq_l$Residue.type == 'Wheat', 1, 0))
mdatseq_l$Residue.type.clove <-  as.factor(ifelse(mdatseq_l$Residue.type == 'Clover', 1, 0))
mdatseq_l$Residue.type.whecl <-  as.factor(ifelse(mdatseq_l$Residue.type == 'Wheat-clover mix', 1, 0))
mdatseq_l$Residue.type.fsm <-  as.factor(ifelse(mdatseq_l$Residue.type == 'Five-species mix', 1, 0))


#### lefse models by residue type for 6 mo 
sub <- which(mdatseq_l$Harvest.day2=="6 mo")
lef_6mo <- SummarizedExperiment(list(counts=abundances_t[,sub]), colData = mdatseq_l[sub,])

lef_6mo_nores2 <- lefser(lef_6mo, groupCol = "Residue.type.nores")
lefserPlot(lef_6mo_nores2, trim.names = FALSE ,colors = c(covercols[1], covercols[1]))
lef_6mo_nores2$Residue.type <- rep("No residue", dim(lef_6mo_nores2)[1])

lef_6mo_wheat2 <- lefser(lef_6mo, groupCol = "Residue.type.wheat")
lefserPlot(lef_6mo_wheat2, trim.names = FALSE ,colors = c(covercols[2], covercols[2]))
lef_6mo_wheat2$Residue.type <- rep("Wheat", dim(lef_6mo_wheat2)[1])

lef_6mo_clover2 <- lefser(lef_6mo, groupCol = "Residue.type.clove", )
lefserPlot(lef_6mo_clover2, trim.names = FALSE ,colors = c(covercols[3], covercols[3]))
lef_6mo_clover2$Residue.type <- rep("Clover", dim(lef_6mo_clover2)[1])

lef_6mo_whecl2 <- lefser(lef_6mo, groupCol = "Residue.type.whecl")
lefserPlot(lef_6mo_whecl2, trim.names = FALSE ,colors = c(covercols[4], covercols[4]))
lef_6mo_whecl2$Residue.type <- rep("Wheat-clover mix", dim(lef_6mo_whecl2)[1])

lef_6mo_fsm2 <- lefser(lef_6mo, groupCol = "Residue.type.fsm")
lefserPlot(lef_6mo_fsm2, trim.names = FALSE ,colors = c(covercols[5], covercols[5]))
lef_6mo_fsm2$Residue.type <- rep("Five-speices mix", dim(lef_6mo_fsm2)[1])

# merge all lefse results together
aa <- rbind(lef_6mo_nores2,lef_6mo_wheat2)
bb <- rbind(lef_6mo_clover2,lef_6mo_whecl2)
cc <- rbind(aa,bb)
dd <- rbind(cc,lef_6mo_fsm2)

# create columns for taxonomic group
dd$kingdom <- rep(NA, dim(dd)[1])
dd$phylum <- rep(NA, dim(dd)[1])
dd$class <- rep(NA, dim(dd)[1])
dd$order <- rep(NA, dim(dd)[1])
dd$family <- rep(NA, dim(dd)[1])
dd$genus <- rep(NA, dim(dd)[1])
dd$species <- rep(NA, dim(dd)[1])
dd$primarylifestyle <- rep(NA, dim(dd)[1])

for(i in 1:dim(dd)[1]){
  split_name <- strsplit(dd$features[i], "[|]")
  dd$kingdom[i] <- split_name[[1]][which(startsWith(split_name[[1]], "k__"))]
  dd$phylum[i] <- ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "p__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "p__"))])
  dd$class[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "c__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "c__"))])
  dd$order[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "o__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "o__"))])
  dd$family[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "f__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "f__"))])
  dd$genus[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "g__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "g__"))])
  dd$species[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "s__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "s__"))])
  dd$primarylifestyle[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "t_"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "t_"))])
}

dd$genus <- gsub(dd$genus, pattern = "g__", replacement="")
dd$species <- gsub(dd$species, pattern = "s__", replacement="")
dd$taxon <- paste(dd$genus, dd$species, sep=" ")
dd$taxon <- gsub(dd$taxon, pattern = " NA", replacement="")
dd$taxon[which(dd$taxon %in% c(" ", "NA"))] <- dd$family[which(dd$taxon %in% c(" ", "NA"))]
dd$taxon[which(dd$taxon %in% c(" ", NA))] <- dd$order[which(dd$taxon %in% c(" ", NA))]
dd$taxon[which(dd$taxon %in% c(" ", NA))] <- dd$class[which(dd$taxon %in% c(" ", NA))]
dd$taxon[which(dd$taxon %in% c(" ", NA))] <- dd$phylum[which(dd$taxon %in% c(" ", NA))]
dd$taxon[which(dd$taxon %in% c(" ", NA))] <- dd$kingdom[which(dd$taxon %in% c(" ", NA))]
dd$primarylifestyle <- gsub(dd$primarylifestyle, pattern="t__", replacement="")

dd$Residue.type <- as.factor(dd$Residue.type)
dd$Residue.type <- factor(dd$Residue.type, levels(dd$Residue.type)[c(3,4,1,5,2)])
pl <- unique(dd$primarylifestyle)

dd <- dd[order(dd$primarylifestyle, dd$taxon),]

# export lefse results 
write.csv(dd, "raw-data/sequence/ITS/Lefse_6mo_rarefied.csv")

# export  abundances for lefse taxa 
abundances_t_lefse <- abundances_t[which(rownames(abundances_t) %in% dd$features),]
write.csv(abundances_t_lefse, "raw-data/sequence/ITS/Abundances by taxon_6mo_Lefse_rarefied.csv")



# create a grid for heatmap with lefse taxa, effect sizes, and residue types
ddgrid <- expand.grid(Residue.type=levels(dd$Residue.type), taxon=dd$taxon)
ddgrid$score <- rep(NA, dim(ddgrid)[1])
ddgrid$primarylifestyle <- rep(NA, dim(ddgrid)[1])
for(i in 1:dim(ddgrid)[1]){
  score <-  dd$scores[which(dd$Residue.type==ddgrid$Residue.type[i] & dd$taxon==ddgrid$taxon[i])]
  ddgrid$score[i] <- ifelse(length(score)==0, NA, score)
  lifestyle <-  dd$primarylifestyle[which(dd$taxon==ddgrid$taxon[i])]
  ddgrid$primarylifestyle[i] <- ifelse(length(lifestyle)==0, NA, lifestyle)
}

ddgrid$primarylifestyle[which(ddgrid$primarylifestyle=="NA")] <- "unspecified"
ddgrid$primarylifestyle <- as.factor(ddgrid$primarylifestyle)
ddgrid$primarylifestyle <- factor(ddgrid$primarylifestyle, levels(ddgrid$primarylifestyle)[rev(c(2,4,6,1,3,5))]) #[rev(c(2,4,7,9,10,1,3,5,6,8))])
ddgrid3 <- ddgrid[order(ddgrid$primarylifestyle, ddgrid$taxon),]
ddgrid3$taxon <- as.factor(ddgrid3$taxon)
ddgrid3$taxon <- factor(ddgrid3$taxon, unique(ddgrid3$taxon))

result <- ddgrid3 %>% 
  dplyr::group_by(primarylifestyle) %>% 
  dplyr::summarise(count = length(unique(taxon)), .groups = 'drop')
result$hal <- result$count/2
result <- result %>% mutate(cumulative = cumsum(count))
result$cumulativehal <- result$cumulative-result$hal

lab <- paste0(gsub(result$primarylifestyle, pattern="_", replacement=" "))
segs <- data.frame(x=5.5, xend=0.5,
                   y=result$cumulative+0.5, yend=result$cumulative+0.5,
                   label=lab, #gsub(lab, pattern=" ", replacement="\n"),
                   ylab=result$cumulativehal)
segs$ybar1 <- c(0.2, segs$yend[c(1:dim(segs)[1]-1)]+0.2)
segs$ybar2 <- c(segs$y-0.2)

ann1 <- ggplot() + 
  geom_text(aes(x=0.1, y=0, label ="B) Fungal taxa at 6 mo"), size = 4, hjust = 0.5) +
  
  theme_void()

ddgrid3$font <- rep(3,dim(ddgrid3)[1])
for(i in 1:dim(ddgrid3)[1]){
  if(str_contains(ddgrid3$taxon[i], "__")==TRUE) {ddgrid3$font[i] <- 1}
}


lefsep <-ggplot(ddgrid3, aes(fill= score)) + 
  geom_tile(color="black", aes(y=reorder(ddgrid3$taxon, -ddgrid3$primarylifestyle), x=Residue.type)) + 
  labs(y="")+
  guides(fill= guide_colorbar(title.position="top", title.hjust = 0.5, title="LEfSE score")) +
  scale_fill_gradient2(low="orchid2", mid = "white", high="seagreen3", guide="colorbar",na.value="white") +
  #scale_fill_gradient2(guide="colorbar",na.value="white", limits = c(-1, 1), low="orchid2", high="seagreen3", mid="white") +
  theme(axis.text.x = element_text(angle = 270, vjust = 0, hjust=0), 
        axis.text.y = element_text(size=7, face="italic"), 
        plot.margin = unit(c(0.1, 0.1, 0.1, 2.1), "cm"),
        axis.title.x=element_blank(), legend.position="top") +
  geom_segment(inherit.aes=FALSE, size=1, data=segs, aes(x=x, y=y, yend=yend, xend=xend)) +
  coord_cartesian(xlim=c(1,5), clip="off")+
  geom_segment(inherit.aes=FALSE, size=1, data=segs, aes(x=-5.75, xend=-5.75, y=ybar1, yend=ybar2)) +
  annotate("text", x = -6, y = segs$ylab+0.5, label = segs$label, size=3, hjust=1, vjust=0.5) #+
#theme_minimal()

#geom_text(inherit.aes=FALSE, data=segs,size=3, fontface=3, aes(x=x, y=y-0.5, label=segs$label), hjust = 1) +

lefsep

figure <- ggpubr::ggarrange(ann1, lefsep, legend="top",
                            ncol = 1, nrow = 2, heights=c(0.1,0.9)) 
figure
png("figures/_lefse/1_ITS_6 mo_rarefied.png", 
    width = 1000, height =1500, res=300, units="px")
figure
dev.off()



#### 4. Pearson correlations: relationships between soil properties 
#### from dbRDA and soil taxa from LEfSe **6 mo**
# Table 3: Relationships between abundances of treatment-responsive microbial taxa identified 
# with LEfSe algorithm and soil C-cycling properties identified with redundancy analysis. Results 
# are shown for two timepoints: 24 hours and six months after glucose addition. Values represent 
# Pearson correlation coefficients, and those with p-values < 0.05 are colored blue (positive) or 
# red (negative). 

# import  abundances for lefse taxa 
abundances_t_lefse <- read.csv("raw-data/sequence/ITS/Abundances by taxon_6mo_Lefse_rarefied.csv")
# order by primary lifestyle
# create columns for taxonomic group
abundances_t_lefse$kingdom <- rep(NA, dim(abundances_t_lefse)[1])
abundances_t_lefse$phylum <- rep(NA, dim(abundances_t_lefse)[1])
abundances_t_lefse$class <- rep(NA, dim(abundances_t_lefse)[1])
abundances_t_lefse$order <- rep(NA, dim(abundances_t_lefse)[1])
abundances_t_lefse$family <- rep(NA, dim(abundances_t_lefse)[1])
abundances_t_lefse$genus <- rep(NA, dim(abundances_t_lefse)[1])
abundances_t_lefse$species <- rep(NA, dim(abundances_t_lefse)[1])
abundances_t_lefse$primarylifestyle <- rep(NA, dim(abundances_t_lefse)[1])
for(i in 1:dim(abundances_t_lefse)[1]){
  split_name <- strsplit(abundances_t_lefse$X[i], "[|]")
  abundances_t_lefse$kingdom[i] <- split_name[[1]][which(startsWith(split_name[[1]], "k__"))]
  abundances_t_lefse$phylum[i] <- ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "p__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "p__"))])
  abundances_t_lefse$class[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "c__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "c__"))])
  abundances_t_lefse$order[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "o__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "o__"))])
  abundances_t_lefse$family[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "f__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "f__"))])
  abundances_t_lefse$genus[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "g__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "g__"))])
  abundances_t_lefse$species[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "s__"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "s__"))])
  abundances_t_lefse$primarylifestyle[i] <-  ifelse(length(split_name[[1]][which(startsWith(split_name[[1]], "t_"))])==0, NA, split_name[[1]][which(startsWith(split_name[[1]], "t_"))])
}
abundances_t_lefse$genus <- gsub(abundances_t_lefse$genus, pattern = "g__", replacement="")
abundances_t_lefse$species <- gsub(abundances_t_lefse$species, pattern = "s__", replacement="")
abundances_t_lefse$taxon <- paste(abundances_t_lefse$genus, abundances_t_lefse$species, sep=" ")
abundances_t_lefse$taxon <- gsub(abundances_t_lefse$taxon, pattern = " NA", replacement="")
abundances_t_lefse$taxon[which(abundances_t_lefse$taxon %in% c(" ", "NA"))] <- abundances_t_lefse$family[which(abundances_t_lefse$taxon %in% c(" ", "NA"))]
abundances_t_lefse$taxon[which(abundances_t_lefse$taxon %in% c(" ", NA))] <- abundances_t_lefse$order[which(abundances_t_lefse$taxon %in% c(" ", NA))]
abundances_t_lefse$taxon[which(abundances_t_lefse$taxon %in% c(" ", NA))] <- abundances_t_lefse$class[which(abundances_t_lefse$taxon %in% c(" ", NA))]
abundances_t_lefse$taxon[which(abundances_t_lefse$taxon %in% c(" ", NA))] <- abundances_t_lefse$phylum[which(abundances_t_lefse$taxon %in% c(" ", NA))]
abundances_t_lefse$taxon[which(abundances_t_lefse$taxon %in% c(" ", NA))] <- abundances_t_lefse$kingdom[which(abundances_t_lefse$taxon %in% c(" ", NA))]
abundances_t_lefse$primarylifestyle <- gsub(abundances_t_lefse$primarylifestyle, pattern="t__", replacement="")
pl <- unique(abundances_t_lefse$primarylifestyle)
abundances_t_lefse <- abundances_t_lefse[order(abundances_t_lefse$primarylifestyle, abundances_t_lefse$taxon),]



# import  soil property data
mdatseq <- read.csv("raw-data/master-dataset-13c_ITS.csv")
mdatseq$Residue.type <- as.factor(mdatseq$Residue.type)
mdatseq$Residue.type <- factor(mdatseq$Residue.type, levels(mdatseq$Residue.type)[rev(c(2,5,1,4,3))])
mdatseq$Moisture.treatment <- factor(mdatseq$Moisture.treatment)
mdatseq$Harvest.day <- as.factor(mdatseq$Harvest.day)
mdatseq$Harvest.day2 <- mdatseq$Harvest.day
levels(mdatseq$Harvest.day2) <- c("24 hr", "6 mo")
mdatseq$Harvest.day2[which(is.na(mdatseq$Harvest.day2))] <- NA
# duplicate 13C measurements from g-50 samples to paired g-00 samples
mdatseq <- mdatseq[order(mdatseq$Glucose.addition.treatment),]
mdatseq[which(mdatseq$Glucose.addition.treatment=="G-00"),c("X13soc", "X13poc", "X13maoc", "X13doc", "X13mbc", "X13co2", "CUE", "co2.native")] <- mdatseq[which(mdatseq$Glucose.addition.treatment=="G-50"),c("X13soc", "X13poc", "X13maoc", "X13doc", "X13mbc", "X13co2", "CUE", "co2.native")]
# order by sequence file names
mdatseq <- mdatseq[order(mdatseq$file.names),]
mdatseq <- mdatseq[-which(is.na(mdatseq$richness)==TRUE),] # sample number 93, ambient moisture, wheat-clover residue, glucose addition, 6 mo harvest
mdatseq$X13residue <- mdatseq$X13soc - mdatseq$X13mbc
mdatseq$CSE <- (mdatseq$X13residue + mdatseq$X13mbc)/(mdatseq$X13residue + mdatseq$X13mbc + mdatseq$X13co2)
mdatseq$CSE[which(mdatseq$CSE<0)] <- NA
mdatseq$X13residue[which(mdatseq$X13residue<0)] <- NA

# subset to timepoint
mdatseq_6mo <- mdatseq[which(mdatseq$Harvest.day2=="6 mo"),]

abundances_t_lefse_t <- as.data.frame(t(abundances_t_lefse))
colnames(abundances_t_lefse_t) <- abundances_t_lefse_t[dim(abundances_t_lefse_t)[1],]
abundances_t_lefse_t2 <- abundances_t_lefse_t[rownames(abundances_t_lefse_t) %in% mdatseq_6mo$file.names,]
abundances_t_lefse_t2$file.names <- rownames(abundances_t_lefse_t2)
unique(mdatseq_6mo$file.names==abundances_t_lefse_t2$file.names) #check

# combine abundances and soil property data
newdat <- merge(mdatseq_6mo, abundances_t_lefse_t2, by="file.names",
                all.x = TRUE, all.y = TRUE)

terms.soil <- c(rownames(summary(myrda_6mo)$biplot)) 
terms.taxa <- colnames(abundances_t_lefse_t2)
terms.taxa <- terms.taxa[-which(terms.taxa %in% c("file.names"))]  



#newdat <- newdat[order(newdat$primarylifestyle, newdat$taxon),]


# create a grid for heatmap with lefse taxa, effect sizes, and residue types
ddgrid2 <- expand.grid(soil=colnames(newdat)[colnames(newdat) %in% terms.soil], taxon=colnames(newdat)[colnames(newdat) %in% terms.taxa])
ddgrid2$cor <- rep(NA, dim(ddgrid2)[1])
ddgrid2$pval <- rep(NA, dim(ddgrid2)[1])
ddgrid2$pvalindic <- rep(NA, dim(ddgrid2)[1])
ddgrid2$primarylifestyle <- rep(NA, dim(ddgrid2)[1])
for(i in 1:dim(ddgrid2)[1]){
  mod <- cor.test(newdat[,which(colnames(newdat)==ddgrid2$soil[i])], as.numeric(newdat[,which(colnames(newdat)==ddgrid2$taxon[i])]))
  ddgrid2$cor[i] <-  mod$estimate
  ddgrid2$pval[i] <- mod$p.value
  if(mod$p.value>0.05/length(ddgrid2$cor)){ddgrid2$cor[i] <- NA}
  if(mod$p.value<0.05/length(ddgrid2$cor)){ddgrid2$pvalindic[i] <- "*"}
  ddgrid2$significance.level[i] <- 0.05/length(ddgrid2$cor)  # this is the bon ferroni adjusted pvalue cutoff
  lifestyle <-  dd$primarylifestyle[which(dd$taxon==ddgrid2$taxon[i])]
  ddgrid2$primarylifestyle[i] <- ifelse(length(lifestyle)==0, NA, lifestyle)
}


ddgrid2$primarylifestyle[which(ddgrid2$primarylifestyle=="NA")] <- "unspecified"
ddgrid2$primarylifestyle <- as.factor(ddgrid2$primarylifestyle)
ddgrid2$primarylifestyle <- factor(ddgrid2$primarylifestyle, levels(ddgrid2$primarylifestyle)[rev(c(2,4,6,1,3,5))])
ddgrid3 <- ddgrid2[order(ddgrid2$primarylifestyle, ddgrid2$taxon),]
ddgrid3$taxon <- as.factor(ddgrid3$taxon)
ddgrid3$taxon <- factor(ddgrid3$taxon, unique(ddgrid3$taxon))

result <- ddgrid3 %>% 
  dplyr::group_by(primarylifestyle) %>% 
  dplyr::summarise(count = length(unique(taxon)), .groups = 'drop')
result <- result %>% mutate(cumulative = cumsum(count))

segs2 <- data.frame(x=length(unique(ddgrid3$soil))+0.5, xend=0,
                    y=result$cumulative+0.5, yend=result$cumulative+0.5,
                    label=paste0(gsub(result$primarylifestyle, pattern="_", replacement=" ")))

ann2 <- ggplot() + 
  geom_text(aes(x=0.1, y=0, label =""), size = 4, hjust = 0.5) +
  theme_void()
labs2 <- c(#expression(paste("CSE")),
  #expression(paste("PO"^"13","C")),
  #expression(paste("POM C:N")), 
  expression(paste("MAOM C:N")), 
  expression(paste("SOC")),  
  expression(paste("MAO"^"13","C")), 
  #expression(paste("DOC")),  
  #expression(paste(" "^"13","CO"[2])), 
  #expression(paste("DO"^"13","C")), 
  #expression(paste("SO"^"13","C")), 
  #expression(paste("Priming")),
  #expression(paste("POM-C")), 
  #expression(paste("MAOM-C")), 
  #expression(paste("DOC")), 
  expression(paste("Cumulative CO"[2])))

p<-ggplot(ddgrid3, aes(fill= cor)) + 
  geom_tile(color="black", aes(y=reorder(ddgrid3$taxon, -ddgrid3$primarylifestyle), x=soil)) + 
  labs(y="")+
  scale_fill_gradient2(guide="colorbar",na.value="white", limits = c(-1, 1), low="orchid2", high="seagreen3", mid="white") +
  scale_x_discrete(labels=labs2) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0, hjust=0), 
        plot.margin = unit(c(0.1, 0.1, 0.4, 0), "cm"),
        axis.title.x=element_blank(), 
        #axis.text.y=element_blank(), 
        legend.position="top") +
  guides(fill=guide_colorbar(title="correlation", title.position="top", title.hjust=0.5)) +
  #geom_text(aes(x=soil, y=reorder(ddgrid3$taxon, -ddgrid3$primarylifestyle), label=pvalindic)) +
  geom_segment(inherit.aes=FALSE, size=1, data=segs2, aes(x=x, y=y, yend=yend, xend=xend)) #+
#geom_text(inherit.aes=FALSE, data=segs2,size=3, fontface=3, aes(x=x, y=y-0.5, label=segs2$label), hjust = 1) 

p
figure <- ggpubr::ggarrange(ann2,  p, 
                            ncol = 1, nrow = 2, heights=c(0.1,0.9)) 
figure
png("figures/_lefse/1_ITS_6 mo_correlations_rarefied.png", 
    width = 1000, height =1500, res=300, units="px")
figure
dev.off()

























