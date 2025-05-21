library(dplyr)
library(ggplot2)
library(nlme)
library(emmeans)
library(multcomp)
library(vegan)
library(ape)
library(ggpubr)



covercols <- c("chocolate", "goldenrod", "forestgreen", "darkcyan", "darkorchid")
mcols <- c("blue", "red")


# All statistical analyses were implemented in R statistical software version 4.3.3 (R Core Team 2022). 
# Diversity indices including richness, evenness, and Simpson’s index were 
# calculated using the R package vegan (Oksanen et al. 2013). 

############## skip to line 160 if diversity indices already calculated

# asv taxonomy
taxa <- read.csv(paste0("raw-data/sequence/ITS/taxonomy/5_taxonomy-assignment-with-ft.csv"))
taxa$asv_number <- as.numeric(taxa$X.1)

# asv abundances
abundances <- readRDS(paste0("raw-data/sequence/ITS/seqtab/5seqtab.rds"))
abundances <- as.data.frame(abundances)


# species abundance curve
colnames(abundances) <- taxa$asv_number
#end <- dim(abundances)[2]
#abunds <- abundances[,c(3:end)]
#rownames(abunds) <- abundances$X.1
abunds2 <- t(abundances)

quantile(colSums(abunds2))/2000
sum(colSums(abunds2))
hist(colSums(abunds2))

# # create curves
# rcp <- rarecurve(abunds2, step = 500,  cex=0.5, tidy = TRUE)
# df <- rcp %>% dplyr::group_by(Site) %>%  dplyr::summarise(maxSample = max(Sample, na.rm=TRUE), maxSpecies = max(Species, na.rm=TRUE))
# df <- df[order(df$maxSpecies),]
# df$Site <- substr(df$Site, start=1, stop=3)
# rcpp <- ggplot(data=rcp, aes(x=Sample, y=Species, group=Site)) + lims(x=c(0,10000)) +
#   geom_line(color="blue") + 
#   theme_bw() + labs(x="Number of samples", y="Number of ASVs", title="")# + 
#   #geom_label(inherit.aes = F,data = as.data.frame(df), fill=NA,aes(label = Site, x=maxSample, y=maxSpecies), na.rm = TRUE, nudge_x = 1000, hjust = 0, size=2, label.size = 0) 
# 
# ggexport(ggarrange(rcpp, ncol=1, nrow=1), height=4500, width=4500, res=400,
#          filename = paste0("Figures/*sequence/ITS_asv abundance curve_rarefied.png"))
# 




#####################################################

# rarefy to this number of reads per sample
rar <- 7798


#### rarefy to rar reads
sampsize <- colSums(abunds2) # gives the number of sequences for each plot
sort(sampsize)[1:20]
sort(sampsize, decreasing = T)[1:20]
remove <- c(sampsize[which(sampsize<rar)]) 
#write.csv(data.frame(sample=names(remove), seqs=remove), paste0("Raw-data/sequence/",sampdate,"/",name,"/seqtab/6removed.csv"))
# abunds3 <- abunds2
# abunds3 <- abunds2[which(sampsize>rar),]
# abunds3 <- abunds3[-which(rownames(abunds3))),]

sampsize <- colSums(abunds2) # gives the number of sequences for each plot
sort(sampsize)[1:20] # check that sample sizes are greater than rar

abunds4 <- rrarefy(t(abunds2), rar)

abunds5 <- as.data.frame(t(abunds4))
abunds5$asv_number <- c(1:dim(abunds5)[1])
min(rowSums(abunds5[,1:length(sampsize)])) # some ASVs now have <10 reads
# remove ASVs that have less than 10 reads
abunds5 <- abunds5[-which(rowSums(abunds5[,1:length(sampsize)])<10),]


# merge taxa and abundance data into one df
taxabund <- merge(taxa, abunds5, 
                  by.x = "asv_number", by.y = "asv_number", 
                  all.x = FALSE, all.y = TRUE) 
rownames(taxabund) <- taxabund$asv_number

write.csv(taxabund, "raw-data/sequence/ITS/taxabund_rarefied.csv")

# read in the taxon data
taxabund <- read.csv("raw-data/sequence/ITS/taxabund_rarefied.csv")
# number of taxa
write.csv(dim(taxabund)[1], "raw-data/sequence/ITS/taxonomy/3_number-taxa.csv")
# number of phyla
tab <- as.data.frame(table(taxabund$Phylum))
tab$pct <- round(tab$Freq/sum(tab$Freq), 2)
write.csv(tab[order(tab$pct, decreasing = T),], "raw-data/sequence/ITS/taxonomy/4a_phyla.csv")
# number of classes
tab <- as.data.frame(table(taxabund$Class))
tab$pct <- round(tab$Freq/sum(tab$Freq), 2)
write.csv(tab[order(tab$pct, decreasing = T),], "raw-data/sequence/ITS/taxonomy/4b_classes.csv")
# number of orders
tab <- as.data.frame(table(taxabund$Order))
tab$pct <- round(tab$Freq/sum(tab$Freq), 2)
write.csv(tab[order(tab$pct, decreasing = T),], "raw-data/sequence/ITS/taxonomy/4c_orders.csv")
# number of families
tab <- as.data.frame(table(taxabund$Family))
tab$pct <- round(tab$Freq/sum(tab$Freq), 2)
write.csv(tab[order(tab$pct, decreasing = T),], "raw-data/sequence/ITS/taxonomy/4d_families.csv")
# number of genera
tab <- as.data.frame(table(taxabund$Genus))
tab$pct <- round(tab$Freq/sum(tab$Freq), 2)
write.csv(tab[order(tab$pct, decreasing = T),], "raw-data/sequence/ITS/taxonomy/4e_genera.csv")
# number of species
taxabund$Genus.Species <- paste(taxabund$Genus, taxabund$Species)
tab <- as.data.frame(table(taxabund$Genus.Species))
tab$pct <- round(tab$Freq/sum(tab$Freq), 2)
write.csv(tab[order(tab$pct, decreasing = T),], "raw-data/sequence/ITS/taxonomy/4f_species.csv")



abundances <- t(taxabund[,14:176])
#####################################################






# divmat_pa: matrix with ASVs as columns and plots as rows (proportional abundance)
# Transform to get proportional abundance (divide by total number of reads within each sample)
divmat_pa <- proportions(as.matrix(abundances),1)
unique(na.omit(rowSums(divmat_pa))) # check to see that all rows sum up to 1
divmat_pa <- data.frame(divmat_pa)
divmat_pa[1:10,1:10]

# divmat_presabs: matrix with ASVs as columns and plots as rows (presence/absence)
divmat_presabs <- divmat_pa
divmat_presabs <- ifelse(divmat_presabs > 0, 1, 0)
range(na.omit(rowSums(divmat_presabs))) # check to see that all rows sum up to integers
abundances[1:10,1:10] 
divmat_presabs[1:10,1:10] # compare to original divmat matrix


# look at summary of reads for each sample
summary(rowSums(abundances))
summary(rowSums(divmat_pa))
summary(rowSums(divmat_presabs))


# calculate diversity indices

diversity_table <- data.frame(sample = rownames(abundances))

########## richness
### (total number of species)
richness <- rowSums(divmat_presabs)
diversity_table$richness <- richness
hist(diversity_table$richness)

########## Shannon:
### (-sum p_i log(b) p_i, where p_i is the proportional abundance of species i and b is the base of the logarithm)
shannon.div <- vegan::diversity(divmat_pa, "shannon")
diversity_table$shannon.div <- shannon.div
hist(diversity_table$shannon.div)

########## Simpson's:
### (1 - sum p_i^2)
simpson.div <- vegan::diversity(divmat_pa, "simpson")
diversity_table$simpson.div <- simpson.div
hist(diversity_table$simpson.div)

########## Inverse Simpson's:
### (1/sum p_i^2)
invsimpson.div <- vegan::diversity(divmat_pa, "invsimpson")
diversity_table$invsimpson.div <- invsimpson.div
hist(diversity_table$invsimpson.div)

########## Evenness:
# ***** must give this function count data, not proportional abundance
### (H/ln(total number of species))
evenness <- microbiome::evenness(t(abundances), "pielou")
diversity_table$evenness <- evenness$pielou
hist(diversity_table$evenness)



# combine diversity data with metadata and export
metadat <- read.csv("raw-data/sequence/metadat.csv")

diversity_dat <- merge(metadat, diversity_table, 
                  by.x = "file.names", by.y = "sample", 
                  all.x = TRUE, all.y = TRUE) 
diversity_dat_ITS <- diversity_dat[which(diversity_dat$library=="ITS"),]

write.csv(diversity_dat_ITS, "raw-data/sequence/ITS-diversity_rarefied.csv")

diversity_dat_ITS2 <- diversity_dat_ITS[,c(1,12,19:23)]


# combine diversity data with master incubation data 
dat <- read.csv("raw-data/master-dataset-13c.csv")

dat <- merge(dat, diversity_dat_ITS2, 
                       by.x = "Sample.number", by.y = "Sample.number", 
                       all.x = TRUE, all.y = TRUE) 
# set groupings
dat$Residue.type <- as.factor(dat$Residue.type)
dat$Residue.type <- factor(dat$Residue.type, levels(dat$Residue.type)[c(2,5,1,4,3)])
levels(dat$Residue.type) <- c("No residue", "Wheat", "Clover", "Wheat-clover mix", "Five-species mix")
dat$Moisture.treatment <- factor(dat$Moisture.treatment)
dat$Harvest.day <- as.factor(dat$Harvest.day)
dat$Harvest.day2 <- dat$Harvest.day
levels(dat$Harvest.day2) <- c("24 hr", "6 mo")

dat2 <- dat[order(dat$Glucose.addition.treatment, dat$Sample.number),]
dat3 <- rbind(dat2[order(as.numeric(dat2$Sample.number[1:160])),], dat2[161:164,])

# export
write.csv(dat3, "raw-data/master-dataset-13c_ITS_rarefied.csv")










####################### diversity models

dat <- read.csv("raw-data/master-dataset-13c_ITS_rarefied.csv")
# set groupings
dat$Residue.type <- as.factor(dat$Residue.type)
dat$Residue.type <- factor(dat$Residue.type, levels(dat$Residue.type)[rev(c(2,5,1,4,3))])
#levels(dat$Residue.type) <- c("No residue", "Wheat", "Clover", "Wheat-clover mix", "Five-species mix")
dat$Moisture.treatment <- factor(dat$Moisture.treatment)
dat$Harvest.day <- as.factor(dat$Harvest.day)
dat$Harvest.day2 <- dat$Harvest.day
levels(dat$Harvest.day2) <- c("24 hr", "6 mo")


# remove baseline values
bdat <- dat[which(dat$Sample.number %in% c("baseline1", "baseline2", "baseline3", "baseline4")),]
dat <- dat[-which(dat$Sample.number %in% c("baseline1", "baseline2", "baseline3", "baseline4")),]



# Response variables were subject to mixed-effects ANOVAs implemented with the nlme 
# package (Pinheiro et al. 2022) with replicate as a random effect to determine residue type and 
# moisture treatment effects: diversity indices of soil fungal and bacterial communities. 
# Models were implemented separately for 
# each timepoint. Where necessary, response variables were transformed to improve residual normality. 
# Fixed effects with p<0.05 were considered significant. Differences among treatments were determined 
# from Bonferroni-corrected Tukey’s post-hoc tests using the emmeans package (Lenth 2022). 

# model: lme
# Main effects: response ~ Moisture.treatment*Residue.type, 
# random: random=~1|Lab.rep
# weights=varIdent(form=~1|Moisture.treatment*Residue.type),
# data=dat[which(dat$Harvest.day2=="timepoint"),]
# na.action="na.omit"

##### anova output dataframe ##### 
p<-10
mod.anova <- data.frame(response=rep(NA,p),
                        shapiro=rep(NA,p),
                        var=rep(NA,p),
                        Moisture.treatment=rep(NA,p), 
                        Residue.type=rep(NA,p), 
                        #Harvest.day2=rep(NA,p), 
                        Moisture.treatment_Residue.type=rep(NA,p))#, 
#Moisture.treatment_Harvest.day2=rep(NA,p), 
# Residue.type_Harvest.day2=rep(NA,p), 
#Moisture.treatment_Residue.type_Harvest.day2=rep(NA,p))
mod.anova2 <- t(mod.anova)


dat %>%
  group_by(Harvest.day2, Moisture.treatment, Residue.type) %>%
  summarise(richness = mean(richness, na.rm = TRUE))


####### richness
n <- 1
mod.anova2[1,n] <- "ITS_richness_24 hr"
hist((dat$richness[which(dat$Harvest.day2=="24 hr")]))
#dat$richness[which(dat$Harvest.day2=="24 hr" & dat$richness<200)] <- NA
# first check if glucose influenced richness
fit0 <- lme((richness) ~ Moisture.treatment*Residue.type*Glucose.addition.treatment, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
mod.out0 <- as.data.frame(car::Anova(fit0)) # default SS for Anova is type-III
mod.out20 <- mod.out0 %>% mutate_if(is.numeric, round, digits=3)
write.csv(mod.out20, "model-output/_uv/anova_3way_ITS_richness_24 hr_rarefied.csv")
# if glucose not significant, remove from model:
fit <- lme((richness) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$richness)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(richness),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
 mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_richness_24 hr_2-way_rarefied.csv")
testlet_richness_24hr <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_richness_24 hr_2-way_mois_rarefied.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_richness_24 hr_2-way_res_rarefied.csv")
testlet_richness_24hr_res <- testlet

n <- 2
mod.anova2[1,n] <- "ITS_richness_6 mo"
hist(dat$richness[which(dat$Harvest.day2=="6 mo")])
# first check if glucose influenced richness
fit0 <- lme((richness) ~ Moisture.treatment*Residue.type*Glucose.addition.treatment, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
mod.out0 <- as.data.frame(car::Anova(fit0)) # default SS for Anova is type-III
mod.out20 <- mod.out0 %>% mutate_if(is.numeric, round, digits=3)
write.csv(mod.out20, "model-output/_uv/anova_3way_ITS_richness_6 mo_rarefied.csv")
# if glucose not significant, remove from model:
fit <- lme(log(richness) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$richness)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(richness),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
 mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_richness_6 mo_2-way_rarefied.csv")
testlet_richness_6mo <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_richness_6 mo_2-way_mois_rarefied.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_richness_6 mo_2-way_res_rarefied.csv")
testlet_richness_6mo_res <- testlet


testlet_richness <- rbind(testlet_richness_24hr, testlet_richness_6mo)
testlet_richness$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_richness$minys <- testlet_richness$response-testlet_richness$SE
#testlet_richness$minys[which(testlet_richness$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_richness, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("Fungal richness"^""))) +   
  facet_wrap(.~Harvest.day2, scales="free_y") + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  #ylim(c(0,max(testlet_richness$response+testlet_richness$SE)+0.2*max(testlet_richness$response+testlet_richness$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=c(testlet_richness$.group), y=c(rep(35,10), rep(30,10))+c(testlet_richness$response)+c(testlet_richness$SE)), 
            position = position_dodge(0.9), size=2)
p

# The two lines we want on the plot
fig_richness3 <- p
# export figure
ggpubr::ggexport(fig_richness3, height=1300, width=2600, filename = "figures/_uv/_diversity/1_ITS_richness_2way_rarefied.png", res = 400)


testlet_richness_res <- rbind(testlet_richness_24hr_res, testlet_richness_6mo_res)
testlet_richness_res$Harvest.day2 <- as.factor(c(rep("24 hr", 5), rep("6 mo", 5)))
testlet_richness_res$minys <- testlet_richness_res$response-testlet_richness_res$SE
#testlet_richness$minys[which(testlet_richness$minys<0)] <- 0


# plot
segs <- data.frame(x=rep(1:5,2)-0.3, xend=rep(1:5,2)+0.3,
                   y=testlet_richness_res$response+testlet_richness_res$SE+40, yend=testlet_richness_res$response+testlet_richness_res$SE+40,
                   Harvest.day2=c(rep("24 hr", 5), rep("6 mo", 5)))
segs$Harvest.day2 <- as.factor(segs$Harvest.day2)


ann1 <- ggplot() + 
  geom_text(aes(x=0.1, y=0, label ="A) Fungal richness at 24 hr"), size = 4, hjust = 0.5) +
  geom_text(aes(x=1.1, y=0, label ="B) Fungal richness at 6 mo"), size = 4, hjust = 0.5) +
  lims(x=c(-0.5,1.5)) +
  theme_void()

p <- ggplot(data=testlet_richness, aes(x=Residue.type, y=response, fill=Moisture.treatment)) +
  theme_minimal() + labs(fill="Moisture treatment") +
  ylab(expression(paste("Soil fungal richness"^""))) +   
  facet_wrap(.~Harvest.day2) + guides(fill="none") + # , scales="free_y" + guides(fill=guide_legend(nrow=1,byrow=TRUE))
  geom_bar(stat = "identity", position = position_dodge(0.8), color="black", width=0.8) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.8)) +
  scale_fill_manual(values = c("blue", "red")) +
  #ylim(c(0,max(testlet_richness$response+testlet_richness$SE)+0.2*max(testlet_richness$response+testlet_richness$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.1, 0.5), "cm"), axis.title.x=element_blank(),
         axis.text.x = element_blank(), # element_text(angle = 290, vjust = 0, hjust=0)
        strip.text.x = element_blank()) +
  geom_segment(inherit.aes=FALSE, data=segs, aes(x=x, y=y, yend=yend, xend=xend)) +
  geom_text(inherit.aes=FALSE, data=testlet_richness_res,  aes(label=testlet_richness_res$.group, x=rep(1:5, 2), y=testlet_richness_res$response+testlet_richness_res$SE+75))
  #           position = position_dodge(0.8), size=2) #+
  # geom_text(aes(label=c(rep(NA,10), testlet_richness$.group[11:20]), y=0.003+c(rep(0,10), testlet_richness$response[11:20])+c(rep(0,10), testlet_richness$SE[11:20])), 
  #           position = position_dodge(0.8), size=2)
p


figure <- ggpubr::ggarrange(ann1, p, legend="none",
                            ncol = 1, nrow = 2, heights=c(0.1,0.9)) 
figure
png("figures/_uv/_diversity/1_ITS_richness_2way_b_rarefied.png", 
    width = 1700, height =900, res=300, units="px")
figure
dev.off()

# # The two lines we want on the plot
# fig_richnessb <- p
# # export figure
# ggpubr::ggexport(fig_richnessb, height=1700, width=1700, filename = "figures/_uv/_diversity/1_ITS_richness_2way_b_rarefied.png", res = 400)







####### shannon.div
n <- 3
mod.anova2[1,n] <- "ITS_shannon.div_24 hr"
hist(dat$shannon.div[which(dat$Harvest.day2=="24 hr")])
# first check if glucose influenced shannon.div
fit0 <- lme((shannon.div) ~ Moisture.treatment*Residue.type*Glucose.addition.treatment, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
mod.out0 <- as.data.frame(car::Anova(fit0)) # default SS for Anova is type-III
mod.out20 <- mod.out0 %>% mutate_if(is.numeric, round, digits=3)
write.csv(mod.out20, "model-output/_uv/anova_3way_ITS_shannon.div_24 hr_rarefied.csv")
### tukey: 
test1 <- emmeans(fit0, ~Moisture.treatment*Glucose.addition.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Glucose.addition.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_shannon.div_24 hr_3-way_rarefied.csv")
# if glucose not significant, remove from model:
fit <- lme(log(shannon.div) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$shannon.div)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(shannon.div),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
 mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_shannon.div_24 hr_2-way_rarefied.csv")
testlet_shannon.div_24hr <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_shannon.div_24 hr_2-way_mois_rarefied.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGF", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_shannon.div_24 hr_2-way_res_rarefied.csv")
testlet_shannon.div_24hr_res <- testlet

n <- 4
mod.anova2[1,n] <- "ITS_shannon.div_6 mo"
hist(dat$shannon.div[which(dat$Harvest.day2=="6 mo")])
# first check if glucose influenced shannon.div
fit0 <- lme((shannon.div) ~ Moisture.treatment*Residue.type*Glucose.addition.treatment, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
mod.out0 <- as.data.frame(car::Anova(fit0)) # default SS for Anova is type-III
mod.out20 <- mod.out0 %>% mutate_if(is.numeric, round, digits=3)
write.csv(mod.out20, "model-output/_uv/anova_3way_ITS_shannon.div_6 mo_rarefied.csv")
# if glucose not significant, remove from model:
fit <- lme(log(shannon.div) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$shannon.div)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(shannon.div),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
 mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_shannon.div_6 mo_2-way_rarefied.csv")
testlet_shannon.div_6mo <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_shannon.div_6 mo_2-way_mois_rarefied.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGF", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_shannon.div_6 mo_2-way_res_rarefied.csv")
testlet_shannon.div_6mo_res <- testlet


testlet_shannon.div <- rbind(testlet_shannon.div_24hr, testlet_shannon.div_6mo)
testlet_shannon.div$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_shannon.div$minys <- testlet_shannon.div$response-testlet_shannon.div$SE
#testlet_shannon.div$minys[which(testlet_shannon.div$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_shannon.div, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("Fungal shannon.div"^""))) +   
  facet_wrap(.~Harvest.day2) + #, scales="free_y"
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  #ylim(c(0,max(testlet_shannon.div$response+testlet_shannon.div$SE)+0.2*max(testlet_shannon.div$response+testlet_shannon.div$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=c(testlet_shannon.div$.group), y=c(rep(0.25,10), rep(0.25,10))+c(testlet_shannon.div$response)+c(testlet_shannon.div$SE)), 
            position = position_dodge(0.9), size=2)
p

# The two lines we want on the plot
fig_shannon.div3 <- p
# export figure
ggpubr::ggexport(fig_shannon.div3, height=1300, width=2600, filename = "figures/_uv/_diversity/1_ITS_shannon.div_2way_rarefied.png", res = 400)


testlet_shannon.div_res <- rbind(testlet_shannon.div_24hr_res, testlet_shannon.div_6mo_res)
testlet_shannon.div_res$Harvest.day2 <- as.factor(c(rep("24 hr", 5), rep("6 mo", 5)))
testlet_shannon.div_res$minys <- testlet_shannon.div_res$response-testlet_shannon.div_res$SE
#testlet_shannon.div$minys[which(testlet_shannon.div$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_shannon.div, aes(x=Residue.type, y=response, fill=Moisture.treatment)) +
  theme_minimal() + labs(fill="Moisture treatment") +
  ylab(expression(paste("Soil fungal shannon.div"^""))) +   
  facet_wrap(.~Harvest.day2, scales="free_y") + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), color="black", width=0.8) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.8)) +
  scale_fill_manual(values = c("blue", "red")) +
  #ylim(c(0,max(testlet_shannon.div$response+testlet_shannon.div$SE)+0.2*max(testlet_shannon.div$response+testlet_shannon.div$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.1, 0.5), "cm"), axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0)) +
  geom_text(inherit.aes=FALSE, data=testlet_shannon.div_res,  aes(label=testlet_shannon.div_res$.group, x=rep(1:5, 2), y=testlet_shannon.div_res$response+testlet_shannon.div_res$SE+0.5))
#           position = position_dodge(0.8), size=2) #+
# geom_text(aes(label=c(rep(NA,10), testlet_shannon.div$.group[11:20]), y=0.003+c(rep(0,10), testlet_shannon.div$response[11:20])+c(rep(0,10), testlet_shannon.div$SE[11:20])), 
#           position = position_dodge(0.8), size=2)
p

# The two lines we want on the plot
fig_shannon.divb <- p
# export figure
ggpubr::ggexport(fig_shannon.divb, height=1700, width=2700, filename = "figures/_uv/_diversity/1_ITS_shannon.div_2way_b_rarefied.png", res = 400)









####### simpson.div
n <- 5
mod.anova2[1,n] <- "ITS_simpson.div_24 hr"
hist(dat$simpson.div[which(dat$Harvest.day2=="24 hr")])
# first check if glucose influenced simpson.div
fit0 <- lme((simpson.div) ~ Moisture.treatment*Residue.type*Glucose.addition.treatment, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
mod.out0 <- as.data.frame(car::Anova(fit0)) # default SS for Anova is type-III
mod.out20 <- mod.out0 %>% mutate_if(is.numeric, round, digits=3)
write.csv(mod.out20, "model-output/_uv/anova_3way_ITS_simpson.div_24 hr_rarefied.csv")
# if glucose not significant, remove from model:
fit <- lme((simpson.div) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$simpson.div)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(simpson.div),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
 mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_simpson.div_24 hr_2-way_rarefied.csv")
testlet_simpson.div_24hr <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_simpson.div_24 hr_2-way_mois_rarefied.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGF", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_simpson.div_24 hr_2-way_res_rarefied.csv")
testlet_simpson.div_24hr_res <- testlet

n <- 6
mod.anova2[1,n] <- "ITS_simpson.div_6 mo"
hist(dat$simpson.div[which(dat$Harvest.day2=="6 mo")])
# first check if glucose influenced simpson.div
fit0 <- lme((simpson.div) ~ Moisture.treatment*Residue.type*Glucose.addition.treatment, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
mod.out0 <- as.data.frame(car::Anova(fit0)) # default SS for Anova is type-III
mod.out20 <- mod.out0 %>% mutate_if(is.numeric, round, digits=3)
write.csv(mod.out20, "model-output/_uv/anova_3way_ITS_simpson.div_6 mo_rarefied.csv")
### tukey: 
test1 <- emmeans(fit0, ~Moisture.treatment*Residue.type*Glucose.addition.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGFG", reversed = TRUE)
testlet <- testlet[order(testlet$Glucose.addition.treatment, testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_simpson.div_6 mo_3-way_rarefied.csv")
# if glucose not significant, remove from model:
fit <- lme((simpson.div) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$simpson.div)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(simpson.div),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
 mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_simpson.div_6 mo_2-way_rarefied.csv")
testlet_simpson.div_6mo <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_simpson.div_6 mo_2-way_mois_rarefied.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGF", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_simpson.div_6 mo_2-way_res_rarefied.csv")
testlet_simpson.div_6mo_res <- testlet


testlet_simpson.div <- rbind(testlet_simpson.div_24hr, testlet_simpson.div_6mo)
testlet_simpson.div$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_simpson.div$minys <- testlet_simpson.div$response-testlet_simpson.div$SE
#testlet_simpson.div$minys[which(testlet_simpson.div$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_simpson.div, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("Fungal simpson.div"^""))) +   
  facet_wrap(.~Harvest.day2) + #, scales="free_y"
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  #ylim(c(0,max(testlet_simpson.div$response+testlet_simpson.div$SE)+0.2*max(testlet_simpson.div$response+testlet_simpson.div$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=c(testlet_simpson.div$.group), y=c(rep(0.03,10), rep(0.03,10))+c(testlet_simpson.div$response)+c(testlet_simpson.div$SE)), 
            position = position_dodge(0.9), size=2)
p

# The two lines we want on the plot
fig_simpson.div3 <- p
# export figure
ggpubr::ggexport(fig_simpson.div3, height=1300, width=2600, filename = "figures/_uv/_diversity/1_ITS_simpson.div_2way_rarefied.png", res = 400)


testlet_simpson.div_res <- rbind(testlet_simpson.div_24hr_res, testlet_simpson.div_6mo_res)
testlet_simpson.div_res$Harvest.day2 <- as.factor(c(rep("24 hr", 5), rep("6 mo", 5)))
testlet_simpson.div_res$minys <- testlet_simpson.div_res$response-testlet_simpson.div_res$SE
testlet_simpson.div_res$.group[which(testlet_simpson.div_res$Harvest.day2=="24 hr")] <- NA
#testlet_simpson.div$minys[which(testlet_simpson.div$minys<0)] <- 0

testlet_simpson.div$.group[which(testlet_simpson.div$Harvest.day2=="6 mo")] <- NA

# plot
p <- ggplot(data=testlet_simpson.div, aes(x=Residue.type, y=response, fill=Moisture.treatment)) +
  theme_minimal() + labs(fill="Moisture treatment") +
  ylab(expression(paste("Soil fungal simpson.div"^""))) +   
  facet_wrap(.~Harvest.day2, scales="free_y") + guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), color="black", width=0.8) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.8)) +
  scale_fill_manual(values = c("blue", "red")) +
  #ylim(c(0,max(testlet_simpson.div$response+testlet_simpson.div$SE)+0.2*max(testlet_simpson.div$response+testlet_simpson.div$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.1, 0.5), "cm"), axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0)) +
  geom_text(inherit.aes=FALSE, data=testlet_simpson.div_res,  aes(label=testlet_simpson.div_res$.group, x=rep(1:5, 2), 
                                                                  y=testlet_simpson.div_res$response+testlet_simpson.div_res$SE+0.05)) +
geom_text(aes(label=testlet_simpson.div$.group, y=0.03+testlet_simpson.div$response+testlet_simpson.div$SE),
          position = position_dodge(0.8), size=2)
p

# The two lines we want on the plot
fig_simpson.divb <- p
# export figure
ggpubr::ggexport(fig_simpson.divb, height=1700, width=2700, filename = "figures/_uv/_diversity/1_ITS_simpson.div_2way_b_rarefied.png", res = 400)










####### invsimpson.div
n <- 7
mod.anova2[1,n] <- "ITS_invsimpson.div_24 hr"
hist(dat$invsimpson.div[which(dat$Harvest.day2=="24 hr")])
# first check if glucose influenced invsimpson.div
fit0 <- lme(log(invsimpson.div) ~ Moisture.treatment*Residue.type*Glucose.addition.treatment, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
mod.out0 <- as.data.frame(car::Anova(fit0)) # default SS for Anova is type-III
mod.out20 <- mod.out0 %>% mutate_if(is.numeric, round, digits=3)
write.csv(mod.out20, "model-output/_uv/anova_3way_ITS_invsimpson.div_24 hr_rarefied.csv")
### tukey: 
test1 <- emmeans(fit0, ~Glucose.addition.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGFG", reversed = TRUE)
testlet <- testlet[order(testlet$Glucose.addition.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_invsimpson.div_24 hr_3-way_rarefied.csv")
# if glucose not significant, remove from model:
fit <- lme(log(invsimpson.div) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$invsimpson.div)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(invsimpson.div),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
 mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_invsimpson.div_24 hr_2-way_rarefied.csv")
testlet_invsimpson.div_24hr <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_invsimpson.div_24 hr_2-way_mois_rarefied.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGF", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_invsimpson.div_24 hr_2-way_res_rarefied.csv")
testlet_invsimpson.div_24hr_res <- testlet

n <- 8
mod.anova2[1,n] <- "ITS_invsimpson.div_6 mo"
hist(dat$invsimpson.div[which(dat$Harvest.day2=="6 mo")])
# first check if glucose influenced invsimpson.div
fit0 <- lme(log(invsimpson.div) ~ Moisture.treatment*Residue.type*Glucose.addition.treatment, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
mod.out0 <- as.data.frame(car::Anova(fit0)) # default SS for Anova is type-III
mod.out20 <- mod.out0 %>% mutate_if(is.numeric, round, digits=3)
write.csv(mod.out20, "model-output/_uv/anova_3way_ITS_invsimpson.div_6 mo_rarefied.csv")
### tukey: 
test1 <- emmeans(fit0, ~Moisture.treatment*Residue.type*Glucose.addition.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGFGHI", reversed = TRUE)
testlet <- testlet[order(testlet$Glucose.addition.treatment, testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_invsimpson.div_6 mo_3-way_rarefied.csv")
# if glucose not significant, remove from model:
fit <- lme(log(invsimpson.div) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$invsimpson.div)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(invsimpson.div),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
 mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_invsimpson.div_6 mo_2-way_rarefied.csv")
testlet_invsimpson.div_6mo <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_invsimpson.div_6 mo_2-way_mois_rarefied.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGF", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_invsimpson.div_6 mo_2-way_res_rarefied.csv")
testlet_invsimpson.div_6mo_res <- testlet


testlet_invsimpson.div <- rbind(testlet_invsimpson.div_24hr, testlet_invsimpson.div_6mo)
testlet_invsimpson.div$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_invsimpson.div$minys <- testlet_invsimpson.div$response-testlet_invsimpson.div$SE
#testlet_invsimpson.div$minys[which(testlet_invsimpson.div$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_invsimpson.div, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("Fungal invsimpson.div"^""))) +   
  facet_wrap(.~Harvest.day2) + #, scales="free_y"
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  #ylim(c(0,max(testlet_invsimpson.div$response+testlet_invsimpson.div$SE)+0.2*max(testlet_invsimpson.div$response+testlet_invsimpson.div$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=c(testlet_invsimpson.div$.group), y=c(rep(10,10), rep(10,10))+c(testlet_invsimpson.div$response)+c(testlet_invsimpson.div$SE)), 
            position = position_dodge(0.9), size=2)
p

# The two lines we want on the plot
fig_invsimpson.div3 <- p
# export figure
ggpubr::ggexport(fig_invsimpson.div3, height=1300, width=2600, filename = "figures/_uv/_diversity/1_ITS_invsimpson.div_2way_rarefied.png", res = 400)


testlet_invsimpson.div_res <- rbind(testlet_invsimpson.div_24hr_res, testlet_invsimpson.div_6mo_res)
testlet_invsimpson.div_res$Harvest.day2 <- as.factor(c(rep("24 hr", 5), rep("6 mo", 5)))
testlet_invsimpson.div_res$minys <- testlet_invsimpson.div_res$response-testlet_invsimpson.div_res$SE
#testlet_invsimpson.div_res$.group[which(testlet_invsimpson.div_res$Harvest.day2=="24 hr")] <- NA
#testlet_invsimpson.div$minys[which(testlet_invsimpson.div$minys<0)] <- 0

#testlet_invsimpson.div$.group[which(testlet_invsimpson.div$Harvest.day2=="6 mo")] <- NA

# plot
p <- ggplot(data=testlet_invsimpson.div, aes(x=Residue.type, y=response, fill=Moisture.treatment)) +
  theme_minimal() + labs(fill="Moisture treatment") +
  ylab(expression(paste("Soil fungal invsimpson.div"^""))) +   
  facet_wrap(.~Harvest.day2) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) + #, scales="free_y"
  geom_bar(stat = "identity", position = position_dodge(0.8), color="black", width=0.8) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.8)) +
  scale_fill_manual(values = c("blue", "red")) +
  #ylim(c(0,max(testlet_invsimpson.div$response+testlet_invsimpson.div$SE)+0.2*max(testlet_invsimpson.div$response+testlet_invsimpson.div$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.1, 0.5), "cm"), axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0)) +
  geom_text(inherit.aes=FALSE, data=testlet_invsimpson.div_res,  aes(label=testlet_invsimpson.div_res$.group, x=rep(1:5, 2), 
                                                                  y=testlet_invsimpson.div_res$response+testlet_invsimpson.div_res$SE+c(40,15,15,15,15, 50,15,15,15,15)))# +
  # geom_text(aes(label=testlet_invsimpson.div$.group, y=0.03+testlet_invsimpson.div$response+testlet_invsimpson.div$SE),
  #           position = position_dodge(0.8), size=2)
p

# The two lines we want on the plot
fig_invsimpson.divb <- p
# export figure
ggpubr::ggexport(fig_invsimpson.divb, height=1700, width=2700, filename = "figures/_uv/_diversity/1_ITS_invsimpson.div_2way_b_rarefied.png", res = 400)









####### evenness
n <- 9
mod.anova2[1,n] <- "ITS_evenness_24 hr"
hist(dat$evenness[which(dat$Harvest.day2=="24 hr")])
# first check if glucose influenced evenness
fit0 <- lme((evenness) ~ Moisture.treatment*Residue.type*Glucose.addition.treatment, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
mod.out0 <- as.data.frame(car::Anova(fit0)) # default SS for Anova is type-III
mod.out20 <- mod.out0 %>% mutate_if(is.numeric, round, digits=3)
write.csv(mod.out20, "model-output/_uv/anova_3way_ITS_evenness_24 hr_rarefied.csv")
### tukey: 
test1 <- emmeans(fit0, ~Glucose.addition.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed = TRUE)
testlet <- testlet[order(testlet$Glucose.addition.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_evenness_24 hr_3-way_rarefied.csv")
# if glucose not significant, remove from model:
fit <- lme((evenness) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$evenness)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(evenness),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
 mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_evenness_24 hr_2-way_rarefied.csv")
testlet_evenness_24hr <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_evenness_24 hr_2-way_mois_rarefied.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGF", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_evenness_24 hr_2-way_res_rarefied.csv")
testlet_evenness_24hr_res <- testlet

n <- 10
mod.anova2[1,n] <- "ITS_evenness_6 mo"
hist(dat$evenness[which(dat$Harvest.day2=="6 mo")])
# first check if glucose influenced evenness
fit0 <- lme(log(evenness) ~ Moisture.treatment*Residue.type*Glucose.addition.treatment, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
mod.out0 <- as.data.frame(car::Anova(fit0)) # default SS for Anova is type-III
mod.out20 <- mod.out0 %>% mutate_if(is.numeric, round, digits=3)
write.csv(mod.out20, "model-output/_uv/anova_3way_ITS_evenness_6 mo_rarefied.csv")
### tukey: 
test1 <- emmeans(fit0, ~Glucose.addition.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed = TRUE)
testlet <- testlet[order(testlet$Glucose.addition.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_evenness_6 mo_3-way_rarefied.csv")
# if glucose not significant, remove from model:
fit <- lme(log(evenness) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$evenness)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(evenness),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
 mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_evenness_6 mo_2-way_rarefied.csv")
testlet_evenness_6mo <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDEFG", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_evenness_6 mo_2-way_mois_rarefied.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEFGF", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/_diversity/ITS_evenness_6 mo_2-way_res_rarefied.csv")
testlet_evenness_6mo_res <- testlet


testlet_evenness <- rbind(testlet_evenness_24hr, testlet_evenness_6mo)
testlet_evenness$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_evenness$minys <- testlet_evenness$response-testlet_evenness$SE
#testlet_evenness$minys[which(testlet_evenness$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_evenness, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("Fungal evenness"^""))) +   
  facet_wrap(.~Harvest.day2) + #, scales="free_y"
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  #ylim(c(0,max(testlet_evenness$response+testlet_evenness$SE)+0.2*max(testlet_evenness$response+testlet_evenness$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=c(testlet_evenness$.group), y=c(rep(0.05,10), rep(0.05,10))+c(testlet_evenness$response)+c(testlet_evenness$SE)), 
            position = position_dodge(0.9), size=2)
p

# The two lines we want on the plot
fig_evenness3 <- p
# export figure
ggpubr::ggexport(fig_evenness3, height=1300, width=2600, filename = "figures/_uv/_diversity/1_ITS_evenness_2way_rarefied.png", res = 400)


testlet_evenness_res <- rbind(testlet_evenness_24hr_res, testlet_evenness_6mo_res)
testlet_evenness_res$Harvest.day2 <- as.factor(c(rep("24 hr", 5), rep("6 mo", 5)))
testlet_evenness_res$minys <- testlet_evenness_res$response-testlet_evenness_res$SE
#testlet_evenness_res$.group[which(testlet_evenness_res$Harvest.day2=="24 hr")] <- NA
#testlet_evenness$minys[which(testlet_evenness$minys<0)] <- 0

#testlet_evenness$.group[which(testlet_evenness$Harvest.day2=="6 mo")] <- NA

# plot
p <- ggplot(data=testlet_evenness, aes(x=Residue.type, y=response, fill=Moisture.treatment)) +
  theme_minimal() + labs(fill="Moisture treatment") +
  ylab(expression(paste("Soil fungal evenness"^""))) +   
  facet_wrap(.~Harvest.day2) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) + #, scales="free_y"
  geom_bar(stat = "identity", position = position_dodge(0.8), color="black", width=0.8) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.8)) +
  scale_fill_manual(values = c("blue", "red")) +
  #ylim(c(0,max(testlet_evenness$response+testlet_evenness$SE)+0.2*max(testlet_evenness$response+testlet_evenness$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.1, 0.5), "cm"), axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0)) +
  geom_text(inherit.aes=FALSE, data=testlet_evenness_res,  aes(label=testlet_evenness_res$.group, x=rep(1:5, 2), 
                                                                     y=testlet_evenness_res$response+testlet_evenness_res$SE+0.07))# +
# geom_text(aes(label=testlet_evenness$.group, y=0.03+testlet_evenness$response+testlet_evenness$SE),
#           position = position_dodge(0.8), size=2)
p

# The two lines we want on the plot
fig_evennessb <- p
# export figure
ggpubr::ggexport(fig_evennessb, height=1700, width=2700, filename = "figures/_uv/_diversity/1_ITS_evenness_2way_b_rarefied.png", res = 400)




write.csv(mod.anova2, "model-output/_uv/anova_2way_ITS_rarefied.csv")














