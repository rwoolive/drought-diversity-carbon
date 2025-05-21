
library(ggpubr)
library(corrplot)
library(ggplot2)
library(dplyr)

# 2023-05-22: initiate phase 1 (combine soil with residues and bring to 60% WHC)
# 2023-05-29: initiate drought
# 2023-06-30: bring all soils to 60% WHC
# 2023-07-10: glucose addition to set 1
# 2023-07-11: harvest set 1 (24 hr)
# 2023-07-12: glucose addition to set 2
# 2024-01-23: harvest set 2 (6 mo)


covercols <- c("chocolate", "goldenrod", "forestgreen", "darkcyan", "darkorchid")
mcols <- c("blue", "red")

# read in data
dat <- read.csv("raw-data/master-dataset-13c_ITS_16S.csv")
# set groupings
dat$Residue.type <- as.factor(dat$Residue.type)
dat$Residue.type <- factor(dat$Residue.type, levels(dat$Residue.type)[rev(c(2,5,1,4,3))])
#levels(dat$Residue.type) <- c("No residue", "Wheat", "Clover", "Wheat-clover mix", "Five-species mix")
dat$Moisture.treatment <- factor(dat$Moisture.treatment)
dat$Harvest.day <- as.factor(dat$Harvest.day)
dat$Harvest.day2 <- dat$Harvest.day
levels(dat$Harvest.day2) <- c("24 hr", "6 mo")

# CSE calculation
dat$X13residue <- dat$X13soc - dat$X13mbc
dat$CSE <- (dat$X13residue + dat$X13mbc)/(dat$X13residue + dat$X13mbc + dat$X13co2)
dat$CSE[which(dat$CSE<0)] <- NA
dat$X13residue[which(dat$X13residue<0)] <- NA

# remove baseline values
bdat <- dat[which(dat$Sample.number %in% c("baseline1", "baseline2", "baseline3")),]
dat <- dat[-which(dat$Sample.number %in% c("baseline1", "baseline2", "baseline3")),]

# # read in data
# dat_decay <- read.csv("model-output/decay_mods.csv")
# # dat_decay <- dat_decay[,1:12]
# dat_decay$sample_number <- as.character(dat_decay$sample_number)
# dat_0 <- left_join(dat, dat_decay, by = join_by(Sample.number == sample_number))
# dat <- dat_0




# treatments: 
treat <- c("Residue.type", "Moisture.treatment", "Harvest.day2")
# 13c measurements of interest: 
c13 <- c("CSE", "X13residue", "X13maoc", "X13poc", "X13co2", "X13mbc", "X13doc", "X13soc", "co2.native")
a0 <- expression(paste("Residue-"^"13","C (", mu,"g g"^-1,") "))
a1 <- expression(paste("MAO"^"13","C (", mu,"g g"^-1,") "))
a2 <- expression(paste("PO"^"13","C (", mu,"g g"^-1,") "))
a3 <- expression(paste(" "^"13","CO"[2], " (", mu,"g g"^-1,") "))
a4 <- expression(paste("MB"^"13","C (", mu,"g g"^-1,") "))
a5 <- expression(paste("DO"^"13","C (", mu,"g g"^-1,") "))
a6 <- expression(paste("SO"^"13","C (", mu,"g g"^-1,") "))
a7 <- expression(paste("Priming (", mu,"g g"^-1,") "))
c13names_24hr <- c("CSE ", a0, a1, a2, a3, a4, a5, a6, a7)

c13names_6mo <- c("CSE ", a0, a1, a2, a3, a4, a5, a6, a7)



# decay <- c("m1_3pool", "m2_3pool", "m3_3pool", "k1_3pool", "k2_3pool", "k3_3pool")
# a7 <- expression(paste("Fast "^"13","C pool size (%)"))
# a8 <- expression(paste("Intermediate "^"13","C pool size (%)"))
# a9 <- expression(paste("Slow "^"13","C pool size (%)"))
# a10 <- expression(paste("Fast "^"13","C pool mineralization rate"))
# a11 <- expression(paste("Intermediate "^"13","C pool mineralization rate"))
# a12 <- expression(paste("Slow "^"13","C pool mineralization rate"))
# decay_names <- c(a7, a8, a9, a10, a11, a12)


# other soil properties
sp <- c("pom.c.g.kg", "pom.n.g.kg", "pom.cn", 
        "maom.c.g.kg", "maom.n.g.kg", "maom.cn", 
        "bulk.soil.C.g.kg", "bulk.soil.N.g.kg", "bulk.soil.cn", 
        "doc", "mbc", "co2","dna", "richness.x", "richness.16S")
spnames_24hr <- c(expression(paste("POM-C (g kg"^-1,") ")), expression(paste("POM-N (g kg"^-1,") ")), expression(paste("POM C:N ")), 
             expression(paste("MAOM-C (g kg"^-1,") ")), expression(paste("MAOM-N (g kg"^-1,") ")), expression(paste("MAOM C:N ")), 
             expression(paste("SOC (g kg"^-1,") ")), expression(paste("Total N (g kg"^-1,") ")), expression(paste("Soil C:N ")), 
             expression(paste("DOC (mg kg"^-1,") ")), expression(paste("MBC (mg kg"^-1,") ")), 
             expression(paste("Cumulative CO"[2]," (", mu,"g g"^-1, ") ")), expression(paste("Microbial DNA (ng ", mu,"L"^-1, ") ")),
             expression(paste("Fungal richness ")), expression(paste("Bacterial richness ")))
spnames_6mo <- c(expression(paste("POM-C (g kg"^-1,") ")), expression(paste("POM-N (g kg"^-1,") ")), expression(paste("POM C:N ")), 
                  expression(paste("MAOM-C (g kg"^-1,") ")), expression(paste("MAOM-N (g kg"^-1,") ")), expression(paste("MAOM C:N ")), 
                  expression(paste("SOC (g kg"^-1,") ")), expression(paste("Total N (g kg"^-1,") ")), expression(paste("Soil C:N ")), 
                  expression(paste("DOC (mg kg"^-1,") ")), expression(paste("MBC (mg kg"^-1,") ")), 
                  expression(paste("Cumulative CO"[2]," (", mu,"g g"^-1, ") ")), expression(paste("Microbial DNA (ng ", mu,"L"^-1, ") ")),
                 expression(paste("Fungal richness ")), expression(paste("Bacterial richness ")))



# complete list of old names
c13_24hr <- paste0(c13, "_24hr")             
c13_6mo <- paste0(c13, "_6mo")             
sp_24hr <- paste0(sp, "_24hr")             
sp_6mo <- paste0(sp, "_6mo")             
all_old_names <- c(c13_24hr, c13_6mo, sp_24hr, sp_6mo) #decay, 
all_old_names_24hr <- c(c13_24hr, sp_24hr) #decay, 
all_old_names_6mo <- c(c13_6mo, sp_6mo) #decay, 


# complete list of new names
all_new_names <- c(c13names_24hr, c13names_6mo, spnames_24hr, spnames_6mo) #decay_names, 
all_new_names_24hr <- c(c13names_24hr, spnames_24hr) #decay_names, 
all_new_names_6mo <- c(c13names_6mo, spnames_6mo) #decay_names, 

#check
data.frame(all_old_names=all_old_names, all_new_names=as.character(all_new_names))


# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))




### PLOT
#layout.show(n=2)

# correlation plot: 

# data for 24 h post glucose addition
cordat1 <- dat[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T1"),which(colnames(dat) %in% c13)]
colnames(cordat1) <- paste0(colnames(cordat1), "_24hr")
cordat2 <- dat[which(dat$Glucose.addition.treatment=="G-00" & dat$Harvest.day=="T1"),which(colnames(dat) %in% sp)]
colnames(cordat2) <- paste0(colnames(cordat2), "_24hr")
cordat3 <- cbind(cordat1, cordat2)

# data for 6 mo post glucose addition
cordat4 <- dat[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T90"),which(colnames(dat) %in% c13)]
colnames(cordat4) <- paste0(colnames(cordat4), "_6mo")
cordat5 <- dat[which(dat$Glucose.addition.treatment=="G-00" & dat$Harvest.day=="T90"),which(colnames(dat) %in% sp)]
colnames(cordat5) <- paste0(colnames(cordat5), "_6mo")
cordat6 <- cbind(cordat4, cordat5)

# change in cn ratio
cordat6$pom.cn_6mo_change <- dat$pom.cn[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T90")] - dat$pom.cn[which(dat$Glucose.addition.treatment=="G-00" & dat$Harvest.day=="T90")]
cordat6$maom.cn_6mo_change <- dat$maom.cn[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T90")] - dat$maom.cn[which(dat$Glucose.addition.treatment=="G-00" & dat$Harvest.day=="T90")]
cordat6$bulk.soil.cn_6mo_change <- dat$bulk.soil.cn[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T90")] - dat$bulk.soil.cn[which(dat$Glucose.addition.treatment=="G-00" & dat$Harvest.day=="T90")]
# change in c 
cordat6$pom.c.g.kg_6mo_change <- dat$pom.c.g.kg[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T90")] - dat$pom.c.g.kg[which(dat$Glucose.addition.treatment=="G-00" & dat$Harvest.day=="T90")]
cordat6$maom.c.g.kg_6mo_change <- dat$maom.c.g.kg[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T90")] - dat$maom.c.g.kg[which(dat$Glucose.addition.treatment=="G-00" & dat$Harvest.day=="T90")]
cordat6$bulk.soil.C.g.kg_6mo_change <- dat$bulk.soil.C.g.kg[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T90")] - dat$bulk.soil.C.g.kg[which(dat$Glucose.addition.treatment=="G-00" & dat$Harvest.day=="T90")]
# change in n
cordat6$pom.n.g.kg_6mo_change <- dat$pom.n.g.kg[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T90")] - dat$pom.n.g.kg[which(dat$Glucose.addition.treatment=="G-00" & dat$Harvest.day=="T90")]
cordat6$maom.n.g.kg_6mo_change <- dat$maom.n.g.kg[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T90")] - dat$maom.n.g.kg[which(dat$Glucose.addition.treatment=="G-00" & dat$Harvest.day=="T90")]
cordat6$bulk.soil.N.g.kg_6mo_change <- dat$bulk.soil.N.g.kg[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T90")] - dat$bulk.soil.N.g.kg[which(dat$Glucose.addition.treatment=="G-00" & dat$Harvest.day=="T90")]

# # data for decay curves
# cordat7 <- dat[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T90"),which(colnames(dat) %in% decay)]

# combine data
cordat8 <- cbind(cordat3, cordat6)
cordat9 <- cordat8 #cbind(cordat8, cordat7)



############### correlations: 24 hr

M<-cor(cordat9)

# matrix of the p-value of the correlation
p.mat <- cor.mtest(cordat9)
p.mat3<- p.mat[c(c13_24hr), all_old_names_24hr] #, decay

# reorder the matrix
m2 <- M[all_old_names_24hr, all_old_names_24hr]
m3 <- M[c(c13_24hr), all_old_names_24hr] # , decay
write.csv(m3, "model-output/correlations_24hr.csv")

# create new names for axes
labs_y_24hr <- all_new_names_24hr[match(c(c13_24hr), rownames(m3))] # , decay
data.frame(a=as.character(labs_y_24hr), b=rownames(m3)) # check

labs_x_24hr <- all_new_names_24hr[match(all_old_names_24hr, colnames(m2))]
data.frame(a=as.character(labs_x_24hr), b=rownames(m2)) # check

# # create new names for x
# xlabs <- spnames[match(sp, colnames(M))]
# data.frame(a=as.character(xlabs), b=colnames(m2)) # check


png("figures/4_correlations/4_corplot_a_24hr.png", width = 4000, height =2300, res=300, units="px")

corrplot(m3, method="color",col=col(200), addgrid.col = "gray",number.cex = 0.7, # order="hclust",
         p.mat = p.mat3, sig.level = 0.05, insig = "blank", cl.pos="r",
         addCoef.col = "black", title="", 
         tl.pos='n', mar = c(10, 15, 1, 1)) # , type = "lower"
mtext(text = "A) 24 hr", side = 3, line = -2, cex=1.5)
text(x=seq(1, dim(m3)[2], length.out=dim(m3)[2]), y=par("usr")[3]+3, labels = labs_x_24hr, srt=-55, adj=0) # x axis labels
mtext(text = labs_y_24hr, side = 2, line = -11.5, at = seq(dim(m3)[1], 1, length.out=dim(m3)[1]), las = 1, adj=1) # y axis labels
dev.off()


############### correlations: 6 mo

M<-cor(cordat9, use = "na.or.complete")

# matrix of the p-value of the correlation
p.mat <- cor.mtest(cordat9)
p.mat3<- p.mat[c(c13_6mo), all_old_names_6mo] #, decay

# reorder the matrix
m2 <- M[all_old_names_6mo, all_old_names_6mo]
m3 <- M[c(c13_6mo), all_old_names_6mo] # , decay
write.csv(m3, "model-output/correlations_6mo.csv")

# create new names for axes
labs_y_6mo <- all_new_names_6mo[match(c(c13_6mo), rownames(m3))] # , decay
data.frame(a=as.character(labs_y_6mo), b=rownames(m3)) # check

labs_x_6mo <- all_new_names_6mo[match(all_old_names_6mo, colnames(m2))]
data.frame(a=as.character(labs_x_6mo), b=rownames(m2)) # check

# # create new names for x
# xlabs <- spnames[match(sp, colnames(M))]
# data.frame(a=as.character(xlabs), b=colnames(m2)) # check


png("figures/4_correlations/4_corplot_a_6mo.png", width = 4000, height =2300, res=300, units="px")

corrplot(m3, method="color",col=col(200), addgrid.col = "gray",number.cex = 0.7, # order="hclust",
         p.mat = p.mat3, sig.level = 0.05, insig = "blank", cl.pos="r",
         addCoef.col = "black", title="", 
         tl.pos='n', mar = c(10, 15, 1, 1)) # , type = "lower"
mtext(text = "B) 6 mo", side = 3, line = -2, cex=1.5)
text(x=seq(1, dim(m3)[2], length.out=dim(m3)[2]), y=par("usr")[3]+3, labels = labs_x_6mo, srt=-55, adj=0) # x axis labels
mtext(text = labs_y_6mo, side = 2, line = -11.5, at = seq(dim(m3)[1], 1, length.out=dim(m3)[1]), las = 1, adj=1) # y axis labels
dev.off()



######### plot the strong correlations (>0.8)

cordat9$Residue.type <- dat$Residue.type[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T1")]
cordat9$Moisture.treatment <- dat$Moisture.treatment[which(dat$Glucose.addition.treatment=="G-50" & dat$Harvest.day=="T1")]




###### CN ratio v. priming

m2 <- ggplot(cordat9, aes(x = pom.cn_6mo, y = co2.native_6mo, size=0.5)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", linetype="dashed") +
  geom_point(alpha=0.7, aes(colour=Residue.type)) +
  scale_colour_manual(values=covercols)+
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) + xlim(c(15, 27)) +
  xlab(expression(paste("Particulate organic matter C:N")))+ 
  ylab(expression(paste("Priming (", mu,"g g"^"-1",")")))+ 
  theme_minimal() + 
  facet_wrap(.~Moisture.treatment) +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"), panel.spacing = unit(1, "lines")) +
  ggpubr::stat_cor(method = "pearson")+ guides(colour = guide_legend(override.aes = list(size=3)))

m2

m1 <- ggplot(cordat9, aes(x = maom.cn_6mo, y = co2.native_6mo, size=0.5)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", linetype="dashed") +
  geom_point(alpha=0.7, aes(colour=Residue.type)) +
  scale_colour_manual(values=covercols)+
  labs(colour="Residue type")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) +
  xlab(expression(paste("Mineral-associated organic matter C:N")))+ 
  ylab(expression(paste("Priming (", mu,"g g"^"-1",")")))+ 
  theme_minimal() + 
  facet_wrap(.~Moisture.treatment) +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"), panel.spacing = unit(1, "lines")) +
  ggpubr::stat_cor(method = "pearson")+ guides(colour = guide_legend(override.aes = list(size=3)))

m1

m3 <- ggplot(cordat9, aes(x = bulk.soil.cn_6mo, y = co2.native_6mo, size=0.5)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", linetype="dashed") +
  geom_point(alpha=0.7, aes(colour=Residue.type)) +
  scale_colour_manual(values=covercols)+
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) +
  xlab(expression(paste("Soil organic matter C:N")))+ 
  ylab(expression(paste("Priming (", mu,"g g"^"-1",")")))+ 
  theme_minimal() + 
  facet_wrap(.~Moisture.treatment) +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"), panel.spacing = unit(1, "lines")) +
  ggpubr::stat_cor(method = "pearson")+ guides(colour = guide_legend(override.aes = list(size=3)))

m3

figure <- ggpubr::ggarrange(m1, m2, m3, common.legend = T,legend="bottom",
                    labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)
png("figures/4_correlations/4_corplot_priming_cn.png", 
    width = 2000, height =3000, res=300, units="px")
figure
dev.off()










###### CN ratio change v. priming

m2 <- ggplot(cordat9, aes(x = pom.cn_6mo_change, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) + #xlim(c(15, 27)) +
  xlab(expression(paste("Change in POM C:N")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m2

m1 <- ggplot(cordat9, aes(x = maom.cn_6mo_change, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) +
  xlab(expression(paste("Change in MAOM C:N")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m1

m3 <- ggplot(cordat9, aes(x = bulk.soil.cn_6mo_change, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) + xlim (c(-1,0.7)) +
  xlab(expression(paste("Change in SOM C:N")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m3

figure <- ggarrange(m1, m2, m3, common.legend = T,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
png("figures/4_correlations/4_corplot_priming_cn_change.png", width = 3000, height =1000, res=300, units="px")
figure
dev.off()





###### nitrogen change v. priming


m2 <- ggplot(cordat9, aes(x = pom.n.g.kg_6mo_change, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) + #xlim(c(15, 27)) +
  xlab(expression(paste("Change in POM-N (g kg)"^"-1",")")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m2

m1 <- ggplot(cordat9, aes(x = maom.n.g.kg_6mo_change, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) +
  xlab(expression(paste("Change in MAOM-N (g kg)"^"-1",")")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m1

m3 <- ggplot(cordat9, aes(x = bulk.soil.N.g.kg_6mo_change, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) +
  xlab(expression(paste("Change in total soil N (g kg)"^"-1",")")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m3

figure <- ggarrange(m1, m2, m3, common.legend = T,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
png("figures/4_correlations/4_corplot_priming_n_change.png", width = 3000, height =1000, res=300, units="px")
figure
dev.off()



###### nitrogen v. priming


m2 <- ggplot(cordat9, aes(x = pom.n.g.kg_6mo, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) + #xlim(c(15, 27)) +
  xlab(expression(paste("POM-N (g kg)"^"-1",")")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m2

m1 <- ggplot(cordat9, aes(x = maom.n.g.kg_6mo, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) +
  xlab(expression(paste("MAOM-N (g kg)"^"-1",")")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m1

m3 <- ggplot(cordat9, aes(x = bulk.soil.N.g.kg_6mo, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) +
  xlab(expression(paste("Total soil N (g kg)"^"-1",")")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m3

figure <- ggarrange(m1, m2, m3, common.legend = T,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
png("figures/4_correlations/4_corplot_priming_n.png", width = 3000, height =1000, res=300, units="px")
figure
dev.off()






###### carbon v. priming


m2 <- ggplot(cordat9, aes(x = pom.c.g.kg_6mo, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) + #xlim(c(15, 27)) +
  xlab(expression(paste("POM-C (g kg)"^"-1",")")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m2

m1 <- ggplot(cordat9, aes(x = maom.c.g.kg_6mo, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) +
  xlab(expression(paste("MAOM-C (g kg)"^"-1",")")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m1

m3 <- ggplot(cordat9, aes(x = bulk.soil.C.g.kg_6mo, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) +
  xlab(expression(paste("Total SOC (g kg)"^"-1",")")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m3

figure <- ggarrange(m1, m2, m3, common.legend = T,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
png("figures/4_correlations/4_corplot_priming_C.png", width = 3000, height =1000, res=300, units="px")
figure
dev.off()




###### carbon change v. priming


m2 <- ggplot(cordat9, aes(x = pom.c.g.kg_6mo_change, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) + #xlim(c(15, 27)) +
  xlab(expression(paste("Change in POM-C (g kg)"^"-1",")")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m2

m1 <- ggplot(cordat9, aes(x = maom.c.g.kg_6mo_change, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) +
  xlab(expression(paste("Change in MAOM-C (g kg)"^"-1",")")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m1

m3 <- ggplot(cordat9, aes(x = bulk.soil.C.g.kg_6mo_change, y = co2.native_6mo, size=0.5, group=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + ylim(c(-150, 150)) + xlim(c(-1.5, 1)) +
  xlab(expression(paste("Change in total SOC (g kg)"^"-1",")")))+ 
  ylab(expression(paste("Priming (", mu,"g g)"^"-1",")")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m3

figure <- ggarrange(m1, m2, m3, common.legend = T,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
png("figures/4_correlations/4_corplot_priming_C_change.png", width = 3000, height =1000, res=300, units="px")
figure
dev.off()





###### maom v. CSE


m2 <- ggplot(cordat9, aes(x = maom.c.g.kg_6mo, y = CSE_6mo, group=Moisture.treatment, size=0.5)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top"),
         group = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + #ylim(c(-150, 150)) + #xlim(c(15, 27)) +
  xlab(expression(paste("MAOM-C (g kg)"^"-1",")")))+ 
  ylab(expression(paste("CSE")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m2

m1 <- ggplot(cordat9, aes(x = maom.n.g.kg_6mo, y = CSE_6mo, group=Moisture.treatment, size=0.5)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top"),
         group = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + #ylim(c(-150, 150)) +
  xlab(expression(paste("MAOM-N (g kg)"^"-1",")")))+ 
  ylab(expression(paste("CSE")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m1

m3 <- ggplot(cordat9, aes(x = maom.cn_6mo, y = CSE_6mo, group=Moisture.treatment, size=0.5)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top"),
         group = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + #ylim(c(-150, 150)) + xlim(c(-1.5, 1)) +
  xlab(expression(paste("MAOM C:N")))+ 
  ylab(expression(paste("CSE")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m3

figure <- ggarrange(m1, m2, m3, common.legend = T,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
png("figures/4_correlations/4_corplot_CSE_maom.png", width = 3500, height =1200, res=300, units="px")
figure
dev.off()



###### pom v. CSE


m2 <- ggplot(cordat9, aes(x = pom.c.g.kg_6mo, y = CSE_6mo, group=Moisture.treatment, size=0.5)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top"),
         group = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + #ylim(c(-150, 150)) + #xlim(c(15, 27)) +
  xlab(expression(paste("POM-C (g kg)"^"-1",")")))+ 
  ylab(expression(paste("CSE")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m2

m1 <- ggplot(cordat9, aes(x = pom.n.g.kg_6mo, y = CSE_6mo, group=Moisture.treatment, size=0.5)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top"),
         group = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + #ylim(c(-150, 150)) +
  xlab(expression(paste("POM-N (g kg)"^"-1",")")))+ 
  ylab(expression(paste("CSE")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson")

m1

m3 <- ggplot(cordat9, aes(x = pom.cn_6mo, y = CSE_6mo, group=Moisture.treatment, size=0.5)) +
  geom_point(alpha=0.7, aes(colour=Residue.type, shape=Moisture.treatment)) +
  scale_colour_manual(values=covercols)+
  scale_shape_manual(values=c(19,21))+
  geom_smooth(method='lm', se = FALSE, size=1, color="black", aes(linetype=Moisture.treatment)) +
  guides(colour = guide_legend(title.position="top"),
         shape = guide_legend(title.position="top"),
         group = guide_legend(title.position="top")) +
  labs(colour="Residue type", shape="Moisture treatment")+
  scale_size(guide = 'none') + #ylim(c(-150, 150)) + xlim(c(-1.5, 1)) +
  xlab(expression(paste("POM C:N")))+ 
  ylab(expression(paste("CSE")))+ 
  theme_minimal() + 
  #facet_wrap(.~Residue.type) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ggpubr::stat_cor(method = "pearson") +
  xlim(c(15,27))

m3

figure <- ggarrange(m1, m2, m3, common.legend = T,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
png("figures/4_correlations/4_corplot_CSE_pom.png", width = 3500, height =1200, res=300, units="px")
figure
dev.off()
