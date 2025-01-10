
library(dplyr)
library(ggplot2)
library(nlme)
library(emmeans)
library(multcomp)

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
dat <- read.csv("raw-data/master-dataset-13c.csv")
# set groupings
dat$Residue.type <- as.factor(dat$Residue.type)
dat$Residue.type <- factor(dat$Residue.type, levels(dat$Residue.type)[c(2,5,1,4,3)])
levels(dat$Residue.type) <- c("No residue", "Wheat", "Clover", "Wheat-clover mix", "Five-species mix")
dat$Moisture.treatment <- factor(dat$Moisture.treatment)
dat$Harvest.day <- as.factor(dat$Harvest.day)
dat$Harvest.day2 <- dat$Harvest.day
levels(dat$Harvest.day2) <- c("24 hr", "6 mo")

# remove baseline values
bdat <- dat[which(dat$Sample.number %in% c("baseline1", "baseline2", "baseline3")),]
dat <- dat[-which(dat$Sample.number %in% c("baseline1", "baseline2", "baseline3")),]





#######################################
# export treatment averages of control (no 13c-labeled glucose, just water) soils
cdat <- dat[which(dat$Glucose.addition.treatment=="G-00"),]
cdatav <- cdat %>%
  group_by(Harvest.day2, Residue.type, Moisture.treatment) %>% 
  summarise(across(
    .cols = is.numeric, 
    .fns = list(Mean = mean, SD = sd), na.rm = TRUE, 
    .names = "{col}_{fn}"
  ))
write.csv(cdatav, "summary-stats/means of no glucose soils at t1 and t2.csv")




##### anova output dataframe ##### 
p<-24+4
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












####### bulk.soil.C.g.kg
n <- 1
mod.anova2[1,n] <- "bulk.soil.C.g.kg_24 hr"
hist(dat$bulk.soil.C.g.kg[which(dat$Harvest.day2=="24 hr")])
fit <- lme(log(bulk.soil.C.g.kg) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$bulk.soil.C.g.kg)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(bulk.soil.C.g.kg),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/bulk.soil.C.g.kg_24 hr_2-way.csv")
testlet_bulk.soil.C.g.kg_24hr <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/bulk.soil.C.g.kg_24 hr_2-way_mois.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/bulk.soil.C.g.kg_24 hr_2-way_res.csv")


n <- 2
mod.anova2[1,n] <- "bulk.soil.C.g.kg_6 mo"
hist(dat$bulk.soil.C.g.kg[which(dat$Harvest.day2=="6 mo")])
fit <- lme(log(bulk.soil.C.g.kg) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$bulk.soil.C.g.kg)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(bulk.soil.C.g.kg),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/bulk.soil.C.g.kg_6 mo_2-way.csv")
testlet_bulk.soil.C.g.kg_6mo <- testlet

### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/bulk.soil.C.g.kg_6mo_2-way_mois.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/bulk.soil.C.g.kg_6mo_2-way_res.csv")



testlet_bulk.soil.C.g.kg <- rbind(testlet_bulk.soil.C.g.kg_24hr, testlet_bulk.soil.C.g.kg_6mo)
testlet_bulk.soil.C.g.kg$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_bulk.soil.C.g.kg$minys <- testlet_bulk.soil.C.g.kg$response-testlet_bulk.soil.C.g.kg$SE
testlet_bulk.soil.C.g.kg$minys[which(testlet_bulk.soil.C.g.kg$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_bulk.soil.C.g.kg, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("bulk.soil.C.g.kg"^""))) +   
  facet_grid(.~Harvest.day2) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_bulk.soil.C.g.kg$response+testlet_bulk.soil.C.g.kg$SE)+0.2*max(testlet_bulk.soil.C.g.kg$response+testlet_bulk.soil.C.g.kg$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) 
p

q <- p +   ylim(c(0,21.77)) +
  annotate("segment", x=0.4, y=11.77, xend=0.5, yend=11.77,
                    col="brown", arrow=arrow(length=unit(0.1, "cm"))) +   
  annotate("segment", x=0.4, y=21.77, xend=0.5, yend=21.77,
           col="royalblue", arrow=arrow(length=unit(0.1, "cm"))) 
q
ggpubr::ggexport(q, height=1300, width=2600, filename = "figures/*pools/1_bulk.soil.C.g.kg_2wayb.png", res = 400)

# The two lines we want on the plot
fig_bulk.soil.C.g.kg3 <- p +
  geom_text(aes(label=.group, y=1+testlet_bulk.soil.C.g.kg$response+testlet_bulk.soil.C.g.kg$SE), 
            position = position_dodge(0.9), size=2.5)
# export figure
ggpubr::ggexport(fig_bulk.soil.C.g.kg3, height=1300, width=2600, filename = "figures/*pools/1_bulk.soil.C.g.kg_2way.png", res = 400)


####### bulk.soil.N.g.kg
n <- 3
mod.anova2[1,n] <- "bulk.soil.N.g.kg_24 hr"
hist(dat$bulk.soil.N.g.kg[which(dat$Harvest.day2=="24 hr")])
fit <- lme(log(bulk.soil.N.g.kg) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$bulk.soil.N.g.kg)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(bulk.soil.N.g.kg),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uva/*pools/bulk.soil.N.g.kg_24 hr_2-way.csv")
testlet_bulk.soil.N.g.kg_24hr <- testlet

n <- 4
mod.anova2[1,n] <- "bulk.soil.N.g.kg_6 mo"
hist(dat$bulk.soil.N.g.kg[which(dat$Harvest.day2=="6 mo")])
fit <- lme((bulk.soil.N.g.kg) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$bulk.soil.N.g.kg)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(bulk.soil.N.g.kg),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/bulk.soil.N.g.kg_6 mo_2-way.csv")
testlet_bulk.soil.N.g.kg_6mo <- testlet



testlet_bulk.soil.N.g.kg <- rbind(testlet_bulk.soil.N.g.kg_24hr, testlet_bulk.soil.N.g.kg_6mo)
testlet_bulk.soil.N.g.kg$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_bulk.soil.N.g.kg$minys <- testlet_bulk.soil.N.g.kg$response-testlet_bulk.soil.N.g.kg$SE
testlet_bulk.soil.N.g.kg$minys[which(testlet_bulk.soil.N.g.kg$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_bulk.soil.N.g.kg, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("bulk.soil.N.g.kg"^""))) +   
  facet_grid(.~Harvest.day2) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_bulk.soil.N.g.kg$response+testlet_bulk.soil.N.g.kg$SE)+0.2*max(testlet_bulk.soil.N.g.kg$response+testlet_bulk.soil.N.g.kg$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=0.2+testlet_bulk.soil.N.g.kg$response+testlet_bulk.soil.N.g.kg$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_bulk.soil.N.g.kg3 <- p
# export figure
ggpubr::ggexport(fig_bulk.soil.N.g.kg3, height=1300, width=2600, filename = "figures/*pools/1_bulk.soil.N.g.kg_2way.png", res = 400)



####### bulk.soil.cn
n <- 5
mod.anova2[1,n] <- "bulk.soil.cn_24 hr"
hist(dat$bulk.soil.cn[which(dat$Harvest.day2=="24 hr")])
fit <- lme((bulk.soil.cn) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$bulk.soil.cn)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(bulk.soil.cn),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/*uv/lsmeans/*pools/bulk.soil.cn_24 hr_2-way.csv")
testlet_bulk.soil.cn_24hr <- testlet

n <- 6
mod.anova2[1,n] <- "bulk.soil.cn_6 mo"
hist(dat$bulk.soil.cn[which(dat$Harvest.day2=="6 mo")])
fit <- lme(log(bulk.soil.cn) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$bulk.soil.cn)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(bulk.soil.cn),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/bulk.soil.cn_6 mo_2-way.csv")
testlet_bulk.soil.cn_6mo <- testlet



testlet_bulk.soil.cn <- rbind(testlet_bulk.soil.cn_24hr, testlet_bulk.soil.cn_6mo)
testlet_bulk.soil.cn$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_bulk.soil.cn$minys <- testlet_bulk.soil.cn$response-testlet_bulk.soil.cn$SE
testlet_bulk.soil.cn$minys[which(testlet_bulk.soil.cn$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_bulk.soil.cn, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("bulk.soil.cn"^""))) +   
  facet_grid(.~Harvest.day2) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_bulk.soil.cn$response+testlet_bulk.soil.cn$SE)+0.2*max(testlet_bulk.soil.cn$response+testlet_bulk.soil.cn$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=0.9+testlet_bulk.soil.cn$response+testlet_bulk.soil.cn$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_bulk.soil.cn3 <- p
# export figure
ggpubr::ggexport(fig_bulk.soil.cn3, height=1300, width=2600, filename = "figures/*pools/1_bulk.soil.cn_2way.png", res = 400)














####### maom.c.g.kg
n <- 1+6
mod.anova2[1,n] <- "maom.c.g.kg_24 hr"
hist(dat$maom.c.g.kg[which(dat$Harvest.day2=="24 hr")])
fit <- lme((maom.c.g.kg) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$maom.c.g.kg)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(maom.c.g.kg),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/maom.c.g.kg_24 hr_2-way.csv")
testlet_maom.c.g.kg_24hr <- testlet

n <- 2+6
mod.anova2[1,n] <- "maom.c.g.kg_6 mo"
hist(dat$maom.c.g.kg[which(dat$Harvest.day2=="6 mo")])
fit <- lme((maom.c.g.kg) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$maom.c.g.kg)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(maom.c.g.kg),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/maom.c.g.kg_6 mo_2-way.csv")
testlet_maom.c.g.kg_6mo <- testlet



testlet_maom.c.g.kg <- rbind(testlet_maom.c.g.kg_24hr, testlet_maom.c.g.kg_6mo)
testlet_maom.c.g.kg$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_maom.c.g.kg$minys <- testlet_maom.c.g.kg$response-testlet_maom.c.g.kg$SE
testlet_maom.c.g.kg$minys[which(testlet_maom.c.g.kg$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_maom.c.g.kg, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("maom.c.g.kg"^""))) +   
  facet_grid(.~Harvest.day2) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_maom.c.g.kg$response+testlet_maom.c.g.kg$SE)+0.2*max(testlet_maom.c.g.kg$response+testlet_maom.c.g.kg$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=0.8+testlet_maom.c.g.kg$response+testlet_maom.c.g.kg$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_maom.c.g.kg3 <- p
# export figure
ggpubr::ggexport(fig_maom.c.g.kg3, height=1300, width=2600, filename = "figures/*pools/1_maom.c.g.kg_2way.png", res = 400)


####### maom.n.g.kg
n <- 3+6
mod.anova2[1,n] <- "maom.n.g.kg_24 hr"
hist(dat$maom.n.g.kg[which(dat$Harvest.day2=="24 hr")])
fit <- lme(log(maom.n.g.kg) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$maom.n.g.kg)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(maom.n.g.kg),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/maom.n.g.kg_24 hr_2-way.csv")
testlet_maom.n.g.kg_24hr <- testlet

n <- 4+6
mod.anova2[1,n] <- "maom.n.g.kg_6 mo"
hist(dat$maom.n.g.kg[which(dat$Harvest.day2=="6 mo")])
fit <- lme((maom.n.g.kg) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$maom.n.g.kg)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(maom.n.g.kg),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/maom.n.g.kg_6 mo_2-way.csv")
testlet_maom.n.g.kg_6mo <- testlet



testlet_maom.n.g.kg <- rbind(testlet_maom.n.g.kg_24hr, testlet_maom.n.g.kg_6mo)
testlet_maom.n.g.kg$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_maom.n.g.kg$minys <- testlet_maom.n.g.kg$response-testlet_maom.n.g.kg$SE
testlet_maom.n.g.kg$minys[which(testlet_maom.n.g.kg$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_maom.n.g.kg, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("maom.n.g.kg"^""))) +   
  facet_grid(.~Harvest.day2) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_maom.n.g.kg$response+testlet_maom.n.g.kg$SE)+0.2*max(testlet_maom.n.g.kg$response+testlet_maom.n.g.kg$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=0.1+testlet_maom.n.g.kg$response+testlet_maom.n.g.kg$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_maom.n.g.kg3 <- p
# export figure
ggpubr::ggexport(fig_maom.n.g.kg3, height=1300, width=2600, filename = "figures/*pools/1_maom.n.g.kg_2way.png", res = 400)



####### maom.cn
n <- 5+6
mod.anova2[1,n] <- "maom.cn_24 hr"
hist(dat$maom.cn[which(dat$Harvest.day2=="24 hr")])
fit <- lme((maom.cn) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$maom.cn)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(maom.cn),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/*uv/lsmeans/*pools/maom.cn_24 hr_2-way.csv")
testlet_maom.cn_24hr <- testlet

n <- 6+6
mod.anova2[1,n] <- "maom.cn_6 mo"
hist(dat$maom.cn[which(dat$Harvest.day2=="6 mo")])
fit <- lme(log(maom.cn) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$maom.cn)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(maom.cn),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEF", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/maom.cn_6 mo_2-way.csv")
testlet_maom.cn_6mo <- testlet



testlet_maom.cn <- rbind(testlet_maom.cn_24hr, testlet_maom.cn_6mo)
testlet_maom.cn$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_maom.cn$minys <- testlet_maom.cn$response-testlet_maom.cn$SE
testlet_maom.cn$minys[which(testlet_maom.cn$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_maom.cn, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("maom.cn"^""))) +   
  facet_grid(.~Harvest.day2) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_maom.cn$response+testlet_maom.cn$SE)+0.2*max(testlet_maom.cn$response+testlet_maom.cn$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=0.9+testlet_maom.cn$response+testlet_maom.cn$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_maom.cn3 <- p
# export figure
ggpubr::ggexport(fig_maom.cn3, height=1300, width=2600, filename = "figures/*pools/1_maom.cn_2way.png", res = 400)
















####### pom.c.g.kg
n <- 1+12
mod.anova2[1,n] <- "pom.c.g.kg_24 hr"
hist(dat$pom.c.g.kg[which(dat$Harvest.day2=="24 hr")])
fit <- lme(log(pom.c.g.kg) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$pom.c.g.kg)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(pom.c.g.kg),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/pom.c.g.kg_24 hr_2-way.csv")
testlet_pom.c.g.kg_24hr <- testlet

n <- 2+12
mod.anova2[1,n] <- "pom.c.g.kg_6 mo"
hist(dat$pom.c.g.kg[which(dat$Harvest.day2=="6 mo")])
fit <- lme(log(pom.c.g.kg) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$pom.c.g.kg)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(pom.c.g.kg),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/pom.c.g.kg_6 mo_2-way.csv")
testlet_pom.c.g.kg_6mo <- testlet



testlet_pom.c.g.kg <- rbind(testlet_pom.c.g.kg_24hr, testlet_pom.c.g.kg_6mo)
testlet_pom.c.g.kg$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_pom.c.g.kg$minys <- testlet_pom.c.g.kg$response-testlet_pom.c.g.kg$SE
testlet_pom.c.g.kg$minys[which(testlet_pom.c.g.kg$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_pom.c.g.kg, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("pom.c.g.kg"^""))) +   
  facet_grid(.~Harvest.day2) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_pom.c.g.kg$response+testlet_pom.c.g.kg$SE)+0.2*max(testlet_pom.c.g.kg$response+testlet_pom.c.g.kg$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=0.4+testlet_pom.c.g.kg$response+testlet_pom.c.g.kg$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_pom.c.g.kg3 <- p
# export figure
ggpubr::ggexport(fig_pom.c.g.kg3, height=1300, width=2600, filename = "figures/*pools/1_pom.c.g.kg_2way.png", res = 400)


####### pom.n.g.kg
n <- 3+12
mod.anova2[1,n] <- "pom.n.g.kg_24 hr"
hist(dat$pom.n.g.kg[which(dat$Harvest.day2=="24 hr")])
fit <- lme(log(pom.n.g.kg) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$pom.n.g.kg)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(pom.n.g.kg),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/pom.n.g.kg_24 hr_2-way.csv")
testlet_pom.n.g.kg_24hr <- testlet

n <- 4+12
mod.anova2[1,n] <- "pom.n.g.kg_6 mo"
hist(dat$pom.n.g.kg[which(dat$Harvest.day2=="6 mo")])
fit <- lme(log(pom.n.g.kg) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$pom.n.g.kg)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(pom.n.g.kg),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/pom.n.g.kg_6 mo_2-way.csv")
testlet_pom.n.g.kg_6mo <- testlet



testlet_pom.n.g.kg <- rbind(testlet_pom.n.g.kg_24hr, testlet_pom.n.g.kg_6mo)
testlet_pom.n.g.kg$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_pom.n.g.kg$minys <- testlet_pom.n.g.kg$response-testlet_pom.n.g.kg$SE
testlet_pom.n.g.kg$minys[which(testlet_pom.n.g.kg$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_pom.n.g.kg, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("pom.n.g.kg"^""))) +   
  facet_grid(.~Harvest.day2) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_pom.n.g.kg$response+testlet_pom.n.g.kg$SE)+0.2*max(testlet_pom.n.g.kg$response+testlet_pom.n.g.kg$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=0.03+testlet_pom.n.g.kg$response+testlet_pom.n.g.kg$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_pom.n.g.kg3 <- p
# export figure
ggpubr::ggexport(fig_pom.n.g.kg3, height=1300, width=2600, filename = "figures/*pools/1_pom.n.g.kg_2way.png", res = 400)



####### pom.cn
n <- 5+12
mod.anova2[1,n] <- "pom.cn_24 hr"
hist(dat$pom.cn[which(dat$Harvest.day2=="24 hr")])
fit <- lme(log(pom.cn) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$pom.cn)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(pom.cn),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/*uv/lsmeans/*pools/pom.cn_24 hr_2-way.csv")
testlet_pom.cn_24hr <- testlet

n <- 6+12
mod.anova2[1,n] <- "pom.cn_6 mo"
hist(dat$pom.cn[which(dat$Harvest.day2=="6 mo")])
fit <- lme(log(pom.cn) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$pom.cn)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(pom.cn),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDEF", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/pom.cn_6 mo_2-way.csv")
testlet_pom.cn_6mo <- testlet



testlet_pom.cn <- rbind(testlet_pom.cn_24hr, testlet_pom.cn_6mo)
testlet_pom.cn$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_pom.cn$minys <- testlet_pom.cn$response-testlet_pom.cn$SE
testlet_pom.cn$minys[which(testlet_pom.cn$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_pom.cn, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("pom.cn"^""))) +   
  facet_grid(.~Harvest.day2) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_pom.cn$response+testlet_pom.cn$SE)+0.2*max(testlet_pom.cn$response+testlet_pom.cn$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=2+testlet_pom.cn$response+testlet_pom.cn$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_pom.cn3 <- p
# export figure
ggpubr::ggexport(fig_pom.cn3, height=1300, width=2600, filename = "figures/*pools/1_pom.cn_2way.png", res = 400)
















####### mbc
n <- 1+18
mod.anova2[1,n] <- "mbc_24 hr"
hist(dat$mbc[which(dat$Harvest.day2=="24 hr")])
fit <- lme((mbc) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$mbc)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(mbc),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/mbc_24 hr_2-way.csv")
testlet_mbc_24hr <- testlet

n <- 2+18
mod.anova2[1,n] <- "mbc_6 mo"
hist(dat$mbc[which(dat$Harvest.day2=="6 mo")])
fit <- lme((mbc) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$mbc)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(mbc),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/mbc_6 mo_2-way.csv")
testlet_mbc_6mo <- testlet



testlet_mbc <- rbind(testlet_mbc_24hr, testlet_mbc_6mo)
testlet_mbc$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_mbc$minys <- testlet_mbc$response-testlet_mbc$SE
testlet_mbc$minys[which(testlet_mbc$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_mbc, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("mbc"^""))) +   
  facet_grid(.~Harvest.day2) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_mbc$response+testlet_mbc$SE)+0.2*max(testlet_mbc$response+testlet_mbc$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=50+testlet_mbc$response+testlet_mbc$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_mbc3 <- p
# export figure
ggpubr::ggexport(fig_mbc3, height=1300, width=2600, filename = "figures/*pools/1_mbc_2way.png", res = 400)


####### doc
n <- 3+18
mod.anova2[1,n] <- "doc_24 hr"
hist(dat$doc[which(dat$Harvest.day2=="24 hr")])
fit <- lme((doc) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$doc)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(doc),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/doc_24 hr_2-way.csv")
testlet_doc_24hr <- testlet

n <- 4+18
mod.anova2[1,n] <- "doc_6 mo"
hist(dat$doc[which(dat$Harvest.day2=="6 mo")])
fit <- lme(log(doc) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$doc)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(doc),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/doc_6 mo_2-way.csv")
testlet_doc_6mo <- testlet



testlet_doc <- rbind(testlet_doc_24hr, testlet_doc_6mo)
testlet_doc$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_doc$minys <- testlet_doc$response-testlet_doc$SE
testlet_doc$minys[which(testlet_doc$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_doc, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("doc"^""))) +   
  facet_grid(.~Harvest.day2) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_doc$response+testlet_doc$SE)+0.2*max(testlet_doc$response+testlet_doc$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=5+testlet_doc$response+testlet_doc$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_doc3 <- p
# export figure
ggpubr::ggexport(fig_doc3, height=1300, width=2600, filename = "figures/*pools/1_doc_2way.png", res = 400)






####### co2
n <- 5+18
mod.anova2[1,n] <- "co2_24 hr"
hist(dat$co2[which(dat$Harvest.day2=="24 hr")])
fit <- lme(log(co2) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$co2)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(co2),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/co2_24 hr_2-way.csv")
testlet_co2_24hr <- testlet

n <- 6+18
mod.anova2[1,n] <- "co2_6 mo"
hist(dat$co2[which(dat$Harvest.day2=="6 mo")])
fit <- lme(log(co2) ~ Moisture.treatment*Residue.type, 
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$co2)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(co2),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/*uv/lsmeans/*pools/co2_6 mo_2-way.csv")
testlet_co2_6mo <- testlet



testlet_co2 <- rbind(testlet_co2_24hr, testlet_co2_6mo)
testlet_co2$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_co2$minys <- testlet_co2$response-testlet_co2$SE
testlet_co2$minys[which(testlet_co2$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_co2, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("co2"^""))) +   
  facet_grid(Harvest.day2~., scales="free") + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  #ylim(c(0,max(testlet_co2$response+testlet_co2$SE)+0.2*max(testlet_co2$response+testlet_co2$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=c(rep(5,10), rep(200,10))+testlet_co2$response+testlet_co2$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_co23 <- p
# export figure
ggpubr::ggexport(fig_co23, height=2000, width=2000, filename = "figures/*pools/1_co2_2way.png", res = 400)







 



# export anova results
write.csv(mod.anova2, "model-output/anova_2way_native.csv")



all_lets <- bind_rows(testlet_bulk.soil.C.g.kg, 
          testlet_bulk.soil.N.g.kg, 
          testlet_bulk.soil.cn, 
          testlet_maom.c.g.kg, 
          testlet_maom.n.g.kg, 
          testlet_maom.cn, 
          testlet_pom.c.g.kg, 
          testlet_pom.n.g.kg, 
          testlet_pom.cn, 
          testlet_mbc, 
          testlet_doc, 
          testlet_co2)
# export tukey results
all_lets$condition <- rep(c("soc", "n", "cn", "maom-c", "maom-n", "maom-cn", "pom-c", "pom-n", "pom-cn","mbc", "doc", "co2"), each=20)
write.csv(all_lets, "model-output/*uv/lsmeans/*pools/*combined_2way_native.csv")
