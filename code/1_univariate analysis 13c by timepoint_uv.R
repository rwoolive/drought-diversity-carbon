# to allow unequal variances:
# weights=varIdent(form=~1|Residue.type) 

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




##### anova output dataframe ##### 
p<-50
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


# calculate recovery of glucose-c

dat$recovery <- dat$X13co2 + dat$X13maoc + dat$X13poc 
recovery_table <- dat %>% 
  dplyr::summarise(pct_mean_recovery = 100*mean(recovery/50, na.rm=T),
            pct_mean_recovery_sd = 100*sd(recovery/50, na.rm=T),
            pct_mean_recovery_se = 100*sd(recovery/50, na.rm=T)/n())
write.csv(recovery_table, "model-output/recovery of glucose-c in CO2, MAOC, POC across all.csv")




####### CUE -- conventional calculation
n <- 1
mod.anova2[1,n] <- "CUE_24 hr"
hist(dat$CUE[which(dat$Harvest.day2=="24 hr")])
fit <- lme(log(CUE) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$CUE)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(CUE),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/CUE_24 hr_2-way.csv")
testlet_CUE_24hr <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/CUE_24 hr_2-way_mois.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/CUE_24 hr_2-way_res.csv")

n <- 2
mod.anova2[1,n] <- "CUE_6 mo"
hist(dat$CUE[which(dat$Harvest.day2=="6 mo")])
fit <- lme(log(CUE) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$CUE)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(CUE),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/CUE_6 mo_2-way.csv")
testlet_CUE_6mo <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/CUE_6 mo_2-way_mois.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/CUE_6 mo_2-way_res.csv")



testlet_CUE <- rbind(testlet_CUE_24hr, testlet_CUE_6mo)
testlet_CUE$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_CUE$minys <- testlet_CUE$response-testlet_CUE$SE
#testlet_CUE$minys[which(testlet_CUE$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_CUE, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("Microbial CUE"^""))) +   
  facet_wrap(.~Harvest.day2, scales="free_y") + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  #ylim(c(0,max(testlet_CUE$response+testlet_CUE$SE)+0.2*max(testlet_CUE$response+testlet_CUE$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=c(rep(NA,10), testlet_CUE$.group[11:20]), y=0.003+c(rep(0,10), testlet_CUE$response[11:20])+c(rep(0,10), testlet_CUE$SE[11:20])), 
            position = position_dodge(0.9), size=2)
p

# The two lines we want on the plot
fig_CUE3 <- p
# export figure
ggpubr::ggexport(fig_CUE3, height=1300, width=2600, filename = "figures/_uv/1_cue_2way.png", res = 400)


# plot
p <- ggplot(data=testlet_CUE, aes(x=Residue.type, y=response, fill=Moisture.treatment)) +
  theme_minimal() + labs(fill="Moisture treatment") +
  ylab(expression(paste("Microbial CUE"^""))) +   
  facet_wrap(.~Harvest.day2) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) + # , scales="free_y"
  geom_bar(stat = "identity", position = position_dodge(0.8), color="black", width=0.8) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.8)) +
  scale_fill_manual(values = c("blue", "red")) +
  #ylim(c(0,max(testlet_CUE$response+testlet_CUE$SE)+0.2*max(testlet_CUE$response+testlet_CUE$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.1, 0.5), "cm"), axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0)) +
  geom_text(aes(label=c(rep(NA,10), testlet_CUE$.group[11:20]), y=0.01+c(rep(0,10), testlet_CUE$response[11:20])+c(rep(0,10), testlet_CUE$SE[11:20])), 
            position = position_dodge(0.8), size=2)
p

# The two lines we want on the plot
fig_CUEb <- p
# export figure
ggpubr::ggexport(fig_CUEb, height=1700, width=2700, filename = "figures/_uv/1_cue_2way_b.png", res = 400)





####### CUE_r --  calculation including microbial residue-C in the denominator
n <- 27
mod.anova2[1,n] <- "CUE_r_24 hr"
dat$X13residue <- dat$X13soc - dat$X13mbc
dat$CUE_r <- dat$X13mbc/(dat$X13mbc + dat$X13co2 + dat$X13residue)
hist(dat$CUE_r[which(dat$Harvest.day2=="24 hr")])
fit <- lme(log(CUE_r) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$CUE_r)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(CUE_r),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/CUE_r_24 hr_2-way.csv")
testlet_CUE_r_24hr <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/CUE_r_24 hr_2-way_mois.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/CUE_r_24 hr_2-way_res.csv")

n <- 28
mod.anova2[1,n] <- "CUE_r_6 mo"
hist(dat$CUE_r[which(dat$Harvest.day2=="6 mo")])
fit <- lme(log(CUE_r) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$CUE_r)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(CUE_r),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/CUE_r_6 mo_2-way.csv")
testlet_CUE_r_6mo <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/CUE_r_6 mo_2-way_mois.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/CUE_r_6 mo_2-way_res.csv")



testlet_CUE_r <- rbind(testlet_CUE_r_24hr, testlet_CUE_r_6mo)
testlet_CUE_r$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_CUE_r$minys <- testlet_CUE_r$response-testlet_CUE_r$SE
#testlet_CUE_r$minys[which(testlet_CUE_r$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_CUE_r, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("Microbial CUE_r"^""))) +   
  facet_wrap(.~Harvest.day2, scales="free_y") + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  #ylim(c(0,max(testlet_CUE_r$response+testlet_CUE_r$SE)+0.2*max(testlet_CUE_r$response+testlet_CUE_r$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=c(rep(NA,10), testlet_CUE_r$.group[11:20]), y=0.003+c(rep(0,10), testlet_CUE_r$response[11:20])+c(rep(0,10), testlet_CUE_r$SE[11:20])), 
            position = position_dodge(0.9), size=2)
p

# The two lines we want on the plot
fig_CUE_r3 <- p
# export figure
ggpubr::ggexport(fig_CUE_r3, height=1300, width=2600, filename = "figures/_uv/1_CUE_r_2way.png", res = 400)


# plot
p <- ggplot(data=testlet_CUE_r, aes(x=Residue.type, y=response, fill=Moisture.treatment)) +
  theme_minimal() + labs(fill="Moisture treatment") +
  ylab(expression(paste("Microbial CUE"[r]))) +   
  facet_wrap(.~Harvest.day2) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) + # , scales="free_y"
  geom_bar(stat = "identity", position = position_dodge(0.8), color="black", width=0.8) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.8)) +
  scale_fill_manual(values = c("blue", "red")) +
  #ylim(c(0,max(testlet_CUE_r$response+testlet_CUE_r$SE)+0.2*max(testlet_CUE_r$response+testlet_CUE_r$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.1, 0.5), "cm"), axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0)) +
  geom_text(aes(label=c(rep(NA,10), testlet_CUE_r$.group[11:20]), y=0.01+c(rep(0,10), testlet_CUE_r$response[11:20])+c(rep(0,10), testlet_CUE_r$SE[11:20])), 
            position = position_dodge(0.8), size=2)
p

# The two lines we want on the plot
fig_CUE_rb <- p
# export figure
ggpubr::ggexport(fig_CUE_rb, height=1700, width=2700, filename = "figures/_uv/1_CUE_r_2way_b.png", res = 400)




####### CSE -- calculation including microbial residue-C in the numerator and denominator
n <- 29
mod.anova2[1,n] <- "CSE_24 hr"
dat$X13residue <- dat$X13soc - dat$X13mbc
dat$CSE <- (dat$X13mbc+ dat$X13residue)/(dat$X13mbc + dat$X13co2 + dat$X13residue)
hist(dat$CSE[which(dat$Harvest.day2=="24 hr")])
fit <- lme((CSE) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$CSE)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(CSE),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/CSE_24 hr_2-way.csv")
testlet_CSE_24hr <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/CSE_24 hr_2-way_mois.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/CSE_24 hr_2-way_res.csv")
testlet_CSE_24hr_res <- testlet

# custom contrast:
test1
A = c(1, 0, 0, 0, 0)
B = c(0, 0.25, 0.25, 0.25, 0.25)
con <- contrast(test1, method = list("A - B" = A - B)) 
write.csv(con, "model-output/_uv/lsmeans/CSE_24 hr_2-way_res_custom contrast.csv")



n <- 30
mod.anova2[1,n] <- "CSE_6 mo"
hist(dat$CSE[which(dat$Harvest.day2=="6 mo")])
dat$CSE[which(dat$CSE<0)] <- NA
fit <- lme((CSE) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$CSE)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(CSE),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/CSE_6 mo_2-way.csv")
testlet_CSE_6mo <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/CSE_6 mo_2-way_mois.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/CSE_6 mo_2-way_res.csv")
testlet_CSE_6mo_res <- testlet


testlet_CSE <- rbind(testlet_CSE_24hr, testlet_CSE_6mo)
testlet_CSE$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_CSE$minys <- testlet_CSE$response-testlet_CSE$SE
#testlet_CSE$minys[which(testlet_CSE$minys<0)] <- 0

testlet_CSE_res <- rbind(testlet_CSE_24hr_res, testlet_CSE_6mo_res)
testlet_CSE_res$Harvest.day2 <- as.factor(c(rep("24 hr", 5), rep("6 mo", 5)))
testlet_CSE_res$minys <- testlet_CSE_res$response-testlet_CSE_res$SE
#testlet_CSE$minys[which(testlet_CSE$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_CSE, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("Microbial CSE"^""))) +   
  facet_wrap(.~Harvest.day2, scales="free_y") + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  #ylim(c(0,max(testlet_CSE$response+testlet_CSE$SE)+0.2*max(testlet_CSE$response+testlet_CSE$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) #+
  # geom_text(inherit.aes = F, data = testlet_CSE_res, 
  #           aes(label = c(testlet_CSE_res$.group[1:5], rep(NA,5)), 
  #               y = 0.003+c(testlet_CSE_res$response[1:5], rep(0,5)) + c(testlet_CSE_res$SE[1:5], rep(0,5)),
  #               x = rep(seq(from = 1, to = 2, length.out = 5), 2)), 
  #           position = position_dodge(0.9), size=2)
p

# The two lines we want on the plot
fig_CSE3 <- p
# export figure
ggpubr::ggexport(fig_CSE3, height=1300, width=2600, filename = "figures/_uv/1_CSE_2way.png", res = 400)


# plot
p <- ggplot(data=testlet_CSE, aes(x=Residue.type, y=response, fill=Moisture.treatment)) +
  theme_minimal() + labs(fill="Moisture treatment") +
  ylab(expression(paste("Microbial CSE"))) +   
  facet_wrap(.~Harvest.day2) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) + # , scales="free_y"
  geom_bar(stat = "identity", position = position_dodge(0.8), color="black", width=0.8) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.8)) +
  scale_fill_manual(values = c("blue", "red")) +
  #ylim(c(0,max(testlet_CSE$response+testlet_CSE$SE)+0.2*max(testlet_CSE$response+testlet_CSE$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.1, 0.5), "cm"), axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0)) +
  geom_text(inherit.aes = F, data = testlet_CSE_res, 
            aes(label = c(testlet_CSE_res$.group[1:5], rep(NA,5)), 
                y = c(rep(0.82, 5), rep(0, 5)), # 0.06+c(testlet_CSE_res$response[1:5], rep(0,5)) + c(testlet_CSE_res$SE[1:5], rep(0,5)),
                x = rep(seq(from = 1, to = 5, length.out = 5), 2)), 
            position = position_dodge(0.9), size=2)
p

# The two lines we want on the plot
fig_CSEb <- p
# export figure
ggpubr::ggexport(fig_CSEb, height=1700, width=2700, filename = "figures/_uv/1_CSE_2way_b.png", res = 400)













####### co2.native-- priming
n <- 15
mod.anova2[1,n] <- "co2.native_24 hr"
hist(dat$co2.native[which(dat$Harvest.day2=="24 hr")])
fit <- lme((co2.native) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),#weights=varIdent(form=~1|Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$co2.native)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(co2.native),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/co2.native_24 hr_2-way.csv")
testlet_co2.native_24hr <- testlet
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/co2.native_24 hr_2-way_res.csv")
testlet_co2.native_24hr_res <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/co2.native_24 hr_2-way_mois.csv")
testlet_co2.native_24hr_mois <- testlet


n <- 16
mod.anova2[1,n] <- "co2.native_6 mo"
hist(dat$co2.native[which(dat$Harvest.day2=="6 mo")])
dat[which(dat$Harvest.day2=="6 mo" & dat$co2.native < -100),]
`%between%` = function(x,range) x>range[1] & x<range[2]
dat$co2.native[which(dat$co2.native %between% c(-250,-100))] <- NA
fit <- lme((co2.native) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),#weights=varIdent(form=~1|Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$co2.native)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(co2.native),funs(mean=mean(., na.rm=T),Std=sd(.)))
mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
# export model stats
mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCD", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/co2.native_6 mo_2-way.csv")
testlet_co2.native_6mo <- testlet
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/co2.native_6 mo_2-way_res.csv")
testlet_co2.native_6mo_res <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/co2.native_6 mo_2-way_mois.csv")
testlet_co2.native_6mo_mois <- testlet


testlet_co2.native <- rbind(testlet_co2.native_24hr, testlet_co2.native_6mo)
testlet_co2.native$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_co2.native$minys <- testlet_co2.native$response-testlet_co2.native$SE
# testlet_co2.native$minys[which(testlet_co2.native$minys<0)] <- 0

testlet_co2.native$.group[1:10] <- rep(NA, 10)


# plot
testlet_co2.native$y_position <- c(rep(0, 10), rep(10, 10)) + testlet_co2.native$response + testlet_co2.native$SE

p <- ggplot(data=testlet_co2.native, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(label=expression(paste("Priming (", mu,"g g"^-1,")"))) +   
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  geom_text(aes(label=.group, y=y_position), position = position_dodge(0.9), size=2.5)+ 
  #scale_y_continuous(breaks=seq(-100,200, 50)) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank(),
        legend.text = element_text(size = 8),  # Change text size
        legend.title = element_text(size = 9)) +
  ggh4x::facet_nested(~Harvest.day2, scales = "free_y", independent="y")  #+
#geom_text(inherit.aes=FALSE, data=testlet_X13maoc,  aes(label=rep(NA,20), x=rep(1:5, 4), y=c(rep(-10,10),rep(10,10)))

#facet_grid(.~Harvest.day2, scales="free_y") + guides(fill=guide_legend(nrow=2,byrow=TRUE)) 

p

# The two lines we want on the plot
fig_co2.native3 <- p
# export figure
ggpubr::ggexport(fig_co2.native3, height=1500, width=2500, filename = "figures/_uv/co2.native_2way.png", res = 400)



testlet_co2.native_res <- rbind(testlet_co2.native_24hr_res, testlet_co2.native_6mo_res)
testlet_co2.native_res$Harvest.day2 <- as.factor(c(rep("24 hr", 5), rep("6 mo", 5)))
testlet_co2.native_res$minys <- testlet_co2.native_res$response-testlet_co2.native_res$SE
# testlet_co2.native$minys[which(testlet_co2.native$minys<0)] <- 0

#testlet_co2.native_res$.group[1:5] <- rep(NA, 5)

# plot
testlet_co2.native_res$y_position <- c(rep(0, 5), rep(0, 5)) + testlet_co2.native_res$response + testlet_co2.native_res$SE + 5
segs <- data.frame(x=rep(1:5,1)-0.3, xend=rep(1:5,2)+0.3,
                   y=c(15,15,15,15,15)-0.7, yend=c(15,15,15,15,15)-0.7,
                   Harvest.day2="24 hr")
segs$Harvest.day2 <- as.factor(segs$Harvest.day2)
  geom_segment(inherit.aes=FALSE, data=segs, aes(x=x, y=y, yend=yend, xend=xend)) +
p <- ggplot(data=testlet_co2.native, aes(x=Residue.type, y=response, fill=Moisture.treatment)) +
  ggh4x::facet_nested(~Harvest.day2)+ # , scales = "free_y", independent="y"
  theme_minimal() + labs(fill="Moisture treatment") +
  ylab(label=expression(paste("Priming (", mu,"g g"^-1,")"))) +   
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_text(aes(label=c(rep(NA,10), testlet_co2.native$.group[11:20]), y=10+c(rep(0,10), testlet_co2.native$response[11:20])+c(rep(0,10), testlet_co2.native$SE[11:20])), 
            position = position_dodge(0.8), size=2) +
  # geom_text(inherit.aes=FALSE, data=testlet_co2.native,size=2.5,
  #           aes(label=testlet_co2.native$.group, x=rep(1:5, 2), y=testlet_co2.native$response+testlet_co2.native$SE))+
  geom_text(inherit.aes=FALSE, data=testlet_co2.native_res,size=2,
            aes(label=c(testlet_co2.native_res$.group[1:5],rep(NA,5)), x= rep(1:5,2), y=c(25,25,25,25,25, rep(0,5)))) +
  geom_segment(inherit.aes=FALSE, data=segs, aes(x=x, y=y, yend=yend, xend=xend)) +
  #scale_y_continuous(breaks=seq(-100,200, 50)) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.1, 0.5), "cm"), axis.title.x=element_blank(),
        legend.text = element_text(size = 8),  # Change text size
        legend.title = element_text(size = 9),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0)) +
        #panel.spacing=unit(0.5, "lines"),
        #axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))) 
  scale_y_continuous(limits = c(-70, 170), breaks = seq(-50, 150, by = 50))
  

p

# The two lines we want on the plot
fig_co2.native3b <- p
# export figure
ggpubr::ggexport(fig_co2.native3b, height=1800, width=2500, filename = "figures/_uv/co2.native_2wayb.png", res = 400)




# export cue and priming figs
fig_CUE4 <- fig_CUE3+
  theme(plot.margin = unit(c(0,0.2,0,0.5), "cm"), legend.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes = list(size=3), title.position = "top")) 
fig_co2.native4 <- fig_co2.native3+
  theme(plot.margin = unit(c(0,0.2,0,0.2), "cm"), legend.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes = list(size=3), title.position = "top")) 

  
fig <- ggpubr::ggarrange(fig_CUE4,fig_co2.native4, nrow=2, labels = "AUTO",
                         common.legend = TRUE, legend="bottom")

ggpubr::ggexport(fig, height=2500, width=2500, 
         filename = "figures/_uv/1_cue and co2.native_2way.png", res = 400)


# export cue and priming figs
fig_CUE4b <- fig_CUEb+
  theme(plot.margin = unit(c(0,0.2,0,0.5), "cm"), legend.title = element_text(hjust = 0.5),
        axis.text.x = element_blank())+
  guides(fill = guide_legend(override.aes = list(size=3), title.position = "top")) 
fig_co2.native4b <- fig_co2.native3b+
  theme(plot.margin = unit(c(0,0.2,0,0.2), "cm"), legend.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes = list(size=3), title.position = "top")) 


fig <- ggpubr::ggarrange(fig_CUE4b, fig_co2.native4b, nrow=2, labels = "AUTO",
                         common.legend = TRUE, legend="bottom", heights=c(0.75,1))

ggpubr::ggexport(fig, height=2500, width=2500, 
                 filename = "figures/_uv/1_cue and co2.native_2wayb.png", res = 400)

# export CUE_r and priming figs
fig_CUE_rb2 <- fig_CUE_rb+
  theme(plot.margin = unit(c(0,0.2,0,0.5), "cm"), legend.title = element_text(hjust = 0.5),
        axis.text.x = element_blank())+
  guides(fill = guide_legend(override.aes = list(size=3), title.position = "top")) 
fig_co2.native4b <- fig_co2.native3b+
  theme(plot.margin = unit(c(0,0.2,0,0.2), "cm"), legend.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes = list(size=3), title.position = "top")) 


fig <- ggpubr::ggarrange(fig_CUE_rb2, fig_co2.native4b, nrow=2, labels = "AUTO",
                         common.legend = TRUE, legend="bottom", heights=c(0.75,1))

ggpubr::ggexport(fig, height=2500, width=2500, 
                 filename = "figures/_uv/1_CUE_r and co2.native_2wayb.png", res = 400)


# export CSE and priming figs
segs2 <- data.frame(x=rep(1:5,1)-0.3, xend=rep(1:5,2)+0.3,
                   y=rep(0.78, 5), yend=rep(0.78,5),
                   Harvest.day2="24 hr")
segs2$Harvest.day2 <- as.factor(segs2$Harvest.day2)

fig_CSEb2 <- fig_CSEb+
  theme(plot.margin = unit(c(0,0.2,0,0.5), "cm"), legend.title = element_text(hjust = 0.5),
        axis.text.x = element_blank())+
  geom_segment(inherit.aes=FALSE, data=segs2, aes(x=x, y=y, yend=yend, xend=xend)) +
  guides(fill = guide_legend(override.aes = list(size=3), title.position = "top")) 
fig_co2.native4b <- fig_co2.native3b+
  theme(plot.margin = unit(c(0,0.2,0,0.2), "cm"), legend.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes = list(size=3), title.position = "top")) 


fig <- ggpubr::ggarrange(fig_CSEb2, fig_co2.native4b, nrow=2, labels = "AUTO",
                         common.legend = TRUE, legend="bottom", heights=c(0.75,1))

ggpubr::ggexport(fig, height=2500, width=2500, 
                 filename = "figures/_uv/1_CSE and co2.native_2wayb.png", res = 400)








####### X13maoc
n <- 3
mod.anova2[1,n] <- "X13maoc_24 hr"
hist(dat$X13maoc[which(dat$Harvest.day2=="24 hr")])
fit <- lme(log(X13maoc) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$X13maoc)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(X13maoc),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/X13maoc_24 hr_2-way.csv")
testlet_X13maoc_24hr <- testlet

n <- 4
mod.anova2[1,n] <- "X13maoc_6 mo"
hist(dat$X13maoc[which(dat$Harvest.day2=="6 mo")])
fit <- lme((X13maoc) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$X13maoc)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(X13maoc),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/X13maoc_6 mo_2-way.csv")
testlet_X13maoc_6mo <- testlet
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/X13maoc_6 mo_2-way_res.csv")
testlet_X13maoc_6mo_res <- testlet

testlet_X13maoc <- rbind(testlet_X13maoc_24hr, testlet_X13maoc_6mo)
testlet_X13maoc$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_X13maoc$minys <- testlet_X13maoc$response-testlet_X13maoc$SE
testlet_X13maoc$minys[which(testlet_X13maoc$minys<0)] <- 0
testlet_X13maoc$.group <- rep(NA,20)

# plot
p <- ggplot(data=testlet_X13maoc, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("X13maoc"^""))) +   
  facet_grid(.~Harvest.day2) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_X13maoc$response+testlet_X13maoc$SE)+0.2*max(testlet_X13maoc$response+testlet_X13maoc$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), #legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=1+testlet_X13maoc$response+testlet_X13maoc$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_X13maoc3 <- p
# export figure
ggpubr::ggexport(fig_X13maoc3, height=1000, width=2500, filename = "figures/_uv/X13maoc_2way.png", res = 400)


testlet_X13maoc_res <- rbind(testlet_X13maoc_6mo_res, testlet_X13maoc_6mo_res)
testlet_X13maoc_res$Harvest.day2 <- as.factor(c(rep("24 hr", 5), rep("6 mo", 5)))
testlet_X13maoc_res$minys <- testlet_X13maoc_res$response-testlet_X13maoc_res$SE
testlet_X13maoc_res$.group[1:5] <- rep(NA,5)


# plot
p <- ggplot(data=testlet_X13maoc, aes(x=Residue.type, y=response, fill=Moisture.treatment)) +
  theme_minimal() + labs(fill="Moisture treatment") +
  ylab(expression(paste("X13maoc"^""))) +   
  facet_grid(.~Harvest.day2) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = c("blue", "red")) +
  ylim(c(0,max(testlet_X13maoc$response+testlet_X13maoc$SE)+0.2*max(testlet_X13maoc$response+testlet_X13maoc$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), #legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0)) +
  # geom_text(aes(label=.group, y=1+testlet_X13maoc$response+testlet_X13maoc$SE), 
  #           position = position_dodge(0.9), size=2.5) +
  geom_text(inherit.aes=FALSE, data=testlet_X13maoc_res,size=2.5,
            aes(label=testlet_X13maoc_res$.group, x=rep(1:5, 2), y=testlet_X13maoc_res$response+testlet_X13maoc_res$SE+2.5))
  
p

# The two lines we want on the plot
fig_X13maoc3b <- p
# export figure
ggpubr::ggexport(fig_X13maoc3b, height=1600, width=2500, filename = "figures/_uv/X13maoc_2way.png", res = 400)





####### X13co2
n <- 5
mod.anova2[1,n] <- "X13co2_24 hr"
hist(dat$X13co2[which(dat$Harvest.day2=="24 hr")])
fit <- lme((X13co2) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$X13co2)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(X13co2),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/X13co2_24 hr_2-way.csv")
testlet_X13co2_24hr <- testlet
### tukey: 
test1 <- emmeans(fit, ~Moisture.treatment)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Moisture.treatment),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/X13co2_24 hr_2-way_mois.csv")
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/X13co2_24 hr_2-way_res.csv")

n <- 6
mod.anova2[1,n] <- "X13co2_6 mo"
hist(dat$X13co2[which(dat$Harvest.day2=="6 mo")])
fit <- lme((X13co2) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$X13co2)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(X13co2),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/X13co2_6 mo_2-way.csv")
testlet_X13co2_6mo <- testlet
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/X13co2_6 mo_2-way_res.csv")


testlet_X13co2 <- rbind(testlet_X13co2_24hr, testlet_X13co2_6mo)
testlet_X13co2$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_X13co2$minys <- testlet_X13co2$response-testlet_X13co2$SE
#testlet_X13co2$minys[which(testlet_X13co2$minys<0)] <- 0


# plot
p <- ggplot(data=testlet_X13co2, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("X13co2"^""))) +   
  facet_grid(.~Harvest.day2) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_X13co2$response+testlet_X13co2$SE)+0.2*max(testlet_X13co2$response+testlet_X13co2$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), #legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=1+testlet_X13co2$response+testlet_X13co2$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_X13co23 <- p
# export figure
ggpubr::ggexport(fig_X13co23, height=1000, width=2500, filename = "figures/_uv/X13co2_2way.png", res = 400)











####### X13poc
n <- 7
mod.anova2[1,n] <- "X13poc_24 hr"
hist(dat$X13poc[which(dat$Harvest.day2=="24 hr")])
fit <- lme((X13poc) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$X13poc)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(X13poc),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/X13poc_24 hr_2-way.csv")
testlet_X13poc_24hr <- testlet

n <- 8
mod.anova2[1,n] <- "X13poc_6 mo"
hist(dat$X13poc[which(dat$Harvest.day2=="6 mo")])
dat[which(dat$Harvest.day2=="6 mo" & dat$X13poc>4),]
fit <- lme((X13poc) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$X13poc)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(X13poc),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/X13poc_6 mo_2-way.csv")
testlet_X13poc_6mo <- testlet
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/X13poc_6 mo_2-way_res.csv")


testlet_X13poc <- rbind(testlet_X13poc_24hr, testlet_X13poc_6mo)
testlet_X13poc$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_X13poc$minys <- testlet_X13poc$response-testlet_X13poc$SE
testlet_X13poc$minys[which(testlet_X13poc$minys<0)] <- 0

testlet_X13poc$.group[1:10] <- rep(NA,10)

# plot
p <- ggplot(data=testlet_X13poc, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("X13poc"^""))) +   
  facet_grid(.~Harvest.day2) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_X13poc$response+testlet_X13poc$SE)+0.2*max(testlet_X13poc$response+testlet_X13poc$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), #legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=0.5+testlet_X13poc$response+testlet_X13poc$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_X13poc3 <- p
# export figure
ggpubr::ggexport(fig_X13poc3, height=1000, width=2500, filename = "figures/_uv/X13poc_2way.png", res = 400)









 ####### X13mbc
 n <- 9
 mod.anova2[1,n] <- "X13mbc_24 hr"
 hist(dat$X13mbc[which(dat$Harvest.day2=="24 hr")])
 fit <- lme((X13mbc) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$X13mbc)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(X13mbc),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
 write.csv(testlet, "model-output/_uv/lsmeans/X13mbc_24 hr_2-way.csv")
 testlet_X13mbc_24hr <- testlet
 
 n <- 10
 mod.anova2[1,n] <- "X13mbc_6 mo"
 hist(dat$X13mbc[which(dat$Harvest.day2=="6 mo")])
 fit <- lme((X13mbc) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$X13mbc)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(X13mbc),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
 write.csv(testlet, "model-output/_uv/lsmeans/X13mbc_6 mo_2-way.csv")
 testlet_X13mbc_6mo <- testlet
 
 
 testlet_X13mbc <- rbind(testlet_X13mbc_24hr, testlet_X13mbc_6mo)
 testlet_X13mbc$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
 testlet_X13mbc$minys <- testlet_X13mbc$response-testlet_X13mbc$SE
 testlet_X13mbc$minys[which(testlet_X13mbc$minys<0)] <- 0
 testlet_X13mbc$.group[1:10] <- rep(NA, 10)
 
 # plot
 p <- ggplot(data=testlet_X13mbc, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
   theme_minimal() + labs(fill="Residue type") +
   ylab(expression(paste("X13mbc"^""))) +   
   facet_grid(.~Harvest.day2) +
   geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
   geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                 width=0.1, position=position_dodge(0.9)) +
   scale_fill_manual(values = covercols) +
   ylim(c(0,max(testlet_X13mbc$response+testlet_X13mbc$SE)+0.2*max(testlet_X13mbc$response+testlet_X13mbc$SE))) +
   theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), #legend.position = "none",
         panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
         panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
         plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
   geom_text(aes(label=.group, y=0.5+testlet_X13mbc$response+testlet_X13mbc$SE), 
             position = position_dodge(0.9), size=2.5)
 p
 
 # The two lines we want on the plot
 fig_X13mbc3 <- p
 # export figure
 ggpubr::ggexport(fig_X13mbc3, height=1000, width=2500, filename = "figures/_uv/X13mbc_2way.png", res = 400)
 




 
 
 ####### X13doc
 n <- 11
 mod.anova2[1,n] <- "X13doc_24 hr"
 dat$X13doc[which(dat$X13doc<0)] <-0
 hist(dat$X13doc[which(dat$Harvest.day2=="24 hr")])
 fit <- lme((X13doc) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$X13doc)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(X13doc),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
 write.csv(testlet, "model-output/_uv/lsmeans/X13doc_24 hr_2-way.csv")
 testlet_X13doc_24hr <- testlet
 ### tukey: 
 test1 <- emmeans(fit, ~Moisture.treatment)
 testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
 testlet <- testlet[order(testlet$Moisture.treatment),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/X13doc_24 hr_2-way_mois.csv")
 ### tukey: 
 test1 <- emmeans(fit, ~Residue.type)
 testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
 testlet <- testlet[order(testlet$Residue.type),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/X13doc_24 hr_2-way_res.csv")
 
 n <- 12
 mod.anova2[1,n] <- "X13doc_6 mo"
 hist(dat$X13doc[which(dat$Harvest.day2=="6 mo")])
 fit <- lme((X13doc) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$X13doc)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(X13doc),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
 write.csv(testlet, "model-output/_uv/lsmeans/X13doc_6 mo_2-way.csv")
 testlet_X13doc_6mo <- testlet
 ### tukey: 
 test1 <- emmeans(fit, ~Residue.type)
 testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
 testlet <- testlet[order(testlet$Residue.type),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/X13doc_6 mo_2-way_res.csv")
 
 
 
 testlet_X13doc <- rbind(testlet_X13doc_24hr, testlet_X13doc_6mo)
 testlet_X13doc$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
 testlet_X13doc$minys <- testlet_X13doc$response-testlet_X13doc$SE
 #testlet_X13doc$minys[which(testlet_X13doc$minys<0)] <- 0
 testlet_X13doc$.group[11:20] <- rep(NA, 10)
 
 # plot
 p <- ggplot(data=testlet_X13doc, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
   theme_minimal() + labs(fill="Residue type") +
   ylab(expression(paste("X13doc"^""))) +   
   facet_grid(.~Harvest.day2) +
   geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
   geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                 width=0.1, position=position_dodge(0.9)) +
   scale_fill_manual(values = covercols) +
   ylim(c(0,max(testlet_X13doc$response+testlet_X13doc$SE)+0.2*max(testlet_X13doc$response+testlet_X13doc$SE))) +
   theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), #legend.position = "none",
         panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
         panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
         plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
   geom_text(aes(label=.group, y=0.5+testlet_X13doc$response+testlet_X13doc$SE), 
             position = position_dodge(0.9), size=2.5)
 p
 
 # The two lines we want on the plot
 fig_X13doc3 <- p
 # export figure
 ggpubr::ggexport(fig_X13doc3, height=1000, width=2500, filename = "figures/_uv/X13doc_2way.png", res = 400)
 

 
 
 ####### X13residue -- glucose in SOC minus glucose in MBC
 n <- 31
 mod.anova2[1,n] <- "X13residue_24 hr"
 dat$X13residue <- dat$X13soc - dat$X13mbc
 hist(dat$X13residue[which(dat$Harvest.day2=="24 hr")])
 fit <- lme((X13residue) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$X13residue)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(X13residue),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
 write.csv(testlet, "model-output/_uv/lsmeans/X13residue_24 hr_2-way.csv")
 testlet_X13residue_24hr <- testlet
 ### tukey: 
 test1 <- emmeans(fit, ~Moisture.treatment)
 testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
 testlet <- testlet[order(testlet$Moisture.treatment),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/X13residue_24 hr_2-way_mois.csv")
 ### tukey: 
 test1 <- emmeans(fit, ~Residue.type)
 testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
 testlet <- testlet[order(testlet$Residue.type),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/X13residue_24 hr_2-way_res.csv")
 testlet_X13residue_24hr_res <- testlet
 
 # custom contrast:
 test1
 A = c(0, 0, 0, 0, 1)
 B = c(0.25, 0.25, 0.25, 0.25, 0)
 con <- contrast(test1, method = list("A - B" = A - B)) 
 write.csv(con, "model-output/_uv/lsmeans/X13residue_24 hr_2-way_res_custom contrast.csv")
 
 
 
 n <- 32
 mod.anova2[1,n] <- "X13residue_6 mo"
 hist(dat$X13residue[which(dat$Harvest.day2=="6 mo")])
 dat$X13residue[which(dat$X13residue<0)] <- NA
 fit <- lme((X13residue) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$X13residue)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(X13residue),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
 write.csv(testlet, "model-output/_uv/lsmeans/X13residue_6 mo_2-way.csv")
 testlet_X13residue_6mo <- testlet
 ### tukey: 
 test1 <- emmeans(fit, ~Moisture.treatment)
 testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
 testlet <- testlet[order(testlet$Moisture.treatment),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/X13residue_6 mo_2-way_mois.csv")
 ### tukey: 
 test1 <- emmeans(fit, ~Residue.type)
 testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
 testlet <- testlet[order(testlet$Residue.type),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/X13residue_6 mo_2-way_res.csv")
 testlet_X13residue_6mo_res <- testlet
 
 
 testlet_X13residue <- rbind(testlet_X13residue_24hr, testlet_X13residue_6mo)
 testlet_X13residue$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
 testlet_X13residue$minys <- testlet_X13residue$response-testlet_X13residue$SE
 #testlet_X13residue$minys[which(testlet_X13residue$minys<0)] <- 0
 
 testlet_X13residue_res <- rbind(testlet_X13residue_24hr_res, testlet_X13residue_6mo_res)
 testlet_X13residue_res$Harvest.day2 <- as.factor(c(rep("24 hr", 5), rep("6 mo", 5)))
 testlet_X13residue_res$minys <- testlet_X13residue_res$response-testlet_X13residue_res$SE
 #testlet_X13residue$minys[which(testlet_X13residue$minys<0)] <- 0
 
 
 # plot
 p <- ggplot(data=testlet_X13residue, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
   theme_minimal() + labs(fill="Residue type") +
   ylab(expression(paste(""^13,"C in microbial residue (", mu,"g g"^-1,")"))) +   
   facet_wrap(.~Harvest.day2, scales="free_y") + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
   geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
   geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                 width=0.1, position=position_dodge(0.9)) +
   scale_fill_manual(values = covercols) +
   #ylim(c(0,max(testlet_X13residue$response+testlet_X13residue$SE)+0.2*max(testlet_X13residue$response+testlet_X13residue$SE))) +
   theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
         panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
         panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
         plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) #+
 # geom_text(inherit.aes = F, data = testlet_X13residue_res, 
 #           aes(label = c(testlet_X13residue_res$.group[1:5], rep(NA,5)), 
 #               y = 0.003+c(testlet_X13residue_res$response[1:5], rep(0,5)) + c(testlet_X13residue_res$SE[1:5], rep(0,5)),
 #               x = rep(seq(from = 1, to = 2, length.out = 5), 2)), 
 #           position = position_dodge(0.9), size=2)
 p
 
 # The two lines we want on the plot
 fig_X13residue3 <- p
 # export figure
 ggpubr::ggexport(fig_X13residue3, height=1300, width=2600, filename = "figures/_uv/1_X13residue_2way.png", res = 400)
 
 
 # plot
 p <- ggplot(data=testlet_X13residue, aes(x=Residue.type, y=response, fill=Moisture.treatment)) +
   theme_minimal() + labs(fill="Moisture treatment") +
   ylab(expression(paste(""^13,"C in microbial residue (", mu,"g g"^-1,")"))) +   
   facet_wrap(.~Harvest.day2) + guides(fill=guide_legend(nrow=1,byrow=TRUE)) + # , scales="free_y"
   geom_bar(stat = "identity", position = position_dodge(0.8), color="black", width=0.8) + 
   geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                 width=0.1, position=position_dodge(0.8)) +
   scale_fill_manual(values = c("blue", "red")) +
   #ylim(c(0,max(testlet_X13residue$response+testlet_X13residue$SE)+0.2*max(testlet_X13residue$response+testlet_X13residue$SE))) +
   theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
         panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
         panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
         plot.margin = unit(c(0.3, 0.3, 0.1, 0.5), "cm"), axis.title.x=element_blank(),
         axis.text.x = element_text(angle = 290, vjust = 0, hjust=0)) +
   geom_text(inherit.aes = F, data = testlet_X13residue_res, 
             aes(label = c(testlet_X13residue_res$.group[1:5], rep(NA,5)), 
                 y = c(rep(28, 5), rep(0, 5)), # 0.06+c(testlet_X13residue_res$response[1:5], rep(0,5)) + c(testlet_X13residue_res$SE[1:5], rep(0,5)),
                 x = rep(seq(from = 1, to = 5, length.out = 5), 2)), 
             position = position_dodge(0.9), size=2)
 p
 
 # The two lines we want on the plot
 fig_X13residueb <- p
 # export figure
 ggpubr::ggexport(fig_X13residueb, height=1700, width=2700, filename = "figures/_uv/1_X13residue_2way_b.png", res = 400)
 



##### FIGURE 2a: stacked barplot of residue-derived CO2, DOC, MBC, POC, and MAOC across residue type, moisture, and day


# co2, poc, maoc
testlet_X13co2$condition <- "13co2"
testlet_X13poc$condition <- "13poc"
testlet_1 <- bind_rows(testlet_X13co2, testlet_X13poc)
testlet_X13maoc$condition <- "13maoc"
testlet_4 <- bind_rows(testlet_1, testlet_X13maoc)

testlet_4$condition <- as.factor(testlet_4$condition)
testlet_4$condition <- factor(testlet_4$condition, levels(testlet_4$condition)[c(1,3,2)])
  
# Stacked
library('tibble')

testlet_X13maoc <- testlet_X13maoc[order(testlet_X13maoc$Harvest.day2, testlet_X13maoc$Moisture.treatment, testlet_X13maoc$Residue.type),]
testlet_X13poc <- testlet_X13poc[order(testlet_X13poc$Harvest.day2, testlet_X13poc$Moisture.treatment, testlet_X13maoc$Residue.type),]
testlet_X13co2 <- testlet_X13co2[order(testlet_X13co2$Harvest.day2, testlet_X13co2$Moisture.treatment, testlet_X13maoc$Residue.type),]


# recovery of glucose-C in CO2, MAOC, and POC
recovery <- data.frame(recovery = testlet_X13maoc$response + testlet_X13poc$response + testlet_X13co2$response)
recovery <- cbind(testlet_X13maoc, recovery)
write.csv(recovery, "model-output/recovery of glucose-c in CO2, MAOC, POC.csv")

stacked_fig <- ggplot(testlet_4, aes(fill=condition, y=response, x=Residue.type)) + 
  theme_bw() + 
  #lims(y=c(0,50)) +
  geom_bar(position="stack", stat="identity", color="black") +
  ggh4x::facet_nested(.~Harvest.day2+Moisture.treatment) +
  #labs(x="Residue type", fill=expression(paste(""^13,"C pool")), y="") +
  labs(x="Residue type", fill="", y="") +
  scale_fill_manual(values=c("red2", "darkolivegreen4",  "navy"), labels=c("Respiration", "Particulate organic matter",  "Mineral-associated organic matter")) +
  theme(panel.spacing=unit(1,"lines"),axis.ticks.x=element_blank(),axis.ticks.y=element_line(), 
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0), legend.position="top", plot.margin = unit(c(0,1,0,0), "cm")) +
  # add maoc letters
  geom_text(inherit.aes=FALSE, data=testlet_X13maoc, color="white", size=3,fontface=2,
            aes(label=testlet_X13maoc$.group, x=rep(1:5, 4), y=testlet_X13maoc$response/2))+
  # # add mbc letters
  # geom_text(inherit.aes=FALSE, data=testlet_X13mbc, color="white", size=3,fontface=2,
  #           aes(label=testlet_X13mbc$.group, x=rep(1:5, 4), y=(testlet_X13mbc$response/2)+testlet_X13maoc$response+c(rep(0,10),rep(-0.1,10))))+
  # add poc letters
  # geom_text(inherit.aes=FALSE, data=testlet_X13poc, color="white", size=3, fontface=2,
  #           aes(label=testlet_X13poc$.group, x=rep(1:5, 4), y=(testlet_X13poc$response/2)+testlet_X13maoc$response)) +
  # # add doc letters
  # geom_text(inherit.aes=FALSE, data=testlet_X13doc, color="white", size=3,fontface=2,
  #           aes(label=testlet_X13doc$.group, x=rep(1:5, 4), y=(testlet_X13doc$response/2)+testlet_X13poc$response+testlet_X13mbc$response+testlet_X13maoc$response+c(rep(0,10),rep(0.1,10)))) +
  # add co2 letters
  geom_text(inherit.aes=FALSE, data=testlet_X13co2, color="white", size=3, fontface=2,
            aes(label=testlet_X13co2$.group, x=rep(1:5, 4), y=(testlet_X13co2$response/2)+testlet_X13poc$response+testlet_X13maoc$response)) +
  scale_y_continuous(name=expression(paste("Glucose-derived carbon (", mu,"g g"^-1,")")),
                     sec.axis = sec_axis(trans=~.*2, breaks = c(0,20,40,60,80), name=expression(paste("Glucose-derived carbon (%)"))))
stacked_fig
ggpubr::ggexport(stacked_fig, height=2000, width=4000, filename = "figures/_uv/2a_stacked.png", res = 400)





# co2, poc, maoc
testlet_X13co2$condition <- "13co2"
testlet_X13poc$condition <- "13poc"
testlet_1 <- bind_rows(testlet_X13co2, testlet_X13poc)
testlet_X13maoc$condition <- "13maoc"
testlet_4 <- bind_rows(testlet_1, testlet_X13maoc)

testlet_4$condition <- as.factor(testlet_4$condition)
testlet_4$condition <- factor(testlet_4$condition, levels(testlet_4$condition)[c(1,3,2)])

# Stacked
library('tibble')
testlet_X13maoc <- testlet_X13maoc[order(testlet_X13maoc$Harvest.day2,testlet_X13maoc$Residue.type),]
testlet_X13poc <- testlet_X13poc[order(testlet_X13poc$Harvest.day2,testlet_X13poc$Residue.type),]
testlet_X13co2 <- testlet_X13co2[order(testlet_X13co2$Harvest.day2,testlet_X13co2$Residue.type),]

segs <- data.frame(x=c(rep(NA, 5), rep(1.5,5))-0.4, xend=c(rep(NA, 5), rep(1.5,5))+0.4,
                   y=c(rep(NA, 5), testlet_X13maoc_res$response[6:10]/2-2), yend=c(rep(NA, 5), testlet_X13maoc_res$response[6:10]/2-2),
                   Harvest.day2=rep(levels(testlet_4$Harvest.day2), each=5), 
                   Residue.type=rep(levels(testlet_4$Residue.type), 2))
segs$Harvest.day2 <- as.factor(segs$Harvest.day2)
segs$Residue.type <- as.factor(segs$Residue.type)
segs$Residue.type <- factor(segs$Residue.type, levels(segs$Residue.type)[c(3,4,1,5,2)])

stacked_fig_b <- ggplot(testlet_4, aes(fill=condition, y=response, x=Moisture.treatment)) + 
  theme_bw() + 
  #lims(y=c(0,50)) +
  geom_bar(position="stack", stat="identity", color="black") +
  ggh4x::facet_nested(.~Harvest.day2+Residue.type) +
  #labs(x="Residue type", fill=expression(paste(""^13,"C pool")), y="") +
  labs(x="Moisture treatment", fill="", y="") +
  scale_fill_manual(values=c( "firebrick2", "olivedrab", "navy"), labels=c("Respiration", "Particulate organic matter",  "Mineral-associated organic matter")) +
  theme(panel.spacing=unit(1,"lines"),axis.ticks.x=element_blank(),axis.ticks.y=element_line(), 
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0), legend.position="top", 
        plot.margin = unit(c(0,1,0,0), "cm"), strip.text.x = element_text(size = 8)) +
  # add maoc letters
  geom_rect(inherit.aes=FALSE, data=segs, aes(xmin=x, ymin=y+1, ymax=yend+4, xmax=xend), fill="black", color="white") +
  geom_text(inherit.aes=FALSE, data=testlet_X13maoc_res, color="white", size=3,fontface=2,
            aes(label=testlet_X13maoc_res$.group, x=rep(1.5, 10), y=0.5+testlet_X13maoc_res$response/2))+
  # # add poc letters
  # geom_text(inherit.aes=FALSE, data=testlet_X13poc, color="white", size=3, fontface=2,
  #           aes(label=testlet_X13poc$.group, x=rep(1:2, 10), y=(testlet_X13poc$response/2)+testlet_X13maoc$response)) +
  # add co2 letters
  geom_text(inherit.aes=FALSE, data=testlet_X13co2, color="white", size=3, fontface=2,
            aes(label=testlet_X13co2$.group, x=rep(1:2, 10), y=(testlet_X13co2$response/2)+testlet_X13poc$response+testlet_X13maoc$response)) +
  scale_y_continuous(name=expression(paste("Glucose-derived carbon (", mu,"g g"^-1,")")),
                     sec.axis = sec_axis(trans=~.*2, breaks = c(0,20,40,60,80), name=expression(paste("Glucose-derived carbon (%)"))))
stacked_fig_b

ggpubr::ggexport(stacked_fig_b, height=2000, width=5000, filename = "figures/_uv/2a_stacked_b.png", res = 400)
















# create a dataset
testlet_X13doc$condition <- "13doc"
testlet_X13mbc$condition <- "13mbc"
testlet_1 <- bind_rows(testlet_X13doc, testlet_X13mbc)
testlet_4 <- testlet_1

testlet_4$condition <- as.factor(testlet_4$condition)
testlet_4$condition <- factor(testlet_4$condition, levels(testlet_4$condition)[c(1,2)])

testlet_X13mbc <- testlet_X13mbc[order(testlet_X13mbc$Harvest.day2, testlet_X13mbc$Moisture.treatment),]
testlet_X13doc <- testlet_X13doc[order(testlet_X13doc$Harvest.day2, testlet_X13doc$Moisture.treatment),]

stacked_fig <- ggplot(testlet_4, aes(fill=condition, y=response, x=Residue.type)) + 
  theme_bw() + 
  lims(y=c(0,15.5)) +
  geom_bar(position="stack", stat="identity", color="black") +
  ggh4x::facet_nested(.~Harvest.day2+Moisture.treatment) + # , scales="free_y", independent="y"
#  labs(x="Residue type", fill=expression(paste(""^13,"C pool")), y="") +
  labs(x="Residue type", fill="Carbon pool", y="") +
  scale_fill_manual(values=c("skyblue", "green3"), labels=c("Dissolved organic matter", "Microbial biomass")) +
  theme(panel.spacing=unit(1,"lines"),axis.ticks.x=element_blank(),axis.ticks.y=element_line(), 
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0), legend.position="top", plot.margin = unit(c(0,1,0,0), "cm")) +
  # # add maoc letters
  # geom_text(inherit.aes=FALSE, data=testlet_X13maoc, color="white", size=3,fontface=2,
  #           aes(label=testlet_X13maoc$.group, x=rep(1:5, 4), y=testlet_X13maoc$response/2))+
  # add mbc letters
  # geom_text(inherit.aes=FALSE, data=testlet_X13mbc, color="black", size=2.3,fontface=1,
  #           aes(label=rep(NA,20), x=rep(1:5, 4), y=rep(c(15,3), each=10)))+
  geom_text(inherit.aes=FALSE, data=testlet_X13mbc, color="black", size=2.3,fontface=1,
            aes(label=testlet_X13mbc$.group, x=rep(1:5, 4), y=(testlet_X13mbc$response/2)))+
  # # add poc letters
  # geom_text(inherit.aes=FALSE, data=testlet_X13poc, color="white", size=3, fontface=2,
  #           aes(label=testlet_X13poc$.group, x=rep(1:5, 4), y=(testlet_X13poc$response/2)+testlet_X13mbc$response+testlet_X13maoc$response)) +
  # add doc letters
  geom_text(inherit.aes=FALSE, data=testlet_X13doc, color="black", size=2.3,fontface=1,
            aes(label=testlet_X13doc$.group, x=rep(1:5, 4), y=(testlet_X13doc$response/2)+testlet_X13mbc$response)) +
  # # add co2 letters
  # geom_text(inherit.aes=FALSE, data=testlet_X13co2, color="white", size=3, fontface=2,
  #           aes(label=testlet_X13co2$.group, x=rep(1:5, 4), y=(testlet_X13co2$response/2)+testlet_X13doc$response+testlet_X13poc$response+testlet_X13mbc$response+testlet_X13maoc$response))
  scale_y_continuous(name=expression(paste("Glucose-derived carbon (", mu,"g g"^-1,")")), breaks = c(0,3,6,9,12,15), 
                     sec.axis = sec_axis(trans=~.*2, breaks = c(0,6,12,18,24,30), name=expression(paste("Glucose-derived carbon (%)"))))

stacked_fig
ggpubr::ggexport(stacked_fig, height=1700, width=4000, filename = "figures/_uv/2a-2_stacked.png", res = 400)

testlet_X13mbc <- testlet_X13mbc[order(testlet_X13mbc$Harvest.day2, testlet_X13mbc$Residue.type, testlet_X13mbc$Moisture.treatment),]
testlet_X13doc <- testlet_X13doc[order(testlet_X13doc$Harvest.day2, testlet_X13doc$Residue.type, testlet_X13doc$Moisture.treatment),]
stacked_figb <- ggplot(testlet_4, aes(fill=condition, y=response, x=Moisture.treatment)) + 
  theme_bw() + 
  lims(y=c(0,15.5)) +
  geom_bar(position="stack", stat="identity", color="black") +
  ggh4x::facet_nested(.~Harvest.day2+Residue.type, ) + # , scales="free_y", independent="y"
  #  labs(x="Residue type", fill=expression(paste(""^13,"C pool")), y="") +
  labs(x="Residue type", fill="Carbon pool", y="") +
  scale_fill_manual(values=c("skyblue", "green3"), labels=c("Dissolved organic matter", "Microbial biomass")) +
  theme(panel.spacing=unit(1,"lines"),axis.ticks.x=element_blank(),axis.ticks.y=element_line(), 
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0), 
        legend.position="top", plot.margin = unit(c(0,1,0,0), "cm"),
        text = element_text(size = 8)) +
  # # add maoc letters
  # geom_text(inherit.aes=FALSE, data=testlet_X13maoc, color="white", size=3,fontface=2,
  #           aes(label=testlet_X13maoc$.group, x=rep(1:5, 4), y=testlet_X13maoc$response/2))+
  # add mbc letters
  # geom_text(inherit.aes=FALSE, data=testlet_X13mbc, color="black", size=2.3,fontface=1,
  #           aes(label=rep(NA,20), x=rep(1:5, 4), y=rep(c(15,3), each=10)))+
  geom_text(inherit.aes=FALSE, data=testlet_X13mbc, color="black", size=2.3,fontface=1,
            aes(label=testlet_X13mbc$.group, x=rep(1:2, 10), y=(testlet_X13mbc$response/2)))+
  # # add poc letters
  # geom_text(inherit.aes=FALSE, data=testlet_X13poc, color="white", size=3, fontface=2,
  #           aes(label=testlet_X13poc$.group, x=rep(1:5, 4), y=(testlet_X13poc$response/2)+testlet_X13mbc$response+testlet_X13maoc$response)) +
  # add doc letters
  geom_text(inherit.aes=FALSE, data=testlet_X13doc, color="black", size=2.3,fontface=1,
            aes(label=testlet_X13doc$.group, x=rep(1:2, 10), y=(testlet_X13doc$response/2)+testlet_X13mbc$response)) +
  # # add co2 letters
  # geom_text(inherit.aes=FALSE, data=testlet_X13co2, color="white", size=3, fontface=2,
  #           aes(label=testlet_X13co2$.group, x=rep(1:5, 4), y=(testlet_X13co2$response/2)+testlet_X13doc$response+testlet_X13poc$response+testlet_X13mbc$response+testlet_X13maoc$response))
  scale_y_continuous(name=expression(paste("Glucose-derived carbon (", mu,"g g"^-1,")")), breaks = c(0,3,6,9,12,15), 
                     sec.axis = sec_axis(trans=~.*2, breaks = c(0,6,12,18,24,30), name=expression(paste("Glucose-derived carbon (%)"))))

stacked_figb
ggpubr::ggexport(stacked_figb, height=1700, width=4000, filename = "figures/_uv/2a-2_stackedb.png", res = 400)







####### X13soc
n <- 13
mod.anova2[1,n] <- "X13soc_24 hr"
#dat$X13soc[which(dat$X13soc<0)] <-0
hist(dat$X13soc[which(dat$Harvest.day2=="24 hr")])
fit <- lme((X13soc) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$X13soc)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(X13soc),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/X13soc_24 hr_2-way.csv")
testlet_X13soc_24hr <- testlet
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/X13soc_24 hr_2-way_res.csv")
testlet_X13soc_24hr_res <- testlet

n <- 14
mod.anova2[1,n] <- "X13soc_6 mo"
hist(dat$X13soc[which(dat$Harvest.day2=="6 mo")])
fit <- lme((X13soc) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
           random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
resid <- residuals(fit)
mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
# check variances: if Var is >5, use a different model that allows variances to differ across groups
resid2 = cbind(dat[which(is.na(dat$X13soc)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
  summarise_at(vars(X13soc),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
write.csv(testlet, "model-output/_uv/lsmeans/X13soc_6 mo_2-way.csv")
testlet_X13soc_6mo <- testlet
### tukey: 
test1 <- emmeans(fit, ~Residue.type)
testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
testlet <- testlet[order(testlet$Residue.type),]
testlet$.group <- gsub(" ", "", testlet$.group)
write.csv(testlet, "model-output/_uv/lsmeans/X13soc_6 mo_2-way_res.csv")
testlet_X13soc_6mo_res <- testlet


testlet_X13soc <- rbind(testlet_X13soc_24hr, testlet_X13soc_6mo)
testlet_X13soc$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
testlet_X13soc$minys <- testlet_X13soc$response-testlet_X13soc$SE
testlet_X13soc$minys[which(testlet_X13soc$minys<0)] <- 0
testlet_X13soc$.group <- rep(NA, 20)

# plot
p <- ggplot(data=testlet_X13soc, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
  theme_minimal() + labs(fill="Residue type") +
  ylab(expression(paste("X13soc"^""))) +   
  facet_grid(.~Harvest.day2) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
  geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                width=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values = covercols) +
  ylim(c(0,max(testlet_X13soc$response+testlet_X13soc$SE)+0.2*max(testlet_X13soc$response+testlet_X13soc$SE))) +
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), #legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
  geom_text(aes(label=.group, y=2+testlet_X13soc$response+testlet_X13soc$SE), 
            position = position_dodge(0.9), size=2.5)
p

# The two lines we want on the plot
fig_X13soc3 <- p
# export figure
ggpubr::ggexport(fig_X13soc3, height=1000, width=2500, filename = "figures/_uv/X13soc_2way.png", res = 400)




##### FIGURE 2c: stacked barplot of residue-derived CO2 and SOC


# create a dataset
testlet_X13co2$condition <- "13co2"
testlet_X13soc$condition <- "13soc"
testlet_5 <- bind_rows(testlet_X13co2, testlet_X13soc)

# Stacked
library('tibble')

testlet_X13co2 <- testlet_X13co2[order(testlet_X13co2$Harvest.day2, testlet_X13co2$Moisture.treatment),]

stacked_fig_c <- ggplot(testlet_5, aes(fill=condition, y=response, x=Residue.type)) + 
  theme_bw() + 
  lims(y=c(0,50)) +
  geom_bar(position="stack", stat="identity", color="black") +
  ggh4x::facet_nested(.~Harvest.day2+Moisture.treatment) +
  labs(x="Residue type", fill=expression(paste(""^13,"C pool")), y="") +
  scale_fill_manual(values=c("red2", "darkorange4"), labels=c("Respiration", "SOC")) +
  theme(panel.spacing=unit(1,"lines"),axis.ticks.x=element_blank(),axis.ticks.y=element_line(), 
        panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
        axis.text.x = element_text(angle = 290, vjust = 0, hjust=0)) +
  # add SOC letters
  geom_text(inherit.aes=FALSE, data=testlet_X13soc, color="white", size=3, fontface=2,
            aes(label=testlet_X13soc$.group, x=rep(1:5, 4), y=testlet_X13soc$response/2))+
  # add co2 letters
  geom_text(inherit.aes=FALSE, data=testlet_X13co2, color="white", size=3, fontface=2,
            aes(label=testlet_X13co2$.group, x=rep(1:5, 4), y=(testlet_X13co2$response/2)+testlet_X13soc$response))
stacked_fig_c
ggpubr::ggexport(stacked_fig_c, height=2000, width=4000, filename = "figures/_uv/2c_stacked.png", res = 400)






          

 
####### m1_3pool
 n <- 17
 mod.anova2[1,n] <- "m1_3pool"
 hist(dat$m1_3pool[which(dat$Harvest.day2=="6 mo")])
 fit <- lme((m1_3pool) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$m1_3pool)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(m1_3pool),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
 write.csv(testlet, "model-output/_uv/lsmeans/m1_3pool_2-way.csv")
 testlet_m1_3pool <- testlet
 
 
 ####### m2_3pool
 n <- 18
 mod.anova2[1,n] <- "m2_3pool"
 hist(dat$m2_3pool[which(dat$Harvest.day2=="6 mo")])
 fit <- lme((m2_3pool) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$m2_3pool)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(m2_3pool),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
 write.csv(testlet, "model-output/_uv/lsmeans/m2_3pool_2-way.csv")
 testlet_m2_3pool <- testlet
 
 
 ####### m3_3pool
 n <- 19
 mod.anova2[1,n] <- "m3_3pool"
 hist(dat$m3_3pool[which(dat$Harvest.day2=="6 mo")])
 fit <- lme((m3_3pool) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$m3_3pool)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(m3_3pool),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
 write.csv(testlet, "model-output/_uv/lsmeans/m3_3pool_2-way.csv")
 testlet_m3_3pool <- testlet
 
 
 
 ##### FIGURE 2a:  barplot of 13c pools
 
 
 # create a dataset
 testlet_m1_3pool$condition <- "Fast pool size (%)"
 testlet_m2_3pool$condition <- "Intermediate pool size (%)"
 testlet_1 <- bind_rows(testlet_m1_3pool, testlet_m2_3pool)
 testlet_m3_3pool$condition <- "Slow pool size (%)"
 testlet_2 <- bind_rows(testlet_1, testlet_m3_3pool)

 testlet_2$condition <- as.factor(testlet_2$condition)

 # Stacked
 library('tibble')
 
 fig <- ggplot(testlet_2, aes(fill=Residue.type, y=response, x=Residue.type)) + 
   theme_bw() + 
   #lims(y=c(0,50)) +
   geom_bar(position=position_dodge(0.9), stat="identity", color="black") +
   geom_errorbar(aes(ymin = response-SE, ymax = response+SE), 
                 width=0.1, position=position_dodge(0.9)) +
   ggh4x::facet_nested(~condition+Moisture.treatment, scales = "free_y", independent="y") +
   labs(x="", fill="", y="", title="A") +
   scale_fill_manual(values=covercols) +
   theme(panel.spacing=unit(1,"lines"),axis.ticks.x=element_blank(),axis.ticks.y=element_line(), 
         panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
         axis.text.x =element_blank(), legend.position = "none") +
   # add letters
   geom_text(color="black", size=2,fontface=1,
             aes(label=.group, x=rep(1:5, 6), y=response+SE+c(rep(2,10), rep(1, 10), rep(5,10)))) +
   # set y lims
   geom_text(color="black", size=0,aes(label=.group, x=rep(1:5, 6), y=c(rep(35,10), rep(25, 10), rep(70,10))))
 fig

 ggpubr::ggexport(fig, height=1100, width=5000, filename = "figures/_uv/5c_barplot-sizes.png", res = 400)
 

 
 
 
 ####### k1_3pool
 n <- 20
 mod.anova2[1,n] <- "k1_3pool"
 hist(dat$k1_3pool[which(dat$Harvest.day2=="6 mo")])
 fit <- lme((k1_3pool) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$k1_3pool)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(k1_3pool),funs(mean=mean(., na.rm=T),Std=sd(.)))
 mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
 # export model stats
 mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
 mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
 mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
 ### tukey: 
 test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
 testlet <- cld(test1, type = "response", Letters = "ABCDEFGH", reversed = TRUE)
 testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/k1_3pool_2-way.csv")
 testlet_k1_3pool <- testlet
 
 
 ####### k2_3pool
 n <- 21
 mod.anova2[1,n] <- "k2_3pool"
 hist(dat$k2_3pool[which(dat$Harvest.day2=="6 mo")])
 fit <- lme((k2_3pool) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$k2_3pool)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(k2_3pool),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
 write.csv(testlet, "model-output/_uv/lsmeans/k2_3pool_2-way.csv")
 testlet_k2_3pool <- testlet
 
 
 ####### k3_3pool
 n <- 22
 mod.anova2[1,n] <- "k3_3pool"
 hist(dat$k3_3pool[which(dat$Harvest.day2=="6 mo")])
 fit <- lme((k3_3pool) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$k3_3pool)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(k3_3pool),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
 write.csv(testlet, "model-output/_uv/lsmeans/k3_3pool_2-way.csv")
 testlet_k3_3pool <- testlet
 
 
 

 

 



 
 
 
 ####### maom.cn.change-- change in maom c:n with glucose addition
 dat$maom.cn.change <- rep(NA, dim(dat)[1])
 dat$maom.cn.change[which(dat$Glucose.addition.treatment=="G-50")] <- dat$maom.cn[which(dat$Glucose.addition.treatment=="G-50")]-dat$maom.cn[which(dat$Glucose.addition.treatment=="G-00")]
 n <- 23
 mod.anova2[1,n] <- "maom.cn.change_24 hr"
 hist(dat$maom.cn.change[which(dat$Harvest.day2=="24 hr")])
 fit <- lme((maom.cn.change) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$maom.cn.change)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(maom.cn.change),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
 write.csv(testlet, "model-output/_uv/lsmeans/maom.cn.change_24 hr_2-way.csv")
 testlet_maom.cn.change_24hr <- testlet
 ### tukey: 
 test1 <- emmeans(fit, ~Residue.type)
 testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
 testlet <- testlet[order(testlet$Residue.type),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/maom.cn.change_24 hr_2-way_res.csv")
 
 n <- 24
 mod.anova2[1,n] <- "maom.cn.change_6 mo"
 hist(dat$maom.cn.change[which(dat$Harvest.day2=="6 mo")])
 fit <- lme((maom.cn.change) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$maom.cn.change)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(maom.cn.change),funs(mean=mean(., na.rm=T),Std=sd(.)))
 mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
 # export model stats
 mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
 mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
 mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
 ### tukey: 
 test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
 testlet <- cld(test1, type = "response", Letters = "ABC", reversed = TRUE)
 testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/maom.cn.change_6 mo_2-way.csv")
 testlet_maom.cn.change_6mo <- testlet
 ### tukey: 
 test1 <- emmeans(fit, ~Residue.type)
 testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
 testlet <- testlet[order(testlet$Residue.type),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/maom.cn.change_6 mo_2-way_res.csv")
 
 
 
 testlet_maom.cn.change <- rbind(testlet_maom.cn.change_24hr, testlet_maom.cn.change_6mo)
 testlet_maom.cn.change$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
 testlet_maom.cn.change$minys <- testlet_maom.cn.change$response-testlet_maom.cn.change$SE
 # testlet_maom.cn.change$minys[which(testlet_maom.cn.change$minys<0)] <- 0
 
 #testlet_maom.cn.change$.group <- rep(NA, 20)
 
 
 # plot
 p <- ggplot(data=testlet_maom.cn.change, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
   theme_minimal() + labs(fill="Residue type") +
   ylab(label=expression(paste("Change in MAOM C:N"))) +   
   facet_grid(.~Harvest.day2) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
   geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
   geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                 width=0.1, position=position_dodge(0.9)) +
   scale_fill_manual(values = covercols) + 
   scale_y_continuous(breaks=seq(-1,1.5, 0.5)) +
   theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
         panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
         panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
         plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
   geom_text(aes(label=.group, y=0.1+testlet_maom.cn.change$response+testlet_maom.cn.change$SE), 
             position = position_dodge(0.9), size=2.5)
 p
 
 # The two lines we want on the plot
 fig_maom.cn.change3 <- p
 # export figure
 ggpubr::ggexport(fig_maom.cn.change3, height=1700, width=2500, filename = "figures/_uv/maom.cn.change_2way.png", res = 400)
 
 
 
 
 
 
 
 
 ####### pom.cn.change-- change in pom c:n with glucose addition
 dat$pom.cn.change <- rep(NA, dim(dat)[1])
 dat$pom.cn.change[which(dat$Glucose.addition.treatment=="G-50")] <- dat$pom.cn[which(dat$Glucose.addition.treatment=="G-50")]-dat$pom.cn[which(dat$Glucose.addition.treatment=="G-00")]
 n <- 25
 mod.anova2[1,n] <- "pom.cn.change_24 hr"
 hist(dat$pom.cn.change[which(dat$Harvest.day2=="24 hr")])
 fit <- lme((pom.cn.change) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="24 hr"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$pom.cn.change)==FALSE & dat$Harvest.day2=="24 hr"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(pom.cn.change),funs(mean=mean(., na.rm=T),Std=sd(.)))
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
 write.csv(testlet, "model-output/_uv/lsmeans/pom.cn.change_24 hr_2-way.csv")
 testlet_pom.cn.change_24hr <- testlet
 ### tukey: 
 test1 <- emmeans(fit, ~Residue.type)
 testlet <- cld(test1, type = "response", Letters = "ABCDE", reversed = TRUE)
 testlet <- testlet[order(testlet$Residue.type),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/pom.cn.change_24 hr_2-way_res.csv")

  
 n <- 26
 mod.anova2[1,n] <- "pom.cn.change_6 mo"
 hist(dat$pom.cn.change[which(dat$Harvest.day2=="6 mo")])
 fit <- lme((pom.cn.change) ~ Moisture.treatment*Residue.type, weights=varIdent(form=~1|Moisture.treatment*Residue.type),
            random=~1|Lab.rep, data=dat[which(dat$Harvest.day2=="6 mo"),], na.action="na.omit")
 resid <- residuals(fit)
 mod.anova2[2,n] <- round(shapiro.test(resid)$p.value, 3) # check for normality: if <0.05, transform the data
 # check variances: if Var is >5, use a different model that allows variances to differ across groups
 resid2 = cbind(dat[which(is.na(dat$pom.cn.change)==FALSE & dat$Harvest.day2=="6 mo"),],resid)
 temp = resid2%>% group_by(Moisture.treatment,Residue.type)%>% 
   summarise_at(vars(pom.cn.change),funs(mean=mean(., na.rm=T),Std=sd(.)))
 mod.anova2[3,n] <- round(max(temp$Std)/min(temp$Std),3) # check that max/min variation across groups is less than five-fold
 # export model stats
 mod.out <- as.data.frame(anova(fit, type = "marginal"))[2:4,] 
 mod.out2 <- mod.out %>% mutate_if(is.numeric, round, digits=3)
 mod.anova2[4:6,n] <- paste0("F=", mod.out2$`F-value`, ", Df=", mod.out2$numDF, ", P=", mod.out2$`p-value`)
 ### tukey: 
 test1 <- emmeans(fit, ~Moisture.treatment*Residue.type)
 testlet <- cld(test1, type = "response", Letters = "ABC", reversed = TRUE)
 testlet <- testlet[order(testlet$Moisture.treatment, testlet$Residue.type),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/pom.cn.change_6 mo_2-way.csv")
 testlet_pom.cn.change_6mo <- testlet
 ### tukey: 
 test1 <- emmeans(fit, ~Residue.type)
 testlet <- cld(test1, type = "response", Letters = "ABC", reversed = TRUE)
 testlet <- testlet[order( testlet$Residue.type),]
 testlet$.group <- gsub(" ", "", testlet$.group)
 write.csv(testlet, "model-output/_uv/lsmeans/pom.cn.change_6 mo_2-way_res.csv")
 
 
 testlet_pom.cn.change <- rbind(testlet_pom.cn.change_24hr, testlet_pom.cn.change_6mo)
 testlet_pom.cn.change$Harvest.day2 <- as.factor(c(rep("24 hr", 10), rep("6 mo", 10)))
 testlet_pom.cn.change$minys <- testlet_pom.cn.change$response-testlet_pom.cn.change$SE
 # testlet_pom.cn.change$minys[which(testlet_pom.cn.change$minys<0)] <- 0
 
 testlet_pom.cn.change$.group <- rep(NA, 20)
 
 
 # plot
 p <- ggplot(data=testlet_pom.cn.change, aes(x=Moisture.treatment, y=response, fill=Residue.type)) +
   theme_minimal() + labs(fill="Residue type") +
   ylab(label=expression(paste("Change in POM C:N"))) +   
   facet_grid(.~Harvest.day2) + guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
   geom_bar(stat = "identity", position = position_dodge(0.9), color="black", width=0.9) + 
   geom_errorbar(aes(ymin = minys, ymax = response+SE), 
                 width=0.1, position=position_dodge(0.9)) +
   scale_fill_manual(values = covercols) + 
   scale_y_continuous(breaks=seq(-3,12, 3)) +
   theme(axis.ticks.x=element_blank(),axis.ticks.y=element_line(), legend.position = "bottom",
         panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
         panel.grid.minor = element_blank(), panel.grid.major.x=element_blank(),
         plot.margin = unit(c(0.3, 0.1, 0.1, 0.5), "cm"), axis.title.x=element_blank()) +
   geom_text(aes(label=.group, y=0.1+testlet_pom.cn.change$response+testlet_pom.cn.change$SE), 
             position = position_dodge(0.9), size=2.5)
 p
 
 # The two lines we want on the plot
 fig_pom.cn.change3 <- p
 # export figure
 ggpubr::ggexport(fig_pom.cn.change3, height=1700, width=2500, filename = "figures/_uv/pom.cn.change_2way.png", res = 400)
 
 
 
 
 
 # export anova results
write.csv(mod.anova2, "model-output/_uv/anova_2way.csv")
