

library(dplyr)
library(ggplot2)

covercols <- c("chocolate", "goldenrod", "forestgreen", "darkcyan", "darkorchid")
mcols <- c("blue", "red")


# 2023-05-22: initiate phase 1 (combine soil with residues and bring to 60% WHC)
# 2023-05-29: initiate drought
# 2023-06-30: bring all soils to 60% WHC
# 2023-07-10: glucose addition to set 1
# 2023-07-11: harvest set 1 (24 hr)
# 2023-07-12: glucose addition to set 2
# 2024-01-23: harvest set 2 (6 mo)

# days on which respiration was measured
days <- read.csv("raw-data/respiration-measurement-days.csv")


### cumulative CO2 data
dat <- read.csv("raw-data/master-dataset-respiration.csv")
dat_glucose <- read.csv("raw-data/master-dataset-respired glucose.csv")


# transform emission data to long format
dat_long <- reshape2::melt(dat, measure.vars=8:201, 
                                     variable.name="Date", value.name="CO2Cumulative")
dat_glucose_long <- reshape2::melt(dat_glucose, measure.vars=7:200, 
                           variable.name="Date", value.name="CO2Glucose")
dat_long$CO2Cumulative_Glucose <- dat_glucose_long$CO2Glucose
dat_long$CO2Cumulative_Native <- dat_long$CO2Cumulative - dat_long$CO2Cumulative_Glucose
dat_long$Date <- as.numeric(gsub("X", "", dat_long$Date))
dat_long$Sample.number.Date <- paste(dat_long$Sample.number, dat_long$Date, sep="_")


### CO2 flux data
dat_flux <- read.csv("raw-data/master-dataset-respiration-flux.csv")
dat_flux_glucose <- read.csv("raw-data/master-dataset-respired glucose-flux.csv")

# transform emission data to long format
dat_flux_long <- reshape2::melt(dat_flux, measure.vars=8:201, 
                           variable.name="Date", value.name="CO2Flux")
dat_flux_glucose_long <- reshape2::melt(dat_flux_glucose, measure.vars=7:200, 
                                   variable.name="Date", value.name="CO2Flux_Glucose")
dat_flux_long$CO2Flux_Glucose <- dat_flux_glucose_long$CO2Flux_Glucose
dat_flux_long$CO2Flux_Native <- dat_flux_long$CO2Flux - dat_flux_long$CO2Flux_Glucose
dat_flux_long$Date <- as.numeric(gsub("X", "", dat_flux_long$Date))
dat_flux_long$Sample.number.Date <- paste(dat_flux_long$Sample.number, dat_flux_long$Date, sep="_")

### merge flux and cumulative data
dat <- merge(dat_long, dat_flux_long[,c("Sample.number.Date", "CO2Flux", "CO2Flux_Glucose", "CO2Flux_Native")], by="Sample.number.Date", all=T)

# set factors
dat$Residue.type <- as.factor(dat$Residue.type)
dat$Residue.type <- factor(dat$Residue.type, levels(dat$Residue.type)[c(2,5,1,4,3)])
dat$Moisture.treatment <- as.factor(dat$Moisture.treatment)
dat$Sample.number <- as.factor(dat$Sample.number)
dat$Date <- as.numeric(dat$Date)
# remove day 0
dat <- dat[-which(dat$Date==0),]


### calculate priming 
# first find the native SOC respired when no plastic is added (cumulative ug C g dry soil-1 day-1)
native_respired <- dat %>%
  filter(Glucose.addition.treatment=="G-00") %>%
  group_by(Date, Moisture.treatment, Residue.type) %>%
  summarise(CO2Flux_Native = mean(CO2Flux_Native, na.rm = TRUE), .groups = "drop")

# calculate the change in native SOC respired when plastic is added (cumulative ug C g dry soil-1 day-1)
dat$priming <- rep(NA, dim(dat)[1])
dat$priming_relative <- rep(NA, dim(dat)[1])
for(i in 1:dim(dat)[1]){
  dat$priming[i] <- dat$CO2Flux_Native[i] - native_respired$CO2Flux_Native[which(native_respired$Residue.type==dat$Residue.type[i] & native_respired$Moisture.treatment==dat$Moisture.treatment[i] & native_respired$Date==dat$Date[i])]
  denom <- native_respired$CO2Flux_Native[which(native_respired$Residue.type==dat$Residue.type[i] & native_respired$Moisture.treatment==dat$Moisture.treatment[i] & native_respired$Date==dat$Date[i])]
  dat$priming_relative[i] <- (dat$CO2Flux_Native[i] - native_respired$CO2Flux_Native[which(native_respired$Residue.type==dat$Residue.type[i] & native_respired$Moisture.treatment==dat$Moisture.treatment[i] & native_respired$Date==dat$Date[i])])/denom
}
dat$priming[which(dat$Glucose.addition.treatment=="G-00")] <- NA
dat$priming_relative[which(dat$Glucose.addition.treatment=="G-00")] <- NA


dat <- dat[order(dat$Sample.number, dat$Date),]


# Calculate cumulative priming for each treatment combination
dat <- dat %>%
  group_by(Sample.number) %>%
  mutate(cumulative_priming = cumsum(priming)) %>%
  ungroup()

# calculate the proportion of glucose-C has been respired
dat$prop_gluc_c_resp <- dat$CO2Cumulative_Glucose/50
dat$prop_gluc_c_resp[which(dat$Residue.type=="NONE")] <- NA



# export data
write.csv(dat, "raw-data/CO2_Long.csv", row.names=FALSE)




##### Plot CO2 emission data: cumulative total CO2
# reorder levels of plastic
cumulative_data <- read.csv("raw-data/CO2_Long.csv")
# set factors
cumulative_data$Residue.type <- as.factor(cumulative_data$Residue.type)
cumulative_data$Residue.type <- factor(cumulative_data$Residue.type, levels(cumulative_data$Residue.type)[c(2,5,1,4,3)])
cumulative_data$Moisture.treatment <- as.factor(cumulative_data$Moisture.treatment)
cumulative_data$Sample.number <- as.factor(cumulative_data$Sample.number)
cumulative_data$Date <- as.numeric(cumulative_data$Date)

# Plot cumulative CO2 emissions for all cups
p <- ggplot(cumulative_data, aes(x = Date, y = CO2Cumulative, group=Sample.number, color = Residue.type)) +
  geom_line(aes(linetype=Moisture.treatment)) +
  scale_color_manual(values=covercols) +
  labs(x = "Days after glucose addition",
       y = bquote("Cumulative CO" [2] * " Emissions (" * mu * "g C g"^-1 * "dry soil)"),
       linetype = "Moisture treatment",
       color = "Residue type") +
  theme_bw()
p

#ggsave(p, file="figures/emissions/CO2-cumulative-sample.png",width = 6, height = 4, dpi = 300)


# Plot cumulative glucose for all cups
p <- ggplot(cumulative_data, aes(x = Date, y = CO2Cumulative_Glucose, group=Sample.number, color = Residue.type)) +
  geom_line(aes(linetype=Moisture.treatment)) +
  scale_color_manual(values=covercols) +
  labs(x = "Days after glucose addition",
       y = bquote("Glucose-derived CO" [2] * "-C Emissions (" * mu * "g C g"^-1 * "dry soil)"),
       linetype = "Moisture treatment",
       color = "Residue type") +
  theme_bw()
p

#ggsave(p, file="figures/emissions/priming-cumulative-sample.png",width = 6, height = 4, dpi = 300)


# Plot cumulative priming for all cups
p <- ggplot(cumulative_data, aes(x = Date, y = cumulative_priming, group=Sample.number, color = Residue.type)) +
  geom_line(aes(linetype=Moisture.treatment)) +
  scale_color_manual(values=covercols) +
  labs(x = "Days after glucose addition",
       y = bquote("Cumulative priming (" * mu * "g C g"^-1 * "dry soil)"),
       linetype = "Moisture treatment",
       color = "Residue type") +
  theme_bw()
p

#ggsave(p, file="figures/emissions/priming-cumulative-sample.png",width = 6, height = 4, dpi = 300)


# Summarize data to calculate treatment averages
cumulative_avg <- cumulative_data %>%
  group_by(Date, Residue.type, Moisture.treatment, Glucose.addition.treatment) %>%
  summarise(mean_CO2 = mean(CO2Cumulative, na.rm = TRUE),
            se_CO2 = sd(CO2Cumulative, na.rm = TRUE)/2,
            mean_CO2_flux = mean(CO2Flux, na.rm = TRUE),
            se_CO2_flux = sd(CO2Flux, na.rm = TRUE)/2,
            mean_CO2_glucose = mean(CO2Cumulative_Glucose, na.rm = TRUE),
            se_CO2_glucose = sd(CO2Cumulative_Glucose, na.rm = TRUE)/2,
            mean_CO2_glucose_flux = mean(CO2Flux_Glucose, na.rm = TRUE),
            se_CO2_glucose_flux = sd(CO2Flux_Glucose, na.rm = TRUE)/2,
            mean_CO2_glucose_prop = mean(prop_gluc_c_resp, na.rm = TRUE),
            se_CO2_glucose_prop = sd(prop_gluc_c_resp, na.rm = TRUE)/2,
           mean_CO2_native = mean(CO2Cumulative_Native, na.rm = TRUE),
            se_CO2_native = sd(CO2Cumulative_Native, na.rm = TRUE)/2,
           mean_CO2_native_flux = mean(CO2Flux_Native, na.rm = TRUE),
           se_CO2_native_flux = sd(CO2Flux_Native, na.rm = TRUE)/2,
           mean_priming = mean(cumulative_priming, na.rm = TRUE),
            se_priming = sd(cumulative_priming, na.rm = TRUE)/2, .groups = "drop")
cumulative_avg <- cumulative_avg[which(cumulative_avg$Date %in% days$Days),]

write.csv(cumulative_avg, "raw-data/CO2_Long_avg.csv")


# Plot cumulative CO2 emissions averages
p <- ggplot(cumulative_avg, aes(x = Date, y = mean_CO2, color = Residue.type, alpha=Glucose.addition.treatment)) +
  geom_errorbar(aes(ymin = mean_CO2 - se_CO2, ymax = mean_CO2 + se_CO2), width=2, size=0.5) +
  #geom_line(aes(linetype=Moisture.treatment), size=0.5) + 
  geom_line(size=0.5) + 
  facet_grid(.~Moisture.treatment) +
  scale_color_manual(values=covercols) +
  scale_alpha_manual(values=c(0.3,1)) +
  labs(x = "Days after glucose addition", y = bquote("Cumulative CO" [2] * " Emissions (" * mu * "g C g"^-1 * "dry soil)"),
       linetype = "Moisture treatment", alpha = "Glucose addition",
       color = "Residue type") +
  theme_bw()
p

ggsave(p, file="figures/emissions/CO2 cumulative.png",width = 8, height = 4, dpi = 300)


# Plot CO2 flux averages
p <- ggplot(cumulative_avg, aes(x = Date, y = mean_CO2_flux, color = Residue.type, alpha=Glucose.addition.treatment)) +
  geom_errorbar(aes(ymin = mean_CO2_flux - se_CO2_flux, ymax = mean_CO2_flux + se_CO2_flux), width=2, size=0.5) +
  #geom_line(aes(linetype=Moisture.treatment), size=0.5) + 
  geom_line(size=0.5) + 
  facet_grid(.~Moisture.treatment) +
  scale_color_manual(values=covercols) +
  scale_alpha_manual(values=c(0.3,1), labels=c("No glucose", "Glucose")) +
  labs(x = "Days after glucose addition", y = bquote("Daily CO" [2] * " Flux (" * mu * "g C g"^-1 * "dry soil)"),
       linetype = "Moisture treatment", alpha = "Glucose addition",
       color = "Residue type") +
  theme_bw()
p

ggsave(p, file="figures/emissions/CO2 flux.png",width = 8, height = 4, dpi = 300)


# Plot cumulative CO2 emissions from glucose averages
p <- ggplot(cumulative_avg, aes(x = Date, y = mean_CO2_glucose, color = Residue.type, alpha=Glucose.addition.treatment)) +
  geom_errorbar(aes(ymin = (mean_CO2_glucose - se_CO2_glucose), ymax = (mean_CO2_glucose + se_CO2_glucose)), width=2, size=0.5) +
  #geom_line(aes(linetype=Moisture.treatment), size=0.5) + 
  geom_line(size=0.5) + 
  facet_grid(.~Moisture.treatment) +
  scale_color_manual(values=covercols) +
  scale_alpha_manual(values=c(0.3,1)) +
  labs(x = "Days after glucose addition",  alpha = "Glucose addition",
       y = bquote("Glucose-derived C Emissions (" * mu * "g C g"^-1 * "dry soil)"),
       #y = bquote("Glucose-derived CO" [2] * "-C Emissions (%)"),
       linetype = "Moisture treatment",
       color = "Residue type") +
  theme_bw()
p

ggsave(p, file="figures/emissions/CO2 cumulative glucose.png",width = 8, height = 4, dpi = 300)

# Plot CO2 flux from glucose averages
p <- ggplot(cumulative_avg, aes(x = Date, y = mean_CO2_glucose_flux, color = Residue.type, alpha=Glucose.addition.treatment)) +
  geom_errorbar(aes(ymin = (mean_CO2_glucose_flux - se_CO2_glucose_flux), ymax = (mean_CO2_glucose_flux + se_CO2_glucose_flux)), width=2, size=0.5) +
  #geom_line(aes(linetype=Moisture.treatment), size=0.5) + 
  geom_line(size=0.5) + 
  facet_grid(.~Moisture.treatment) +
  scale_color_manual(values=covercols) +
  scale_alpha_manual(values=c(0.3,1), labels=c("No glucose", "Glucose")) +
  labs(x = "Days after glucose addition",  alpha = "Glucose addition",
       y = bquote("Daily Glucose-derived C Flux (" * mu * "g C g"^-1 * "dry soil)"),
       #y = bquote("Glucose-derived CO" [2] * "-C Emissions (%)"),
       linetype = "Moisture treatment",
       color = "Residue type") +
  theme_bw() 
p

ggsave(p, file="figures/emissions/CO2 flux glucose.png",width = 8, height = 4, dpi = 300)














 # set groupings
dat$Date <- as.Date(dat$Date, "%m/%d/%y")
dat$Residue.type <- as.factor(dat$Residue.type)
dat$Residue.type <- factor(dat$Residue.type, levels(dat$Residue.type)[c(2,5,1,4,3)])


# remove values prior to incubation initiation
dat <- dat[-which(dat$Date=="2023-05-22" & dat$Time=="1"),]


#calculate means
mean_dat <- dat %>% 
  group_by(Date, Time, Moisture.treatment) %>% 
  summarise(mean_WHC=mean(WHC),
            n_WHC=n())


#create mean line by moisture treatment
mplot <- ggplot(data=na.omit(mean_dat), 
       aes(x=Date, y=mean_WHC, group=Moisture.treatment, 
                                   fill=Moisture.treatment)) +
  theme_bw() +
  # geom_jitter(data = dat, width=0.25, height=0, alpha=0.25, colour="black", shape=21,
  #             mapping = aes(x = Date, y = WHC, fill = Moisture.treatment)) +
  geom_line(size=1, aes(color=Moisture.treatment), alpha=0.75)+
  #geom_point(shape=21, colour="black", size=3, alpha=1) + 
  labs(y="Water-holding capacity (%)", x="", color="Moisture\ntreatment") +
  scale_fill_manual(values=c("blue", "red")) +
  scale_color_manual(values=c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_date(date_breaks = "1 week",  date_labels = "%d %b %Y", 
               limit=c(as.Date("2023-05-22"),as.Date("2023-07-09"))) + 
  scale_y_continuous(sec.axis = sec_axis(~ . * 0.2482859/60, name = expression(paste("Gravimetric moisture (g g"^-1,")"))))
  
ggpubr::ggexport(mplot, height=1200, width=2000, filename="figures/0_moisture-treatments.png", res=400)


#calculate means
mean_dat3 <- dat %>% 
  group_by(Date, Time, Moisture.treatment) %>% 
  summarise(mean_WHC=mean(WHC, na.rm=T),
            n_WHC=n())

#create mean line by moisture treatment (revised)
tsize <- 2
mplot <- ggplot(data=(mean_dat3), 
                aes(x=Date, y=mean_WHC, group=Moisture.treatment, 
                    fill=Moisture.treatment)) +
  theme_bw() +
  # fill in the line
  geom_segment(x= as.Date("2023-06-19"), xend=as.Date("2023-06-21"), y=60, yend=58.5, color="blue", size=1, alpha=0.5) +
  geom_segment(x= as.Date("2023-06-21"), xend=as.Date("2023-06-23"), y=60, yend=58.8, color="blue", size=1, alpha=0.5) +
  geom_segment(x= as.Date("2023-06-29"), xend=as.Date("2023-06-30"), y=60, yend=58.8, color="blue", size=1, alpha=0.5) +
  geom_line(size=1, aes(color=Moisture.treatment), alpha=0.75)+
  labs(y="Water-holding capacity (%)", x="", color="Moisture\ntreatment") +
  scale_fill_manual(values=c("blue", "red")) +
  scale_color_manual(values=c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # x axis
    scale_x_date(date_breaks = "1 week",  date_labels = "%d %b %Y", 
               limit=c(as.Date("2023-05-20"),as.Date("2023-08-21"))) + 
  # y axis
  scale_y_continuous(sec.axis = sec_axis(~ . * 0.2482859/60, name = expression(paste("Gravimetric moisture (g g"^-1,")"))),
                     limits = c(15,95)) +
  # indicate residue addition
  geom_segment(x= as.Date("2023-05-22"), xend=as.Date("2023-05-22"), y=68, yend=62, arrow=arrow(length=unit(0.1, "cm"))) +
  annotate("text", label="Residue", x=as.Date("2023-05-22"), y=73, fontface="italic", size=tsize) +
  annotate("text", label="addition", x=as.Date("2023-05-22"), y=70, fontface="italic", size=tsize) +
  # indicate pre-incubation period
  geom_segment(x= as.Date("2023-05-22"), xend=as.Date("2023-05-29"), y=55, yend=55, linetype="dashed") +
  annotate("text", label="Pre-incubation", x=as.Date("2023-05-26"), y=50, fontface="italic", size=tsize) +
  annotate("text", label="(7 days)", x=as.Date("2023-05-26"), y=47, fontface="italic", size=tsize) +
  # indicate phase I
  geom_segment(x= as.Date("2023-05-29"), xend=as.Date("2023-06-30"), y=65, yend=65) +
  annotate("text", label="Phase I: Moisture manipulation", x=as.Date("2023-06-15"), y=73, fontface="italic", size=tsize) +
  annotate("text", label="(30 days)", x=as.Date("2023-06-15"), y=70, fontface="italic", size=tsize) +
  # indicate recovery period
  geom_segment(x= as.Date("2023-06-30"), xend=as.Date("2023-07-10"), y=55, yend=55, linetype="dashed") +
  annotate("text", label="Recovery", x=as.Date("2023-07-06"), y=50, fontface="italic", size=tsize) +
  annotate("text", label="(10 days)", x=as.Date("2023-07-06"), y=47, fontface="italic", size=tsize) +
  # indicate glucose addition to set 1
  geom_segment(x= as.Date("2023-07-10"), xend=as.Date("2023-07-10"), y=82, yend=68, arrow=arrow(length=unit(0.1, "cm"))) +
  annotate("text", label="Glucose addition (set 1)", x=as.Date("2023-07-10"), y=83, fontface="italic", size=tsize, hjust=0) +
  # indicate harvest set 1
  geom_segment(x= as.Date("2023-07-11"), xend=as.Date("2023-07-11"), y=78, yend=68, arrow=arrow(length=unit(0.1, "cm"))) +
  annotate("text", label="Harvest (set 1)", x=as.Date("2023-07-11"), y=79, fontface="italic", size=tsize, hjust=0) +
  # indicate glucose addition to set 2
  geom_segment(x= as.Date("2023-07-12"), xend=as.Date("2023-07-12"), y=74, yend=68, arrow=arrow(length=unit(0.1, "cm"))) +
  annotate("text", label="Glucose addition (set 2)", x=as.Date("2023-07-12"), y=75, fontface="italic", size=tsize, hjust=0) +
  # indicate phase II
  geom_segment(x= as.Date("2023-07-10"), xend=as.Date("2023-08-21"), y=86, yend=86, arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", label="Phase II: Glucose addition and tracking", x=as.Date("2023-07-10"), y=93, fontface="italic", size=tsize, hjust=0) +
  annotate("text", label="(six months)", x=as.Date("2023-07-10"), y=90, fontface="italic", size=tsize, hjust=0) 

  
mplot
ggpubr::ggexport(mplot, height=1600, width=2500, filename="figures/0_moisture-treatments_revised.png", res=400)


# 2023-05-22: initiate phase 1 (combine soil with residues and bring to 60% WHC)
# 2023-05-29: initiate drought
# 2023-06-30: bring all soils to 60% WHC
# 2023-07-10: glucose addition to set 1
# 2023-07-11: harvest set 1 (24 hr)
# 2023-07-12: glucose addition to set 2
# 2024-01-23: harvest set 2 (6 mo)



#calculate means
mean_dat2 <- dat %>% 
  group_by(Date, Time, Moisture.treatment, Residue.type) %>% 
  summarise(mean_WHC=mean(WHC),
            n_WHC=n())


#create mean line by moisture treatment and facet by residue type
mplot2 <- ggplot(data=na.omit(mean_dat2), 
                aes(x=Date, y=mean_WHC, group=Moisture.treatment, 
                    fill=Moisture.treatment)) +
  facet_wrap(.~Residue.type, nrow=3)+
  theme_bw() +
  # geom_jitter(data = dat, width=0.25, height=0, alpha=0.25, colour="black", shape=21,
  #             mapping = aes(x = Date, y = WHC, fill = Moisture.treatment)) +
  geom_line(size=2, aes(color=Moisture.treatment), alpha=0.75)+
  #geom_point(shape=21, colour="black", size=3, alpha=1) + 
  labs(y="Water-holding capacity (%)", x="", color="Moisture\ntreatment") +
  scale_fill_manual(values=c("blue", "red")) +
  scale_color_manual(values=c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust=1)) +
  scale_x_date(date_breaks = "1 week",  date_labels = "%d %b %Y", 
               limit=c(as.Date("2023-05-22"),as.Date("2023-07-09"))) + 
  scale_y_continuous(sec.axis = sec_axis(~ . * 0.2482859/60, name = expression(paste("Gravimetric moisture (g g"^-1,")"))))

ggpubr::ggexport(mplot2, height=1200*3, width=2000*2, filename="figures/moisture-treatments2.png", res=400)
