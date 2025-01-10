

library(dplyr)
library(ggplot2)

# 2023-05-22: initiate phase 1 (combine soil with residues and bring to 60% WHC)
# 2023-05-29: initiate drought
# 2023-06-30: bring all soils to 60% WHC
# 2023-07-10: glucose addition to set 1
# 2023-07-11: harvest set 1 (24 hr)
# 2023-07-12: glucose addition to set 2
# 2024-01-23: harvest set 2 (6 mo)



dat <- read.csv("raw-data/master-dataset-moisture.csv")
# set groupings
dat$Date <- as.Date(dat$Date, "%m/%d/%y")
dat$Residue.type <- as.factor(dat$Residue.type)
dat$Residue.type <- factor(dat$Residue.type, levels(dat$Residue.type)[c(2,5,1,4,3)])


# remove values prior to incubation initiation
dat <- dat[-which(dat$Date=="2023-05-22" & dat$Time=="1"),]


#calculate means
mean_dat <- dat %>% 
  group_by(Date, Time, Moisture.treatment) %>% 
  summarise(mean_WHC=mean(WHC))


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
mean_dat2 <- dat %>% 
  group_by(Date, Time, Moisture.treatment, Residue.type) %>% 
  summarise(mean_WHC=mean(WHC))


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
