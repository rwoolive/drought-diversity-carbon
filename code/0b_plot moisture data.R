

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
