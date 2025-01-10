## R version 4.3.3 (2024-02-29) -- "Angel Food Cake"


library(ggplot2)
library(hms)
library(ggpubr)

setwd("/Users/rachelwooliver/Library/CloudStorage/GoogleDrive-rwoolive@utk.edu/My Drive/*1_UTK-postdoc_2020-/6_Incubation_cover crop diversity")


#########################
# MAOM 6 month harvest
#########################



# metadata for unlabeled maom fractions from 6 month harvest
labs <- read.csv("entered data/raw data/2024_03-07-rachel incubation unlabeled maom.csv")
path <- c("entered data/raw data/Picarro real time log/2024/03/07", "entered data/raw data/Picarro real time log/2024/03/08")
filenames <- list.files(path, pattern="*.dat", full.names=TRUE)
ldf <- lapply(filenames, read.table, header=TRUE)
dat <- ldf[[1]]
for(i in 2:37){dat <- rbind(dat, ldf[[i]])}
dat$TIME2 <- as.hms(dat$TIME)
dat$DATE2 <- as.POSIXct(dat$DATE)
dat$DATETIME <- as.POSIXct(paste(dat$DATE2, dat$TIME2), format="%Y-%m-%d %H:%M:%S") - 4*60*60 # subtract 4 hours due to computer error
p <- ggplot(dat[40000:90000,], aes(DATETIME, Delta_Raw)) +
  geom_point() + lims(y=c(-60,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
p
q <- ggplot(dat[40000:90000,], aes(DATETIME, X12CO2)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
q
ggexport(ggarrange(p,q, nrow=2), height=1500, width=10000, filename="usda-diversity-incubation/figures/hand-integration/0a_01picarro-maom-unlabeled-6mo.png")



# metadata for labeled maom fractions from 6 month harvest
labs <- read.csv("entered data/raw data/2024_03-12-rachel incubation labeled maom.csv")
path <- c("entered data/raw data/Picarro real time log/2024/03/12", "entered data/raw data/Picarro real time log/2024/03/13")
filenames <- list.files(path, pattern="*.dat", full.names=TRUE)
ldf <- lapply(filenames, read.table, header=TRUE)
dat <- ldf[[1]]
for(i in 2:37){dat <- rbind(dat, ldf[[i]])}
dat$TIME2 <- as.hms(dat$TIME)
dat$DATE2 <- as.POSIXct(dat$DATE)
dat$DATETIME <- as.POSIXct(paste(dat$DATE2, dat$TIME2), format="%Y-%m-%d %H:%M:%S") - 4*60*60 # subtract 4 hours due to computer error
p <- ggplot(dat[35000:85000,], aes(DATETIME, Delta_Raw)) +
  geom_point() + lims(y=c(-60,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
p
q <- ggplot(dat[35000:85000,], aes(DATETIME, X12CO2)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
q
ggexport(ggarrange(p,q, nrow=2), height=1500, width=10000, filename="usda-diversity-incubation/figures/hand-integration/0a_02picarro-maom-labeled-6mo.png")




#########################
# POM 6 month harvest
#########################



# metadata for unlabeled pom fractions from 6 month harvest
labs <- read.csv("entered data/raw data/2024_03-06-rachel incubation unlabeled pom.csv")
path <- c("entered data/raw data/Picarro real time log/2024/03/06", "entered data/raw data/Picarro real time log/2024/03/07")
filenames <- list.files(path, pattern="*.dat", full.names=TRUE)
ldf <- lapply(filenames, read.table, header=TRUE)
dat <- ldf[[1]]
for(i in 2:37){dat <- rbind(dat, ldf[[i]])}
dat$TIME2 <- as.hms(dat$TIME)
dat$DATE2 <- as.POSIXct(dat$DATE)
dat$DATETIME <- as.POSIXct(paste(dat$DATE2, dat$TIME2), format="%Y-%m-%d %H:%M:%S") - 4*60*60 # subtract 4 hours due to computer error
p <- ggplot(dat[40000:90000,], aes(DATETIME, Delta_Raw)) +
  geom_point() + lims(y=c(-60,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
p
q <- ggplot(dat[40000:90000,], aes(DATETIME, X12CO2)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
q
ggexport(ggarrange(p,q, nrow=2), height=1500, width=10000, filename="usda-diversity-incubation/figures/hand-integration/0a_03picarro-pom-unlabeled-6mo.png")



# metadata for labeled pom fractions from 6 month harvest
labs <- read.csv("entered data/raw data/2024_03-11-rachel incubation labeled pom.csv")
path <- c("entered data/raw data/Picarro real time log/2024/03/11", "entered data/raw data/Picarro real time log/2024/03/12")
filenames <- list.files(path, pattern="*.dat", full.names=TRUE)
ldf <- lapply(filenames, read.table, header=TRUE)
dat <- ldf[[1]]
for(i in 2:37){dat <- rbind(dat, ldf[[i]])}
dat$TIME2 <- as.hms(dat$TIME)
dat$DATE2 <- as.POSIXct(dat$DATE)
dat$DATETIME <- as.POSIXct(paste(dat$DATE2, dat$TIME2), format="%Y-%m-%d %H:%M:%S") - 4*60*60 # subtract 4 hours due to computer error
p <- ggplot(dat[35000:90000,], aes(DATETIME, Delta_Raw)) +
  geom_point() + lims(y=c(-60,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
p
q <- ggplot(dat[35000:90000,], aes(DATETIME, X12CO2)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
q
ggexport(ggarrange(p,q, nrow=2), height=1500, width=10000, filename="usda-diversity-incubation/figures/hand-integration/0a_04picarro-pom-labeled-6mo.png")








#########################
# SOM 6 month harvest
#########################



# metadata for unlabeled som fractions from 6 month harvest
labs <- read.csv("entered data/raw data/2024_03-04-rachel incubation unlabeled soil.csv")
path <- c("entered data/raw data/Picarro real time log/2024/03/04", "entered data/raw data/Picarro real time log/2024/03/05")
filenames <- list.files(path, pattern="*.dat", full.names=TRUE)
ldf <- lapply(filenames, read.table, header=TRUE)
dat <- ldf[[1]]
for(i in 2:37){dat <- rbind(dat, ldf[[i]])}
dat$TIME2 <- as.hms(dat$TIME)
dat$DATE2 <- as.POSIXct(dat$DATE)
dat$DATETIME <- as.POSIXct(paste(dat$DATE2, dat$TIME2), format="%Y-%m-%d %H:%M:%S") - 4*60*60 # subtract 4 hours due to computer error
p <- ggplot(dat[40000:80000,], aes(DATETIME, Delta_Raw)) +
  geom_point() + lims(y=c(-60,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
p
q <- ggplot(dat[40000:80000,], aes(DATETIME, X12CO2)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
q
ggexport(ggarrange(p,q, nrow=2), height=1500, width=10000, filename="usda-diversity-incubation/figures/hand-integration/0a_05picarro-som-unlabeled-6mo.png")



# metadata for labeled som fractions from 6 month harvest
labs <- read.csv("entered data/raw data/2024_03-05-rachel incubation labeled soil.csv")
path <- c("entered data/raw data/Picarro real time log/2024/03/05", "entered data/raw data/Picarro real time log/2024/03/06")
filenames <- list.files(path, pattern="*.dat", full.names=TRUE)
ldf <- lapply(filenames, read.table, header=TRUE)
dat <- ldf[[1]]
for(i in 2:37){dat <- rbind(dat, ldf[[i]])}
dat$TIME2 <- as.hms(dat$TIME)
dat$DATE2 <- as.POSIXct(dat$DATE)
dat$DATETIME <- as.POSIXct(paste(dat$DATE2, dat$TIME2), format="%Y-%m-%d %H:%M:%S") - 4*60*60 # subtract 4 hours due to computer error
p <- ggplot(dat[40000:100000,], aes(DATETIME, Delta_Raw)) +
  geom_point() + lims(y=c(-60,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
p
q <- ggplot(dat[40000:100000,], aes(DATETIME, X12CO2)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
q
ggexport(ggarrange(p,q, nrow=2), height=1500, width=10000, filename="usda-diversity-incubation/figures/hand-integration/0a_06picarro-som-labeled-6mo.png")









#########################
# MAOM 24 hour harvest
#########################



# metadata for unlabeled maom fractions 
labs <- read.csv("entered data/raw data/2023_11_15- avi microplastic incubation soils 155-160, rachel incubation unlabeled soils.csv")
path <- c("entered data/raw data/Picarro real time log/2023/11/15", "entered data/raw data/Picarro real time log/2023/11/16")
filenames <- list.files(path, pattern="*.dat", full.names=TRUE)
ldf <- lapply(filenames, read.table, header=TRUE)
dat <- ldf[[1]]
for(i in 2:37){dat <- rbind(dat, ldf[[i]])}
dat$TIME2 <- as.hms(dat$TIME)
dat$DATE2 <- as.POSIXct(dat$DATE)
dat$DATETIME <- as.POSIXct(paste(dat$DATE2, dat$TIME2), format="%Y-%m-%d %H:%M:%S") - 4*60*60 # subtract 4 hours due to computer error
p <- ggplot(dat[55000:109000,], aes(DATETIME, Delta_Raw)) +
  geom_point() + lims(y=c(-60,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
p
q <- ggplot(dat[55000:109000,], aes(DATETIME, X12CO2)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
q
ggexport(ggarrange(p,q, nrow=2), height=1500, width=10000, filename="usda-diversity-incubation/figures/hand-integration/0a_07picarro-maom-unlabeled-24hr.png")



# metadata for labeled maom fractions
labs <- read.csv("entered data/raw data/2023_11_16-rachel incubation labeled soils.csv")
path <- c("entered data/raw data/Picarro real time log/2023/11/16", "entered data/raw data/Picarro real time log/2023/11/17")
filenames <- list.files(path, pattern="*.dat", full.names=TRUE)
ldf <- lapply(filenames, read.table, header=TRUE)
dat <- ldf[[1]]
for(i in 2:37){dat <- rbind(dat, ldf[[i]])}
dat$TIME2 <- as.hms(dat$TIME)
dat$DATE2 <- as.POSIXct(dat$DATE)
dat$DATETIME <- as.POSIXct(paste(dat$DATE2, dat$TIME2), format="%Y-%m-%d %H:%M:%S") - 4*60*60 # subtract 4 hours due to computer error
p <- ggplot(dat[30000:100000,], aes(DATETIME, Delta_Raw)) +
  geom_point() + lims(y=c(-60,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
p
q <- ggplot(dat[30000:100000,], aes(DATETIME, X12CO2)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
q
ggexport(ggarrange(p,q, nrow=2), height=1500, width=10000, filename="usda-diversity-incubation/figures/hand-integration/0a_08picarro-maom-labeled-24hr.png")




#########################
# POM 24 hour harvest
#########################



# metadata for unlabeled pom fractions 
labs <- read.csv("entered data/raw data/2024_02_01_rachel incubation unlabeled pom and sand fraction.csv")
path <- c("entered data/raw data/Picarro real time log/2024/02/01", "entered data/raw data/Picarro real time log/2024/02/02")
filenames <- list.files(path, pattern="*.dat", full.names=TRUE)
ldf <- lapply(filenames, read.table, header=TRUE)
dat <- ldf[[1]]
for(i in 2:37){dat <- rbind(dat, ldf[[i]])}
dat$TIME2 <- as.hms(dat$TIME)
dat$DATE2 <- as.POSIXct(dat$DATE)
dat$DATETIME <- as.POSIXct(paste(dat$DATE2, dat$TIME2), format="%Y-%m-%d %H:%M:%S") - 4*60*60 # subtract 4 hours due to computer error
p <- ggplot(dat[55000:107000,], aes(DATETIME, Delta_Raw)) +
  geom_point() + lims(y=c(-60,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
p
q <- ggplot(dat[55000:107000,], aes(DATETIME, X12CO2)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
q
ggexport(ggarrange(p,q, nrow=2), height=1500, width=10000, filename="usda-diversity-incubation/figures/hand-integration/0a_09picarro-pom-unlabeled-24hr.png")



# metadata for labeled pom fractions 
labs <- read.csv("entered data/raw data/2024_02_02_rachel incubation labeled pom and sand fraction and binsiya soils.csv")
path <- c("entered data/raw data/Picarro real time log/2024/02/02", "entered data/raw data/Picarro real time log/2024/02/03")
filenames <- list.files(path, pattern="*.dat", full.names=TRUE)
ldf <- lapply(filenames, read.table, header=TRUE)
dat <- ldf[[1]]
for(i in 2:37){dat <- rbind(dat, ldf[[i]])}
dat$TIME2 <- as.hms(dat$TIME)
dat$DATE2 <- as.POSIXct(dat$DATE)
dat$DATETIME <- as.POSIXct(paste(dat$DATE2, dat$TIME2), format="%Y-%m-%d %H:%M:%S") - 4*60*60 # subtract 4 hours due to computer error
p <- ggplot(dat[25000:108000,], aes(DATETIME, Delta_Raw)) +
  geom_point() + lims(y=c(-60,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
p
q <- ggplot(dat[25000:108000,], aes(DATETIME, X12CO2)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
q
ggexport(ggarrange(p,q, nrow=2), height=1500, width=10000, filename="usda-diversity-incubation/figures/hand-integration/0a_10picarro-pom-labeled-24hr.png")

path <- c("entered data/raw data/Picarro real time log/2024/02/03", "entered data/raw data/Picarro real time log/2024/02/04")
filenames <- list.files(path, pattern="*.dat", full.names=TRUE)
ldf <- lapply(filenames, read.table, header=TRUE)
dat <- ldf[[1]]
for(i in 2:37){dat <- rbind(dat, ldf[[i]])}
dat$TIME2 <- as.hms(dat$TIME)
dat$DATE2 <- as.POSIXct(dat$DATE)
dat$DATETIME <- as.POSIXct(paste(dat$DATE2, dat$TIME2), format="%Y-%m-%d %H:%M:%S") - 4*60*60 # subtract 4 hours due to computer error
p <- ggplot(dat, aes(DATETIME, Delta_Raw)) +
  geom_point() + lims(y=c(-60,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
p
q <- ggplot(dat, aes(DATETIME, X12CO2)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
q
ggexport(ggarrange(p,q, nrow=2), height=1500, width=10000, filename="usda-diversity-incubation/figures/hand-integration/0a_11picarro-pom-labeled-24hr part 2.png")








#########################
# SOM 24 hr harvest
#########################



# metadata for unlabeled som fractions 
labs <- read.csv("entered data/raw data/2023_08_17_RACHEL 1-44 UNLABELED RP TEST 10 PICARRO.csv")
path <- c("entered data/raw data/Picarro real time log/2023/08/16", "entered data/raw data/Picarro real time log/2023/08/17")
filenames <- list.files(path, pattern="*.dat", full.names=TRUE)
ldf <- lapply(filenames, read.table, header=TRUE)
dat <- ldf[[1]]
for(i in 2:37){dat <- rbind(dat, ldf[[i]])}
dat$TIME2 <- as.hms(dat$TIME)
dat$DATE2 <- as.POSIXct(dat$DATE)
dat$DATETIME <- as.POSIXct(paste(dat$DATE2, dat$TIME2), format="%Y-%m-%d %H:%M:%S") - 4*60*60 # subtract 4 hours due to computer error
p <- ggplot(dat[50000:100000,], aes(DATETIME, Delta_Raw)) +
  geom_point() + lims(y=c(-60,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
p
q <- ggplot(dat[50000:100000,], aes(DATETIME, X12CO2)) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
q
ggexport(ggarrange(p,q, nrow=2), height=1500, width=10000, filename="usda-diversity-incubation/figures/hand-integration/0a_12picarro-som-unlabeled-24hr.png")



# metadata for labeled som fractions
labs <- read.csv("entered data/raw data/2023_08_17_RACHEL 1-42 LABELED 10 PICARRO CFGDS2373_IsotopicData_20230818_210501.csv")
path <- c("entered data/raw data/Picarro real time log/2023/08/18", "entered data/raw data/Picarro real time log/2023/08/19")
filenames <- list.files(path, pattern="*.dat", full.names=TRUE)
ldf <- lapply(filenames, read.table, header=TRUE)
dat <- ldf[[1]]
for(i in 2:37){dat <- rbind(dat, ldf[[i]])}
dat$TIME2 <- as.hms(dat$TIME)
dat$DATE2 <- as.POSIXct(dat$DATE)
dat$DATETIME <- as.POSIXct(paste(dat$DATE2, dat$TIME2), format="%Y-%m-%d %H:%M:%S") - 4*60*60 # subtract 4 hours due to computer error
p <- ggplot(dat[c(60000:100000),], aes(DATETIME, Delta_Raw)) +
  geom_point() + lims(y=c(-60,0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p
q <- ggplot(dat[c(60000:100000),], aes(DATETIME, X12CO2)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
q
ggexport(ggarrange(p,q, nrow=2), height=1500, width=10000, filename="usda-diversity-incubation/figures/hand-integration/0a_13picarro-som-labeled-24hr.png")


