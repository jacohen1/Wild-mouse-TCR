####initial code####
rm(list=ls())

#set working directory
#setwd("O:/Documents/Empirical/Wild_mouse_tcr/Data/tcrs_wild_mice/")
setwd("C:/Users/JCohen94/Desktop/Wild_mouse_tcr/Data/tcrs_wild_mice/")

##load necessary packages
library(immunarch)
library(stringr)
library(igraph)
library(ggplot2)
library(showtext)
font_add(family = "Calibri", regular = "C:/Windows/Fonts/calibri.ttf")
showtext_auto()
showtext_opts(dpi = 300)
library(dplyr)
library(reshape2)
library(network)
library(sna)
library(GGally)

####Read in data####
##read in mouse metadata file
mouse.data <- read.csv("../Mouse_metadata_tcr_mice_Jacob_2_rem.csv", header = TRUE)

##data to plot diversity
simpsons <- as.data.frame(rbind(cbind("cd4a", mouse.data$CD4.alpha.simpsons, mouse.data$Age_days, mouse.data$ID, mouse.data$Site),
                                cbind("cd4b", mouse.data$CD4.beta.simpsons, mouse.data$Age_days, mouse.data$ID, mouse.data$Site),
                                cbind("cd8a", mouse.data$CD8.alpha.simpsons, mouse.data$Age_days, mouse.data$ID, mouse.data$Site),
                                cbind("cd8b", mouse.data$CD8.beta.simpsons, mouse.data$Age_days, mouse.data$ID, mouse.data$Site)))

shannons <- as.data.frame(rbind(cbind("cd4a", mouse.data$CD4.alpha.shannons, mouse.data$Age_days, mouse.data$ID, mouse.data$Site),
                                cbind("cd4b", mouse.data$CD4.beta.shannons, mouse.data$Age_days, mouse.data$ID, mouse.data$Site),
                                cbind("cd8a", mouse.data$CD8.alpha.shannons, mouse.data$Age_days, mouse.data$ID, mouse.data$Site),
                                cbind("cd8b", mouse.data$CD8.beta.shannons, mouse.data$Age_days, mouse.data$ID, mouse.data$Site)))

colnames(simpsons) <- c("type", "value", "age", "id", "site")
colnames(shannons) <- c("type", "value", "age", "id", "site")

##data for unique aa sequences
unique.data <- as.data.frame(rbind(cbind("cd4a", mouse.data$CD4.alpha.unique.seq, mouse.data$Age_days, mouse.data$ID, mouse.data$Site),
                                   cbind("cd4b", mouse.data$CD4.beta.unique.seq, mouse.data$Age_days, mouse.data$ID, mouse.data$Site),
                                   cbind("cd8a", mouse.data$CD8.alpha.unique.seq, mouse.data$Age_days, mouse.data$ID, mouse.data$Site),
                                   cbind("cd8b", mouse.data$CD8.beta.unique.seq, mouse.data$Age_days, mouse.data$ID, mouse.data$Site)))

colnames(unique.data) <- c("type", "value", "age", "id", "site")

#read in repertoire data
rep.cd4a <- read.csv("../CD4_alpha_Repertoire.csv", header = TRUE)
rep.cd4b <- read.csv("../CD4_beta_Repertoire.csv", header = TRUE)
rep.cd8a <- read.csv("../CD8_alpha_Repertoire.csv", header = TRUE)
rep.cd8b <- read.csv("../CD8_beta_Repertoire.csv", header = TRUE)

#read in TCR frequency data
tcr.cd4a <- read.csv("../CD4_alpha_TCR_freq.csv", header = TRUE)
tcr.cd4b <- read.csv("../CD4_beta_TCR_freq.csv", header = TRUE)
tcr.cd8a <- read.csv("../CD8_alpha_TCR_freq.csv", header = TRUE)
tcr.cd8b <- read.csv("../CD8_beta_TCR_freq.csv", header = TRUE)

####Clean up data####
#remove unneeded mice from repertoire data
rep.cd4a <- subset(rep.cd4a, select = -c(N63, N42, N24, N01, N27, N11, N12, N09, N35),
                   !(rep.cd4a$X %in% c("N63", "N42", "N24", "N01", "N27", "N11", "N12", "N09", "N35")))
rep.cd4b <- subset(rep.cd4b, select = -c(N63, N42, N24, N01, N27, N11, N12, N09, N35),
                   !(rep.cd4b$X %in% c("N63", "N42", "N24", "N01", "N27", "N11", "N12", "N09", "N35")))
rep.cd8a <- subset(rep.cd8a, select = -c(N17, N42, N22, N01, N21, N11, N09, N27, N12, N51, N10, N35, N36),
                   !(rep.cd8a$X %in% c("N17", "N42", "N22", "N01", "N21", "N11", "N09", "N27", "N12", "N51", "N10", "N35", "N36")))
rep.cd8b <- subset(rep.cd8b, select = -c(N17, N42, N22, N01, N21, N11, N09, N27, N12, N51, N10, N35, N36),
                   !(rep.cd8b$X %in% c("N17", "N42", "N22", "N01", "N21", "N11", "N09", "N27", "N12", "N51", "N10", "N35", "N36")))

#change from matrix format to long format
cd4a.net <- as.data.frame.table(`row.names<-`(as.matrix(rep.cd4a)[,-1],
                                                   rep.cd4a$X))
cd4b.net <- as.data.frame.table(`row.names<-`(as.matrix(rep.cd4b)[,-1],
                                                   rep.cd4b$X))
cd8a.net <- as.data.frame.table(`row.names<-`(as.matrix(rep.cd8a)[,-1],
                                                   rep.cd8a$X))
cd8b.net <- as.data.frame.table(`row.names<-`(as.matrix(rep.cd8b)[,-1],
                                                   rep.cd8b$X))

##get mouse IDs
cd4.names <- rep.cd4a$X
cd8.names <- rep.cd8a$X

##drop mouse IDs from the heatmap data
rep.cd4a$X <- 1:62
colnames(rep.cd4a) <- c("X", 1:62)
rep.cd4b$X <- 1:62
colnames(rep.cd4b) <- c("X", 1:62)
rep.cd8a$X <- 1:58
colnames(rep.cd8a) <- c("X", 1:58)
rep.cd8b$X <- 1:58
colnames(rep.cd8b) <- c("X", 1:58)

#change from matrix format to long format
rep.cd4a.long <- as.data.frame.table(`row.names<-`(as.matrix(rep.cd4a)[,-1],
                                                   rep.cd4a$X))
rep.cd4b.long <- as.data.frame.table(`row.names<-`(as.matrix(rep.cd4b)[,-1],
                                                   rep.cd4b$X))
rep.cd8a.long <- as.data.frame.table(`row.names<-`(as.matrix(rep.cd8a)[,-1],
                                                   rep.cd8a$X))
rep.cd8b.long <- as.data.frame.table(`row.names<-`(as.matrix(rep.cd8b)[,-1],
                                                   rep.cd8b$X))

####Plot repertoire size by age####
cd4a.size.plot <- ggplot(data = mouse.data, aes(x = Age_days, y = CD4.alpha.count, color = Site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "CD4 alpha repertoire size") +
  scale_color_manual(name = "Site", labels = c("Lab", "Nott.", "Wirral"),
                     values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

cd4b.size.plot <- ggplot(data = mouse.data, aes(x = Age_days, y = CD4.beta.count, color = Site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "CD4 beta repertoire size") +
  scale_color_manual(name = "Site", labels = c("Lab", "Nott.", "Wirral"),
                     values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

cd8a.size.plot <- ggplot(data = mouse.data, aes(x = Age_days, y = CD8.alpha.count, color = Site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "CD8 alpha repertoire size") +
  scale_color_manual(name = "Site", labels = c("Lab", "Nott.", "Wirral"),
                     values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

cd8b.size.plot <- ggplot(data = mouse.data, aes(x = Age_days, y = CD8.beta.count, color = Site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "CD8 beta repertoire size") +
  scale_color_manual(name = "Site", labels = c("Lab", "Nott.", "Wirral"),
                     values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

ggsave("../../Results/Exploratory_plots/CD4_alpha_age_repertoire_size.png", cd4a.size.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/CD4_beta_age_repertoire_size.png", cd4b.size.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/CD8_alpha_age_repertoire_size.png", cd8a.size.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/CD8_beta_age_repertoire_size.png", cd8b.size.plot, bg = "transparent", units = "cm", width = 30, height = 20)

####Plot TCR count results####
##CD4 results
cd4.tcr.plot <- ggplot(data = mouse.data, aes(x = CD4.alpha.count, y = CD4.beta.count, color = Site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "CD4 alpha repertoire size") +
  scale_y_continuous(name = "CD4 beta repertoire size") +
  scale_color_manual(name = "Site", labels = c("Lab", "Nott.", "Wirral"),
                     values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(0, 304000), xlim = c(0, 304000))

cd4.tcr.plot.log <- ggplot(data = mouse.data, aes(x = log(CD4.alpha.count), y = log(CD4.beta.count), color = Site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Log, CD4 alpha repertoire size") +
  scale_y_continuous(name = "Log, CD4 beta repertoire size") +
  scale_color_manual(name = "Site", labels = c("Lab", "Nott.", "Wirral"),
                     values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(0, 14), xlim = c(0, 14))

##CD8 results
cd8.tcr.plot <- ggplot(data = mouse.data, aes(x = CD8.alpha.count, y = CD8.beta.count, color = Site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "CD8 alpha repertoire size", breaks = c(0, 100000, 200000)) +
  scale_y_continuous(name = "CD8 beta repertoire size", breaks = c(0, 100000, 200000)) +
  scale_color_manual(name = "Site", labels = c("Lab", "Nott.", "Wirral"),
                     values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(0, 238000), xlim = c(0, 238000))

cd8.tcr.plot.log <- ggplot(data = mouse.data, aes(x = log(CD8.alpha.count), y = log(CD8.beta.count), color = Site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Log, CD8 repertoire size") +
  scale_y_continuous(name = "Log, CD8 repertoire size") +
  scale_color_manual(name = "Site", labels = c("Lab", "Nott.", "Wirral"),
                     values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(4, 13), xlim = c(4, 13))

##Save plots
ggsave("../../Results/Exploratory_plots/CD4_TCR_count_rem.png", cd4.tcr.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/CD4_TCR_count_rem_log.png", cd4.tcr.plot.log, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/CD8_TCR_count_rem.png", cd8.tcr.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/CD8_TCR_count_rem_log.png", cd8.tcr.plot.log, bg = "transparent", units = "cm", width = 30, height = 20)

####Plot alpha coefficient against mouse age####
CD4.alpha.powerlaw.plot <- ggplot(data = mouse.data, aes(x = Age_days, y = CD4.alpha.powerlaw, color = Site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "CD4 alpha power law coefficient") +
  scale_color_manual(name = "Site", labels = c("Lab", "Nott.", "Wirral"),
                     values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(0, 6), xlim = c(0, 350))

CD4.beta.powerlaw.plot <- ggplot(data = mouse.data, aes(x = Age_days, y = CD4.beta.powerlaw, color = Site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "CD4 beta power law coefficient") +
  scale_color_manual(name = "Site", labels = c("Lab", "Nott.", "Wirral"),
                     values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(0, 6), xlim = c(0, 350))

CD8.alpha.powerlaw.plot <- ggplot(data = mouse.data, aes(x = Age_days, y = CD8.alpha.powerlaw, color = Site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "CD8 alpha power law coefficient") +
  scale_color_manual(name = "Site", labels = c("Lab", "Nott.", "Wirral"),
                     values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(0, 6), xlim = c(0, 350))

CD8.beta.powerlaw.plot <- ggplot(data = mouse.data, aes(x = Age_days, y = CD8.beta.powerlaw, color = Site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "CD8 beta power law coefficient") +
  scale_color_manual(name = "Site", labels = c("Lab", "Nott.", "Wirral"),
                     values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(0, 6), xlim = c(0, 350))

#Save plots
ggsave("../../Results/Exploratory_plots/CD4_alpha_power_law_coefficient.png", CD4.alpha.powerlaw.plot, bg = "transparent", units = "cm", width = 30, height = 23)
ggsave("../../Results/Exploratory_plots/CD4_beta_power_law_coefficient.png", CD4.beta.powerlaw.plot, bg = "transparent", units = "cm", width = 30, height = 23)
ggsave("../../Results/Exploratory_plots/CD8_alpha_power_law_coefficient.png", CD8.alpha.powerlaw.plot, bg = "transparent", units = "cm", width = 30, height = 23)
ggsave("../../Results/Exploratory_plots/CD8_beta_power_law_coefficient.png", CD8.beta.powerlaw.plot, bg = "transparent", units = "cm", width = 30, height = 23)

####Plot no. unique aa sequences####
unique.plot <- ggplot(unique.data, aes(x = type, y = as.numeric(value), color = type)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("CD4 alpha", "CD4 beta", "CD8 alpha", "CD8 beta")) +
  scale_y_continuous(name = "No. unique AA sequences") +
  ggtitle("") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

#separate by cd4a/b and cd8a/b
unique.cd4a.plot <- ggplot(unique.data[which(unique.data$type == "cd4a"), ], aes(x = site, y = as.numeric(value), color = site)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("Lab", "Nott.", "Wirral")) +
  scale_y_continuous(name = "No. unique AA sequences") +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD4 alpha") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(5000, 67500))

unique.cd4b.plot <- ggplot(unique.data[which(unique.data$type == "cd4b"), ], aes(x = site, y = as.numeric(value), color = site)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("Lab", "Nott.", "Wirral")) +
  scale_y_continuous(name = "No. unique AA sequences") +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD4 beta") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(5000, 67500))

unique.cd8a.plot <- ggplot(unique.data[which(unique.data$type == "cd8a"), ], aes(x = site, y = as.numeric(value), color = site)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("Lab", "Nott.", "Wirral")) +
  scale_y_continuous(name = "No. unique AA sequences") +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD8 alpha") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(5000, 67500))

unique.cd8b.plot <- ggplot(unique.data[which(unique.data$type == "cd8b"), ], aes(x = site, y = as.numeric(value), color = site)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("Lab", "Nott.", "Wirral")) +
  scale_y_continuous(name = "No. unique AA sequences") +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD8 beta") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(5000, 67500))

##by age
unique.age.plot <- ggplot(unique.data, aes(x = as.numeric(age), y = as.numeric(value), color = type, shape = site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "No. unique AA sequences") +
  ggtitle("") +
  scale_color_discrete(labels = c("CD4 alpha", "CD4 beta", "CD8 alpha", "CD8 beta"), name = "") +
  scale_shape_discrete(labels = c("Lab", "Nott.", "Wirral"), name = "") +
  guides(color = guide_legend(override.aes = list(size = 6)), shape = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

#separate by cd4a/b and cd8a/b
unique.age.cd4a.plot <- ggplot(unique.data[which(unique.data$type == "cd4a"), ], aes(x = as.numeric(age), y = as.numeric(value), color = site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "No. unique AA sequences") +
  ggtitle("CD4 alpha") +
  scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(0, 67500), xlim = c(0, 315))

unique.age.cd4b.plot <- ggplot(unique.data[which(unique.data$type == "cd4b"), ], aes(x = as.numeric(age), y = as.numeric(value), color = site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "No. unique AA sequences") +
  ggtitle("CD4 beta") +
  scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(0, 67500), xlim = c(0, 315))

unique.age.cd8a.plot <- ggplot(unique.data[which(unique.data$type == "cd8a"), ], aes(x = as.numeric(age), y = as.numeric(value), color = site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "No. unique AA sequences") +
  ggtitle("CD8 alpha") +
  scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(0, 67500), xlim = c(0, 315))

unique.age.cd8b.plot <- ggplot(unique.data[which(unique.data$type == "cd8b"), ], aes(x = as.numeric(age), y = as.numeric(value), color = site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "No. unique AA sequences") +
  ggtitle("CD8 beta") +
  scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(0, 67500), xlim = c(0, 315))

ggsave("../../Results/Exploratory_plots/Unique_aa.png", unique.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Unique_aa_cd4a.png", unique.cd4a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Unique_aa_cd4b.png", unique.cd4b.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Unique_aa_cd8a.png", unique.cd8a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Unique_aa_cd8b.png", unique.cd8b.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Unique_aa_age.png", unique.age.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Unique_aa_age_cd4a.png", unique.age.cd4a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Unique_aa_age_cd4b.png", unique.age.cd4b.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Unique_aa_age_cd8a.png", unique.age.cd8a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Unique_aa_age_cd8b.png", unique.age.cd8b.plot, bg = "transparent", units = "cm", width = 30, height = 20)

####Plot TCR frequencies####
##CD4 alpha
tcr.cd4a.hist <- ggplot(tcr.cd4a, aes(x = count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "white") +
  scale_x_continuous(name = "No. total repeats") +
  scale_y_continuous(name = "Frequency") +
  ggtitle("CD4 alpha") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

tcr.cd4a.hist2 <- ggplot(tcr.cd4a, aes(x = count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "white") +
  scale_x_continuous(name = "No. total repeats") +
  scale_y_log10(name = "Log10, Frequency") +
  ggtitle("CD4 alpha") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

##CD4 beta
tcr.cd4b.hist <- ggplot(tcr.cd4b, aes(x = count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "white") +
  scale_x_continuous(name = "No. total repeats") +
  scale_y_continuous(name = "Frequency") +
  ggtitle("CD4 beta") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

tcr.cd4b.hist2 <- ggplot(tcr.cd4b, aes(x = count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "white") +
  scale_x_continuous(name = "No. total repeats") +
  scale_y_log10(name = "Log10, Frequency") +
  ggtitle("CD4 beta") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

##CD8 alpha
tcr.cd8a.hist <- ggplot(tcr.cd8a, aes(x = count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "white") +
  scale_x_continuous(name = "No. total repeats") +
  scale_y_continuous(name = "Frequency") +
  ggtitle("CD8 alpha") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

tcr.cd8a.hist2 <- ggplot(tcr.cd8a, aes(x = count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "white") +
  scale_x_continuous(name = "No. total repeats") +
  scale_y_log10(name = "Log10, Frequency") +
  ggtitle("CD8 alpha") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

##cd8 beta
tcr.cd8b.hist <- ggplot(tcr.cd8b, aes(x = count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "white") +
  scale_x_continuous(name = "No. total repeats") +
  scale_y_continuous(name = "Frequency") +
  ggtitle("CD8 beta") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

tcr.cd8b.hist2 <- ggplot(tcr.cd8b, aes(x = count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "white") +
  scale_x_continuous(name = "No. total repeats") +
  scale_y_log10(name = "Log10, Frequency") +
  ggtitle("CD8 beta") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

##Save plots
ggsave("../../Results/Exploratory_plots/CD4_alpha_tcr_freq.png", tcr.cd4a.hist, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/CD4_alpha_tcr_freq_log.png", tcr.cd4a.hist2, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/CD4_beta_tcr_freq.png", tcr.cd4b.hist, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/CD4_beta_tcr_freq_log.png", tcr.cd4b.hist2, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/CD8_alpha_tcr_freq.png", tcr.cd8a.hist, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/CD8_alpha_tcr_freq_log.png", tcr.cd8a.hist2, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/CD8_beta_tcr_freq.png", tcr.cd8b.hist, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/CD8_beta_tcr_freq_log.png", tcr.cd8b.hist2, bg = "transparent", units = "cm", width = 30, height = 20)

####Plot shared aa sequence results####
cd4a.heatmap <- ggplot(rep.cd4a.long, aes(x = as.numeric(Var1), y = as.numeric(Var2), fill = Freq)) +
  geom_tile() +
  scale_x_continuous(name = "", breaks = c(1:62), labels = cd4.names) +
  scale_y_continuous(name = "", breaks = c(1:62), labels = cd4.names) +
  scale_fill_gradientn(colours = c("#ffffd4", "#fed98e", "#fe9929", "#d95f0e", "#993404"), name = "Shared\nsequences", limits = c(110, 12000)) +
  ggtitle("CD4 alpha") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 20, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 20, color = "black", angle = 90, vjust = 0.5, hjust = 1),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_fixed()

cd4b.heatmap <- ggplot(rep.cd4b.long, aes(x = as.numeric(Var1), y = as.numeric(Var2), fill = Freq)) +
  geom_tile() +
  scale_x_continuous(name = "", breaks = c(1:62), labels = cd4.names) +
  scale_y_continuous(name = "", breaks = c(1:62), labels = cd4.names) +
  scale_fill_gradientn(colours = c("#ffffd4", "#fed98e", "#fe9929", "#d95f0e", "#993404"), name = "Shared\nsequences", limits = c(110, 12000)) +
  ggtitle("CD4 beta") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 20, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 20, color = "black", angle = 90, vjust = 0.5, hjust = 1),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_fixed()

cd8a.heatmap <- ggplot(rep.cd8a.long, aes(x = as.numeric(Var1), y = as.numeric(Var2), fill = Freq)) +
  geom_tile() +
  scale_x_continuous(name = "", breaks = c(1:58), labels = cd8.names) +
  scale_y_continuous(name = "", breaks = c(1:58), labels = cd8.names) +
  scale_fill_gradientn(colours = c("#ffffd4", "#fed98e", "#fe9929", "#d95f0e", "#993404"), name = "Shared\nsequences", limits = c(110, 12000)) +
  ggtitle("CD8 alpha") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 20, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 20, color = "black", angle = 90, vjust = 0.5, hjust = 1),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_fixed()

cd8b.heatmap <- ggplot(rep.cd8b.long, aes(x = as.numeric(Var1), y = as.numeric(Var2), fill = Freq)) +
  geom_tile() +
  scale_x_continuous(name = "", breaks = c(1:58), labels = cd8.names) +
  scale_y_continuous(name = "", breaks = c(1:58), labels = cd8.names) +
  scale_fill_gradientn(colours = c("#ffffd4", "#fed98e", "#fe9929", "#d95f0e", "#993404"), name = "Shared\nsequences", limits = c(110, 12000)) +
  ggtitle("CD8 beta") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 20, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 20, color = "black", angle = 90, vjust = 0.5, hjust = 1),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_fixed()

ggsave("../../Results/Exploratory_plots/CD4_alpha_shared_sequences.png", cd4a.heatmap, bg = "transparent", units = "cm", width = 60, height = 40)
ggsave("../../Results/Exploratory_plots/CD4_beta_shared_sequences.png", cd4b.heatmap, bg = "transparent", units = "cm", width = 60, height = 40)
ggsave("../../Results/Exploratory_plots/CD8_alpha_shared_sequences.png", cd8a.heatmap, bg = "transparent", units = "cm", width = 60, height = 40)
ggsave("../../Results/Exploratory_plots/CD8_beta_shared_sequences.png", cd8b.heatmap, bg = "transparent", units = "cm", width = 60, height = 40)

####Plot simpsons diversity index####
simpsons.plot <- ggplot(simpsons, aes(x = type, y = log(1/as.numeric(value)), color = type)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("CD4 alpha", "CD4 beta", "CD8 alpha", "CD8 beta")) +
  scale_y_continuous(name = "Log(1/Simpson's)") +
  ggtitle("") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

#separate by cd4a/b and cd8a/b
simpsons.cd4a.plot <- ggplot(simpsons[which(simpsons$type == "cd4a"), ], aes(x = site, y = log(1/as.numeric(value)), color = site)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("Lab", "Nott.", "Wirral")) +
  scale_y_continuous(name = "Log(1/Simpson's)") +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD4 alpha") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(6, 10.5))

simpsons.cd4b.plot <- ggplot(simpsons[which(simpsons$type == "cd4b"), ], aes(x = site, y = log(1/as.numeric(value)), color = site)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("Lab", "Nott.", "Wirral")) +
  scale_y_continuous(name = "Log(1/Simpson's)") +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD4 beta") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(6, 10.5))

simpsons.cd8a.plot <- ggplot(simpsons[which(simpsons$type == "cd8a"), ], aes(x = site, y = log(1/as.numeric(value)), color = site)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("Lab", "Nott.", "Wirral")) +
  scale_y_continuous(name = "Log(1/Simpson's)") +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD8 alpha") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(6, 10.5))

simpsons.cd8b.plot <- ggplot(simpsons[which(simpsons$type == "cd8b"), ], aes(x = site, y = log(1/as.numeric(value)), color = site)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("Lab", "Nott.", "Wirral")) +
  scale_y_continuous(name = "Log(1/Simpson's)") +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD8 beta") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(6, 10.5))

##by age
simpsons.age.plot <- ggplot(simpsons, aes(x = as.numeric(age), y = log(1/as.numeric(value)), color = type)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Log(1/Simpson's)") +
  ggtitle("") +
  scale_color_discrete(labels = c("CD4 alpha", "CD4 beta", "CD8 alpha", "CD8 beta"), name = "") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

#separate by cd4a/b and cd8a/b
simpsons.age.cd4a.plot <- ggplot(simpsons[which(simpsons$type == "cd4a"), ], aes(x = as.numeric(age), y = log(1/as.numeric(value)), color = site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Log(1/Simpson's)") +
  scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD4 alpha") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(6.25, 10.5))

simpsons.age.cd4b.plot <- ggplot(simpsons[which(simpsons$type == "cd4b"), ], aes(x = as.numeric(age), y = log(1/as.numeric(value)), color = site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Log(1/Simpson's)") +
  scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD4 beta") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(6.25, 10.5))

simpsons.age.cd8a.plot <- ggplot(simpsons[which(simpsons$type == "cd8a"), ], aes(x = as.numeric(age), y = log(1/as.numeric(value)), color = site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Log(1/Simpson's)") +
  scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD8 alpha") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(6.25, 10.5))

simpsons.age.cd8b.plot <- ggplot(simpsons[which(simpsons$type == "cd8b"), ], aes(x = as.numeric(age), y = log(1/as.numeric(value)), color = site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Log(1/Simpson's)") +
  scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD8 beta") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(6.25, 10.5))

##Save plots
ggsave("../../Results/Exploratory_plots/Simpsons_diversity.png", simpsons.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Simpsons_diversity_cd4a.png", simpsons.cd4a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Simpsons_diversity_cd4b.png", simpsons.cd4b.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Simpsons_diversity_cd8a.png", simpsons.cd8a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Simpsons_diversity_cd8b.png", simpsons.cd8b.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Simpsons_diversity_age.png", simpsons.age.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Simpsons_diversity_age_cd4a.png", simpsons.age.cd4a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Simpsons_diversity_age_cd4b.png", simpsons.age.cd4b.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Simpsons_diversity_age_cd8a.png", simpsons.age.cd8a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Simpsons_diversity_age_cd8b.png", simpsons.age.cd8b.plot, bg = "transparent", units = "cm", width = 30, height = 20)

####Plot shannons diversity index####
shannons.plot <- ggplot(shannons, aes(x = type, y = as.numeric(value), color = type)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("CD4 alpha", "CD4 beta", "CD8 alpha", "CD8 beta")) +
  scale_y_continuous(name = "Shannon's") +
  ggtitle("") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

#separate by cd4a/b and cd8a/b
shannons.cd4a.plot <- ggplot(shannons[which(shannons$type == "cd4a"), ], aes(x = site, y = as.numeric(value), color = site)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("Lab", "Nott.", "Wirral")) +
  scale_y_continuous(name = "Shannon's") +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD4 alpha") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(8, 10.5))

shannons.cd4b.plot <- ggplot(shannons[which(shannons$type == "cd4b"), ], aes(x = site, y = as.numeric(value), color = site)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("Lab", "Nott.", "Wirral")) +
  scale_y_continuous(name = "Shannon's") +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD4 beta") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(8, 10.5))

shannons.cd8a.plot <- ggplot(shannons[which(shannons$type == "cd8a"), ], aes(x = site, y = as.numeric(value), color = site)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("Lab", "Nott.", "Wirral")) +
  scale_y_continuous(name = "Shannon's") +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD8 alpha") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(8, 10.5))

shannons.cd8b.plot <- ggplot(shannons[which(shannons$type == "cd8b"), ], aes(x = site, y = as.numeric(value), color = site)) +
  geom_boxplot() +
  scale_x_discrete(name = "", labels = c("Lab", "Nott.", "Wirral")) +
  scale_y_continuous(name = "Shannon's") +
  scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD8 beta") +
  guides(color = "none") +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(8, 10.5))

##by age
shannons.age.plot <- ggplot(shannons, aes(x = as.numeric(age), y = as.numeric(value), color = type)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Shannon's") +
  ggtitle("") +
  scale_color_discrete(labels = c("CD4 alpha", "CD4 beta", "CD8 alpha", "CD8 beta"), name = "") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

#separate by cd4a/b and cd8a/b
shannons.age.cd4a.plot <- ggplot(shannons[which(shannons$type == "cd4a"), ], aes(x = as.numeric(age), y = as.numeric(value), color = site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Shannon's") +
  scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD4 alpha") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(8, 10.5))

shannons.age.cd4b.plot <- ggplot(shannons[which(shannons$type == "cd4b"), ], aes(x = as.numeric(age), y = as.numeric(value), color = site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Shannon's") +
  scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD4 beta") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(8, 10.5))

shannons.age.cd8a.plot <- ggplot(shannons[which(shannons$type == "cd8a"), ], aes(x = as.numeric(age), y = as.numeric(value), color = site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Shannon's") +
  scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD8 alpha") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(8, 10.5))

shannons.age.cd8b.plot <- ggplot(shannons[which(shannons$type == "cd8b"), ], aes(x = as.numeric(age), y = as.numeric(value), color = site)) +
  geom_point(size = 4) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Shannon's") +
  scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("CD8 beta") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(#panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    #legend.background = element_rect(fill='transparent'),
    plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm')) +
  coord_cartesian(ylim = c(8, 10.5))

#Save plots
ggsave("../../Results/Exploratory_plots/Shannons_diversity.png", shannons.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Shannons_diversity_cd4a.png", shannons.cd4a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Shannons_diversity_cd4b.png", shannons.cd4b.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Shannons_diversity_cd8a.png", shannons.cd8a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Shannons_diversity_cd8b.png", shannons.cd8b.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Shannons_diversity_age.png", shannons.age.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Shannons_diversity_age_cd4a.png", shannons.age.cd4a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Shannons_diversity_age_cd4b.png", shannons.age.cd4b.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Shannons_diversity_age_cd8a.png", shannons.age.cd8a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave("../../Results/Exploratory_plots/Shannons_diversity_age_cd8b.png", shannons.age.cd8b.plot, bg = "transparent", units = "cm", width = 30, height = 20)

####Testing network plotting####
nodes <- read.csv("../Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("../Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)

#rename columns
colnames(cd4a.net) <- c("from", "to", "weight")
#make weight column into numeric
cd4a.net$weight <- as.numeric(cd4a.net$weight)
#remove NA values
cd4a.net <- cd4a.net[-which(is.na(cd4a.net$weight)),]
#create data frame for nodes
cd4a.nodes <- as.data.frame(unique(cd4a.net$from))
colnames(cd4a.nodes) <- c("id")
cd4a.nodes$location <- c(rep("Lab", 5), rep("Nottingham", 52), rep("Wirral", 4), "Lab")
cd4a.nodes$location.num <- c(rep(1, 5), rep(2, 52), rep(3, 4), 1)

net <- graph_from_data_frame(d = cd4a.net, vertices = cd4a.nodes, directed = F)

#Set edge width based on weight:
E(net)$width <- E(net)$weight/1000

#change arrow size and edge color:
E(net)$arrow.size <- 0.2
E(net)$edge.color <- "gray80"

#Generate colors based on media type:
colrs <- c("gray50", "tomato", "gold")
V(net)$color <- colrs[V(net)$location.num]

#Setting them to NA will render no labels:
V(net)$label <- NA

plot(net, layout = layout_with_fr, edge.width = E(net)$width, edge.color=adjustcolor(1, alpha.f = 0.15))
