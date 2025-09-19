# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Project and code information ------------------------------------------------
##
## Project: Wild Mouse TCR
##
## Purpose of script: Explore and visualise data
##
## Author: Jacob Cohen
##
## Date Created: 2025-05-09
##
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

## Clear workspace -------------------------------------------------------------
message("Clear workspace")

rm(list=ls())


## Load packages ---------------------------------------------------------------
message("Load packages")

packages <- c('here', 'ggplot2', 'stringr', 'ggbeeswarm', 'emmeans', 'showtext', 'patchwork',
              'DHARMa', 'dplyr', 'tidyr', 'vegan', 'ggpubr', 'gplots', 'modeest', 'scales')
lapply(packages, library, character.only = TRUE)

#additional code for showtext package
font_add(family = "Calibri", regular = "C:/Windows/Fonts/calibri.ttf")
showtext_auto()
showtext_opts(dpi = 300)


## Load data -------------------------------------------------------------------
message("Load data")

#read in updated mouse metadata file
mouse.data <- read.csv(here("Data", "Mouse_data.csv"), header = TRUE)

#read in repertoire data
rep.cd4a <- read.csv(here("Data", "CD4_alpha_Repertoire.csv"), header = TRUE)
rep.cd4b <- read.csv(here("Data", "CD4_beta_Repertoire.csv"), header = TRUE)
rep.cd8a <- read.csv(here("Data", "CD8_alpha_Repertoire.csv"), header = TRUE)
rep.cd8b <- read.csv(here("Data", "CD8_beta_Repertoire.csv"), header = TRUE)

#read in TCR frequency data
tcr.cd4a <- read.csv(here("Data", "CD4_alpha_TCR_freq.csv"), header = TRUE)
tcr.cd4b <- read.csv(here("Data", "CD4_beta_TCR_freq.csv"), header = TRUE)
tcr.cd8a <- read.csv(here("Data", "CD8_alpha_TCR_freq.csv"), header = TRUE)
tcr.cd8b <- read.csv(here("Data", "CD8_beta_TCR_freq.csv"), header = TRUE)

#read in PCOA data
share.cd4a <- read.csv(here("Data", "CD4_alpha_Repertoire_PCOA.csv"), header = TRUE)
share.cd4b <- read.csv(here("Data", "CD4_beta_Repertoire_PCOA.csv"), header = TRUE)
share.cd8a <- read.csv(here("Data", "CD8_alpha_Repertoire_PCOA.csv"), header = TRUE)
share.cd8b <- read.csv(here("Data", "CD8_beta_Repertoire_PCOA.csv"), header = TRUE)


## Clean data ------------------------------------------------------------------
message("Clean data")

#keep version which includes the lab mice
mouse.data.lab <- mouse.data

#use version without the lab mice
mouse.data <- mouse.data[-c(1:6),]

#remove unneeded mice from repertoire data
rep.cd4a <- subset(rep.cd4a, select = -c(N63, N42, N24, N01, N27, N11, N12, N09, N35),
                   !(rep.cd4a$X %in% c("N63", "N42", "N24", "N01", "N27", "N11", "N12", "N09", "N35")))
rep.cd4b <- subset(rep.cd4b, select = -c(N63, N42, N24, N01, N27, N11, N12, N09, N35),
                   !(rep.cd4b$X %in% c("N63", "N42", "N24", "N01", "N27", "N11", "N12", "N09", "N35")))
rep.cd8a <- subset(rep.cd8a, select = -c(N17, N42, N22, N01, N21, N11, N09, N27, N12, N51, N10, N35, N36),
                   !(rep.cd8a$X %in% c("N17", "N42", "N22", "N01", "N21", "N11", "N09", "N27", "N12", "N51", "N10", "N35", "N36")))
rep.cd8b <- subset(rep.cd8b, select = -c(N17, N42, N22, N01, N21, N11, N09, N27, N12, N51, N10, N35, N36),
                   !(rep.cd8b$X %in% c("N17", "N42", "N22", "N01", "N21", "N11", "N09", "N27", "N12", "N51", "N10", "N35", "N36")))

#save these data frames for heatmap plotting
heat.cd4a <- rep.cd4a
rownames(heat.cd4a) <- heat.cd4a[, 1]
heat.cd4a <- heat.cd4a[, -1]
heat.cd4a <- as.matrix(heat.cd4a)

heat.cd4b <- rep.cd4b
rownames(heat.cd4b) <- heat.cd4b[, 1]
heat.cd4b <- heat.cd4b[, -1]
heat.cd4b <- as.matrix(heat.cd4b)

heat.cd8a <- rep.cd8a
rownames(heat.cd8a) <- heat.cd8a[, 1]
heat.cd8a <- heat.cd8a[, -1]
heat.cd8a <- as.matrix(heat.cd8a)

heat.cd8b <- rep.cd8b
rownames(heat.cd8b) <- heat.cd8b[, 1]
heat.cd8b <- heat.cd8b[, -1]
heat.cd8b <- as.matrix(heat.cd8b)

# #change from matrix format to long format
# cd4a.net <- as.data.frame.table(`row.names<-`(as.matrix(rep.cd4a)[,-1],
#                                               rep.cd4a$X))
# cd4b.net <- as.data.frame.table(`row.names<-`(as.matrix(rep.cd4b)[,-1],
#                                               rep.cd4b$X))
# cd8a.net <- as.data.frame.table(`row.names<-`(as.matrix(rep.cd8a)[,-1],
#                                               rep.cd8a$X))
# cd8b.net <- as.data.frame.table(`row.names<-`(as.matrix(rep.cd8b)[,-1],
#                                               rep.cd8b$X))

#get mouse IDs
cd4.names <- rep.cd4a$X
cd8.names <- rep.cd8a$X

#drop mouse IDs from the heatmap data
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

#remove unneeded data
rm(rep.cd4a, rep.cd4b, rep.cd8a, rep.cd8b)


## Calculate TCR repertoire richness -------------------------------------------
message("Calculate TCR repertoire richness")

##add columns for unique seq/repertoire size
mouse.data$CD4.alpha.unique.div.rep <- mouse.data$CD4.alpha.unique.seq / mouse.data$CD4.alpha.count
mouse.data$CD4.beta.unique.div.rep <- mouse.data$CD4.beta.unique.seq / mouse.data$CD4.beta.count
mouse.data$CD8.alpha.unique.div.rep <- mouse.data$CD8.alpha.unique.seq / mouse.data$CD8.alpha.count
mouse.data$CD8.beta.unique.div.rep <- mouse.data$CD8.beta.unique.seq / mouse.data$CD8.beta.count

#data for unique unique seq/repertoire size
unique.data <- as.data.frame(rbind(cbind("cd4a", mouse.data$CD4.alpha.unique.div.rep, mouse.data$Age_days, mouse.data$ID, mouse.data$Site, mouse.data$Sex),
                                   cbind("cd4b", mouse.data$CD4.beta.unique.div.rep, mouse.data$Age_days, mouse.data$ID, mouse.data$Site, mouse.data$Sex),
                                   cbind("cd8a", mouse.data$CD8.alpha.unique.div.rep, mouse.data$Age_days, mouse.data$ID, mouse.data$Site, mouse.data$Sex),
                                   cbind("cd8b", mouse.data$CD8.beta.unique.div.rep, mouse.data$Age_days, mouse.data$ID, mouse.data$Site, mouse.data$Sex)))

#adjust structure of the data
colnames(unique.data) <- c("type", "value", "age", "id", "site", "sex")
unique.data$type <- as.factor(unique.data$type)
unique.data$value <- as.numeric(unique.data$value)
unique.data$age <- as.numeric(unique.data$age)
unique.data$id <- as.factor(unique.data$id)
unique.data$site <- as.factor(unique.data$site)
unique.data$sex <- as.factor(unique.data$sex)


## Generate diversity metric data frame ----------------------------------------
message("Generate diversity metric data frame")

##data for diversity data analysis with lab mice
div.data.lab <- as.data.frame(rbind(cbind("cd4a", mouse.data.lab$CD4.alpha.simpsons, mouse.data.lab$CD4.alpha.shannons, mouse.data.lab$Age_days, mouse.data.lab$ID, mouse.data.lab$Site, mouse.data.lab$Sex),
                                    cbind("cd4b", mouse.data.lab$CD4.beta.simpsons, mouse.data.lab$CD4.beta.shannons, mouse.data.lab$Age_days, mouse.data.lab$ID, mouse.data.lab$Site, mouse.data.lab$Sex),
                                    cbind("cd8a", mouse.data.lab$CD8.alpha.simpsons, mouse.data.lab$CD8.alpha.shannons, mouse.data.lab$Age_days, mouse.data.lab$ID, mouse.data.lab$Site, mouse.data.lab$Sex),
                                    cbind("cd8b", mouse.data.lab$CD8.beta.simpsons, mouse.data.lab$CD8.beta.shannons, mouse.data.lab$Age_days, mouse.data.lab$ID, mouse.data.lab$Site, mouse.data.lab$Sex)))

#adjust structure of the data
colnames(div.data.lab) <- c("type", "simpsons", "shannons", "age", "id", "site", "sex")
div.data.lab$type <- as.factor(div.data.lab$type)
div.data.lab$simpsons <- as.numeric(div.data.lab$simpsons)
div.data.lab$shannons <- as.numeric(div.data.lab$shannons)
div.data.lab$age <- as.numeric(div.data.lab$age)
div.data.lab$id <- as.factor(div.data.lab$id)
div.data.lab$site <- as.factor(div.data.lab$site)
div.data.lab$sex <- as.factor(div.data.lab$sex)

##data for diversity data analysis without lab mice
div.data <- as.data.frame(rbind(cbind("cd4a", mouse.data$CD4.alpha.simpsons, mouse.data$CD4.alpha.shannons, mouse.data$Age_days, mouse.data$ID, mouse.data$Site, mouse.data$Sex),
                                cbind("cd4b", mouse.data$CD4.beta.simpsons, mouse.data$CD4.beta.shannons, mouse.data$Age_days, mouse.data$ID, mouse.data$Site, mouse.data$Sex),
                                cbind("cd8a", mouse.data$CD8.alpha.simpsons, mouse.data$CD8.alpha.shannons, mouse.data$Age_days, mouse.data$ID, mouse.data$Site, mouse.data$Sex),
                                cbind("cd8b", mouse.data$CD8.beta.simpsons, mouse.data$CD8.beta.shannons, mouse.data$Age_days, mouse.data$ID, mouse.data$Site, mouse.data$Sex)))

#adjust structure of the data
colnames(div.data) <- c("type", "simpsons", "shannons", "age", "id", "site", "sex")
div.data$type <- as.factor(div.data$type)
div.data$simpsons <- as.numeric(div.data$simpsons)
div.data$shannons <- as.numeric(div.data$shannons)
div.data$age <- as.numeric(div.data$age)
div.data$id <- as.factor(div.data$id)
div.data$site <- as.factor(div.data$site)
div.data$sex <- as.factor(div.data$sex)


## Statistical tests for TCR repertoire richness -----------------------------------
message("Statistical tests TCR repertoire richness")

#glm
unique.mod <- glm(value ~ site + age + type + sex + age*sex, data = unique.data, family = "Gamma")
summary(unique.mod)

#test fit
unique.output <- simulateResiduals(fittedModel = unique.mod)
plot(unique.output)
#DHARMa says the fit is good

#Post hoc test for pairwise receptor type comparisons
unique.post <- emmeans(unique.mod, specs = pairwise ~ type, weights = "proportional")
summary(unique.post)

#Post hoc test for pairwise site comparisons
unique.post2 <- emmeans(unique.mod, specs = pairwise ~ site, weights = "proportional")
summary(unique.post2)

## Plot TCR repertoire richness ------------------------------------------------
message("Plot TCR repertoire richness")

unique.plot <- ggplot(unique.data, aes(x = type, y = as.numeric(value), color = type)) +
  geom_boxplot(linewidth = 2, outlier.size = 3) +
  scale_x_discrete(name = "", labels = c("CD4 alpha", "CD4 beta", "CD8 alpha", "CD8 beta")) +
  scale_y_continuous(name = "TCR richness") +
  geom_bracket(xmin = c(1, 2, 3, 1), xmax = c(2, 3, 4, 4), label = "*", y.position = c(0.55, 0.515, 0.55, 0.6), color = "black", inherit.aes = FALSE,
               size = 1, label.size = 15, tip.length = 0.01, vjust = 0.5) +
  ggtitle("") +
  guides(color = "none") +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 0.65))

ggsave(here("Results", "Updated_plots", "Final_plots", "TCR_richness.png"), unique.plot, bg = "transparent", units = "cm", width = 30, height = 20)


## Exploratory plots for TCR repertoire richness -------------------------------
message("Exploratory plots for TCR repertoire richness")

# unique.cd4a.plot <- ggplot(unique.data[which(unique.data$type == "cd4a"), ], aes(x = site, y = as.numeric(value), color = site, shape = sex)) +
#   geom_beeswarm(size = 4, cex = 3) +
#   scale_x_discrete(name = "", labels = c("Nott.", "Wirral")) +
#   scale_y_continuous(name = "TCR repertoire richness") +
#   scale_color_manual(values = c('#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD4 alpha") +
#   guides(color = "none", shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(0, 0.5))
# 
# unique.cd4b.plot <- ggplot(unique.data[which(unique.data$type == "cd4b"), ], aes(x = site, y = as.numeric(value), color = site, shape = sex)) +
#   geom_beeswarm(size = 4, cex = 3) +
#   scale_x_discrete(name = "", labels = c("Nott.", "Wirral")) +
#   scale_y_continuous(name = "TCR repertoire richness") +
#   scale_color_manual(values = c('#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD4 beta") +
#   guides(color = "none", shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(0, 0.5))
# 
# unique.cd8a.plot <- ggplot(unique.data[which(unique.data$type == "cd8a"), ], aes(x = site, y = as.numeric(value), color = site, shape = sex)) +
#   geom_beeswarm(size = 4, cex = 3) +
#   scale_x_discrete(name = "", labels = c("Nott.", "Wirral")) +
#   scale_y_continuous(name = "TCR repertoire richness") +
#   scale_color_manual(values = c('#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD8 alpha") +
#   guides(color = "none", shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(0, 0.5))
# 
# unique.cd8b.plot <- ggplot(unique.data[which(unique.data$type == "cd8b"), ], aes(x = site, y = as.numeric(value), color = site, shape = sex)) +
#   geom_beeswarm(size = 4, cex = 3) +
#   scale_x_discrete(name = "", labels = c("Nott.", "Wirral")) +
#   scale_y_continuous(name = "TCR repertoire richness") +
#   scale_color_manual(values = c('#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD8 beta") +
#   guides(color = "none", shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(0, 0.5))
# 
# ggsave(here("Results", "Updated_plots", "Unique_aa_div_rep_size_cd4a.png"), unique.cd4a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "Unique_aa_div_rep_size_cd4b.png"), unique.cd4b.plot, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "Unique_aa_div_rep_size_cd8a.png"), unique.cd8a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "Unique_aa_div_rep_size_cd8b.png"), unique.cd8b.plot, bg = "transparent", units = "cm", width = 30, height = 20)


## Exploratory plots for powerlaw distribution ---------------------------------
message("xploratory plots for powerlaw distribution")

CD4.alpha.powerlaw.plot <- ggplot(data = mouse.data, aes(x = Age_days, y = CD4.alpha.powerlaw, color = Site)) +
  geom_point(size = 8) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "CD4 alpha power law coefficient") +
  scale_color_manual(name = "Site", labels = c("Nottingham", "Wirral"),
                     values = c('#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.key.height= unit(1.5, 'cm'),
    legend.key.width= unit(1.5, 'cm'),
    legend.position = "bottom", 
    legend.justification = "left",
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 6), xlim = c(0, 350))

CD4.beta.powerlaw.plot <- ggplot(data = mouse.data, aes(x = Age_days, y = CD4.beta.powerlaw, color = Site)) +
  geom_point(size = 8) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "CD4 beta power law coefficient") +
  scale_color_manual(name = "Site", labels = c("Nottingham", "Wirral"),
                     values = c('#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.key.height= unit(1.5, 'cm'),
    legend.key.width= unit(1.5, 'cm'),
    legend.position = "bottom", 
    legend.justification = "left",
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 6), xlim = c(0, 350))

CD8.alpha.powerlaw.plot <- ggplot(data = mouse.data, aes(x = Age_days, y = CD8.alpha.powerlaw, color = Site)) +
  geom_point(size = 8) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "CD8 alpha power law coefficient") +
  scale_color_manual(name = "Site", labels = c("Nottingham", "Wirral"),
                     values = c('#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.key.height= unit(1.5, 'cm'),
    legend.key.width= unit(1.5, 'cm'),
    legend.position = "bottom", 
    legend.justification = "left",
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 6), xlim = c(0, 350))

CD8.beta.powerlaw.plot <- ggplot(data = mouse.data, aes(x = Age_days, y = CD8.beta.powerlaw, color = Site)) +
  geom_point(size = 8) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "CD8 beta power law coefficient") +
  scale_color_manual(name = "Site", labels = c("Nottingham", "Wirral"),
                     values = c('#d95f02', '#7570b3')) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.key.height= unit(1.5, 'cm'),
    legend.key.width= unit(1.5, 'cm'),
    legend.position = "bottom", 
    legend.justification = "left",
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 6), xlim = c(0, 350))

# #Save plots
# ggsave(here("Results", "Updated_plots", "CD4_alpha_power_law_coefficient.png"), CD4.alpha.powerlaw.plot, bg = "transparent", units = "cm", width = 30, height = 23)
# ggsave(here("Results", "Updated_plots", "CD4_beta_power_law_coefficient.png"), CD4.beta.powerlaw.plot, bg = "transparent", units = "cm", width = 30, height = 23)
# ggsave(here("Results", "Updated_plots", "CD8_alpha_power_law_coefficient.png"), CD8.alpha.powerlaw.plot, bg = "transparent", units = "cm", width = 30, height = 23)
# ggsave(here("Results", "Updated_plots", "CD8_beta_power_law_coefficient.png"), CD8.beta.powerlaw.plot, bg = "transparent", units = "cm", width = 30, height = 23)

#patchwork plots
cd4.powerlaw.plot <- CD4.alpha.powerlaw.plot + CD4.beta.powerlaw.plot +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(plot.tag = element_text(size = 42, face = "bold"), legend.position = "bottom", legend.justification = "left"))

cd8.powerlaw.plot <- CD8.alpha.powerlaw.plot + CD8.beta.powerlaw.plot +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(plot.tag = element_text(size = 42, face = "bold"), legend.position = "bottom", legend.justification = "left"))

#save patchwork
ggsave(here("Results", "Updated_plots", "Final_plots", "CD4_powerlaw_plot.png"), cd4.powerlaw.plot, bg = "transparent", units = "cm", width = 60, height = 30)
ggsave(here("Results", "Updated_plots", "Final_plots", "CD8_powerlaw_plot.png"), cd8.powerlaw.plot, bg = "transparent", units = "cm", width = 60, height = 30)


## Statistical tests for Simpson's diversity -----------------------------------
message("Statistical tests for Simpson's diversity")

##EXCLUDING LAB MICE
simpsons.mod <- glm(log(1/simpsons) ~ site + age + type + sex + age*sex, data = div.data, family = "gaussian")
summary(simpsons.mod)

#test fit
simpsons.output <- simulateResiduals(fittedModel = simpsons.mod)
plot(simpsons.output)
#DHARMa says the fit is not great, but better than a gamma distribution (and also has a better AIC)

#Post hoc test for pairwise receptor type comparisons
simpsons.post <- emmeans(simpsons.mod, specs = pairwise ~ type, weights = "proportional")
summary(simpsons.post)

#Post hoc test for pairwise site comparisons
simpsons.post2 <- emmeans(simpsons.mod, specs = pairwise ~ site, weights = "proportional")
summary(simpsons.post2)

#Post hoc test for pairwise site + type comparisons
simpsons.post3 <- emmeans(simpsons.mod, specs = pairwise ~ site + type, weights = "proportional")
summary(simpsons.post3)

#Post hoc test for sex comparison
simpsons.post4 <- emmeans(simpsons.mod, specs = pairwise ~ sex, weights = "proportional")
summary(simpsons.post4)
#marginally non-significant


##INCLUDING LAB MICE
simpsons.mod.lab <- glm(log(1/simpsons) ~ site + age + type + sex + age*sex, data = div.data.lab, family = "gaussian")
summary(simpsons.mod.lab)

#test fit
simpsons.output.lab <- simulateResiduals(fittedModel = simpsons.mod.lab)
plot(simpsons.output.lab)
#DHARMa says the fit is not great, but better than a gamma distribution (and also has a better AIC)

#Post hoc test for pairwise receptor type comparisons
simpsons.post.lab <- emmeans(simpsons.mod.lab, specs = pairwise ~ type, weights = "proportional")
summary(simpsons.post.lab)

#Post hoc test for pairwise site comparisons
simpsons.post2.lab <- emmeans(simpsons.mod.lab, specs = pairwise ~ site, weights = "proportional")
summary(simpsons.post2.lab)

#Post hoc test for pairwise site + type comparisons
simpsons.post3.lab <- emmeans(simpsons.mod.lab, specs = pairwise ~ site + type, weights = "proportional")
summary(simpsons.post3.lab)

#Post hoc test for sex comparison
simpsons.post4.lab <- emmeans(simpsons.mod.lab, specs = pairwise ~ sex, weights = "proportional")
summary(simpsons.post4.lab)


## Plot Simpson's diversity ----------------------------------------------------
message("Plot Simpson's diversity")

##EXCLUDING LAB MICE
simpsons.plot <- ggplot(div.data, aes(x = type, y = log(1/simpsons), color = type)) +
  geom_boxplot(linewidth = 2, outlier.size = 3) +
  scale_x_discrete(name = "", labels = c("CD4 alpha", "CD4 beta", "CD8 alpha", "CD8 beta")) +
  scale_y_continuous(name = "Simpson's index") +
  geom_bracket(xmin = c(1, 2, 2, 1), xmax = c(3, 3, 4, 4), label = "*", y.position = c(10.7, 10.1, 11.3, 11.9), color = "black", inherit.aes = FALSE,
               size = 1, label.size = 15, tip.length = 0.01, vjust = 0.5) +
  ggtitle("") +
  guides(color = "none") +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(5.8, 11.9))

ggsave(here("Results", "Updated_plots", "Final_plots", "Simpsons_diversity.png"), simpsons.plot, bg = "transparent", units = "cm", width = 30, height = 20)


## Exploratory plots for Simpson's diversity -----------------------------------
message("Exploratory plots for Simpson's diversity")

# ##INCLUDING LAB MICE
# simpsons.plot.lab <- ggplot(div.data.lab, aes(x = type, y = log(1/simpsons), color = type)) +
#   geom_boxplot(linewidth = 2, outlier.size = 3) +
#   scale_x_discrete(name = "", labels = c("CD4 alpha", "CD4 beta", "CD8 alpha", "CD8 beta")) +
#   scale_y_continuous(name = "Log(1/Simpson's)") +
#   geom_bracket(xmin = c(1, 2, 2), xmax = c(3, 3, 4), label = "*", y.position = c(10.1, 10.7, 11.3), color = "black", inherit.aes = FALSE,
#                size = 1, label.size = 15, tip.length = 0.01, vjust = 0.5) +
#   ggtitle("") +
#   guides(color = "none") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(5.8, 11.3))
# 
# ggsave(here("Results", "Updated_plots", "Simpsons_diversity_lab.png"), simpsons.plot.lab, bg = "transparent", units = "cm", width = 30, height = 20)
# 
# #separate by cd4a/b and cd8a/b
# #EXCLUDING LAB MICE
# simpsons.age.cd4a.plot <- ggplot(div.data[which(div.data$type == "cd4a"), ], aes(x = age, y = log(1/simpsons), color = site, shape = sex)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Age (days)") +
#   scale_y_continuous(name = "Log(1/Simpson's)") +
#   scale_color_manual(labels = c("Nott.", "Wirral"), name = "", values = c('#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD4 alpha") +
#   guides(color = guide_legend(override.aes = list(size = 8)), shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(6.25, 10.5))
# 
# simpsons.age.cd4b.plot <- ggplot(div.data[which(div.data$type == "cd4b"), ], aes(x = age, y = log(1/simpsons), color = site, shape = sex)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Age (days)") +
#   scale_y_continuous(name = "Log(1/Simpson's)") +
#   scale_color_manual(labels = c("Nott.", "Wirral"), name = "", values = c('#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD4 beta") +
#   guides(color = guide_legend(override.aes = list(size = 8)), shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(6.25, 10.5))
# 
# simpsons.age.cd8a.plot <- ggplot(div.data[which(div.data$type == "cd8a"), ], aes(x = age, y = log(1/simpsons), color = site, shape = sex)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Age (days)") +
#   scale_y_continuous(name = "Log(1/Simpson's)") +
#   scale_color_manual(labels = c("Nott.", "Wirral"), name = "", values = c('#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD8 alpha") +
#   guides(color = guide_legend(override.aes = list(size = 8)), shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(6.25, 10.5))
# 
# simpsons.age.cd8b.plot <- ggplot(div.data[which(div.data$type == "cd8b"), ], aes(x = age, y = log(1/simpsons), color = site, shape = sex)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Age (days)") +
#   scale_y_continuous(name = "Log(1/Simpson's)") +
#   scale_color_manual(labels = c("Nott.", "Wirral"), name = "", values = c('#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD8 beta") +
#   guides(color = guide_legend(override.aes = list(size = 8)), shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(6.25, 10.5))
# 
# ##Save plots
# ggsave(here("Results", "Updated_plots", "Simpsons_diversity_age_cd4a.png"), simpsons.age.cd4a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "Simpsons_diversity_age_cd4b.png"), simpsons.age.cd4b.plot, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "Simpsons_diversity_age_cd8a.png"), simpsons.age.cd8a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "Simpsons_diversity_age_cd8b.png"), simpsons.age.cd8b.plot, bg = "transparent", units = "cm", width = 30, height = 20)
# 
# #separate by cd4a/b and cd8a/b
# #INCLUDING LAB MICE
# simpsons.age.cd4a.plot.lab <- ggplot(div.data.lab[which(div.data.lab$type == "cd4a"), ], aes(x = age, y = log(1/simpsons), color = site, shape = sex)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Age (days)") +
#   scale_y_continuous(name = "Log(1/Simpson's)") +
#   scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD4 alpha") +
#   guides(color = guide_legend(override.aes = list(size = 8)), shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(6.25, 10.5))
# 
# simpsons.age.cd4b.plot.lab <- ggplot(div.data.lab[which(div.data.lab$type == "cd4b"), ], aes(x = age, y = log(1/simpsons), color = site, shape = sex)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Age (days)") +
#   scale_y_continuous(name = "Log(1/Simpson's)") +
#   scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD4 beta") +
#   guides(color = guide_legend(override.aes = list(size = 8)), shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(6.25, 10.5))
# 
# simpsons.age.cd8a.plot.lab <- ggplot(div.data.lab[which(div.data.lab$type == "cd8a"), ], aes(x = age, y = log(1/simpsons), color = site, shape = sex)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Age (days)") +
#   scale_y_continuous(name = "Log(1/Simpson's)") +
#   scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD8 alpha") +
#   guides(color = guide_legend(override.aes = list(size = 8)), shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(6.25, 10.5))
# 
# simpsons.age.cd8b.plot.lab <- ggplot(div.data.lab[which(div.data.lab$type == "cd8b"), ], aes(x = age, y = log(1/simpsons), color = site, shape = sex)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Age (days)") +
#   scale_y_continuous(name = "Log(1/Simpson's)") +
#   scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD8 beta") +
#   guides(color = guide_legend(override.aes = list(size = 8)), shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(6.25, 10.5))
# 
# ##Save plots
# ggsave(here("Results", "Updated_plots", "Simpsons_diversity_age_cd4a_lab.png"), simpsons.age.cd4a.plot.lab, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "Simpsons_diversity_age_cd4b_lab.png"), simpsons.age.cd4b.plot.lab, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "Simpsons_diversity_age_cd8a_lab.png"), simpsons.age.cd8a.plot.lab, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "Simpsons_diversity_age_cd8b_lab.png"), simpsons.age.cd8b.plot.lab, bg = "transparent", units = "cm", width = 30, height = 20)


## Statistical tests for Shannon's diversity -----------------------------------
message("Statistical test for Shannon's diversity")

##EXCLUDING LAB MICE
shannons.mod <- glm(shannons ~ site + age + type + sex + age*sex, data = div.data, family = "gaussian")
summary(shannons.mod)

#test fit
shannons.output <- simulateResiduals(fittedModel = shannons.mod)
plot(shannons.output)
#DHARMa says the fit is good

#Post hoc test for pairwise receptor type comparisons
shannons.post <- emmeans(shannons.mod, specs = pairwise ~ type, weights = "proportional")
summary(shannons.post)

#Post hoc test for pairwise site comparisons
shannons.post2 <- emmeans(shannons.mod, specs = pairwise ~ site, weights = "proportional")
summary(shannons.post2)

#Post hoc test for pairwise site + type comparisons
shannons.post3 <- emmeans(shannons.mod, specs = pairwise ~ site + type, weights = "proportional")
summary(shannons.post3)

#Post hoc test for sex comparison
shannons.post4 <- emmeans(shannons.mod, specs = pairwise ~ sex, weights = "proportional")
summary(shannons.post4)
#marginally non-signficant 

##INCLUDING LAB MICE
shannons.mod.lab <- glm(shannons ~ site + age + type + sex, data = div.data.lab, family = "gaussian")
summary(shannons.mod.lab)

#test fit
shannons.output.lab <- simulateResiduals(fittedModel = shannons.mod.lab)
plot(shannons.output.lab)
#DHARMa says the fit is good

#Post hoc test for pairwise receptor type comparisons
shannons.post.lab <- emmeans(shannons.mod.lab, specs = pairwise ~ type, weights = "proportional")
summary(shannons.post.lab)

#Post hoc test for pairwise site comparisons
shannons.post2.lab <- emmeans(shannons.mod.lab, specs = pairwise ~ site, weights = "proportional")
summary(shannons.post2.lab)

#Post hoc test for pairwise site + type comparisons
shannons.post3.lab <- emmeans(shannons.mod.lab, specs = pairwise ~ site + type, weights = "proportional")
summary(shannons.post3.lab)

#Post hoc test for sex comparison
shannons.post4.lab <- emmeans(shannons.mod.lab, specs = pairwise ~ sex, weights = "proportional")
summary(shannons.post4.lab)
#marginally non-signficant 


## Plot Shannon's diversity ----------------------------------------------------
message("Plot Shannon's diversity")

##EXCLUDING LAB MICE
shannons.plot <- ggplot(div.data, aes(x = type, y = shannons, color = type)) +
  geom_boxplot(linewidth = 2, outlier.size = 3) +
  scale_x_discrete(name = "", labels = c("CD4 alpha", "CD4 beta", "CD8 alpha", "CD8 beta")) +
  scale_y_continuous(name = "Shannon's index") +
  geom_bracket(xmin = c(2), xmax = c(3), label = "*", y.position = c(10.4), color = "black", inherit.aes = FALSE,
               size = 1, label.size = 15, tip.length = 0.01, vjust = 0.5) +
  ggtitle("") +
  guides(color = "none") +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(8, 10.4))

ggsave(here("Results", "Updated_plots", "Final_plots", "Shannons_diversity.png"), shannons.plot, bg = "transparent", units = "cm", width = 30, height = 20)


#separate by cd4a/b and cd8a/b
shannons.age.cd4a.plot <- ggplot(div.data[which(div.data$type == "cd4a"), ], aes(x = age, y = shannons, color = site, shape = sex)) +
  geom_point(size = 10, alpha = 0.8) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Shannon's index") +
  scale_color_manual(labels = c("Nottingham", "Wirral"), name = "Site", values = c('#d95f02', '#7570b3')) +
  scale_shape_manual(values = c("circle", "triangle"), name = "Sex") +
  ggtitle("") +
  guides(color = guide_legend(order = 1, override.aes = list(size = 10, alpha = 1)), shape = guide_legend(order = 2, override.aes = list(size = 10, fill = "black", alpha = 1))) +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.key.height= unit(1.5, 'cm'),
    legend.key.width= unit(1.5, 'cm'),
    legend.position = "bottom", 
    legend.justification = "left",
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(8, 10.5))

shannons.age.cd4b.plot <- ggplot(div.data[which(div.data$type == "cd4b"), ], aes(x = age, y = shannons, color = site, shape = sex)) +
  geom_point(size = 10, alpha = 0.8) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "") +
  scale_color_manual(labels = c("Nottingham", "Wirral"), name = "Site", values = c('#d95f02', '#7570b3')) +
  scale_shape_manual(values = c("circle", "triangle"), name = "Sex") +
  ggtitle("") +
  guides(color = guide_legend(order = 1, override.aes = list(size = 10, alpha = 1)), shape = guide_legend(order = 2, override.aes = list(size = 10, fill = "black", alpha = 1))) +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.key.height= unit(1.5, 'cm'),
    legend.key.width= unit(1.5, 'cm'),
    legend.position = "bottom", 
    legend.justification = "left",
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(8, 10.5))

shannons.age.cd4b.plot2 <- ggplot(div.data[which(div.data$type == "cd4b"), ], aes(x = age, y = shannons, color = site, shape = sex)) +
  geom_point(size = 10, alpha = 0.8) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Shannon's index") +
  scale_color_manual(labels = c("Nottingham", "Wirral"), name = "Site", values = c('#d95f02', '#7570b3')) +
  scale_shape_manual(values = c("circle", "triangle"), name = "Sex") +
  ggtitle("") +
  guides(color = guide_legend(order = 1, override.aes = list(size = 10, alpha = 1)), shape = guide_legend(order = 2, override.aes = list(size = 10, fill = "black", alpha = 1))) +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.key.height= unit(1.5, 'cm'),
    legend.key.width= unit(1.5, 'cm'),
    legend.position = "bottom", 
    legend.justification = "left",
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(8, 10.5))

shannons.age.cd8a.plot <- ggplot(div.data[which(div.data$type == "cd8a"), ], aes(x = age, y = shannons, color = site, shape = sex)) +
  geom_point(size = 10, alpha = 0.8) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "Shannon's index") +
  scale_color_manual(labels = c("Nottingham", "Wirral"), name = "Site", values = c('#d95f02', '#7570b3')) +
  scale_shape_manual(values = c("circle", "triangle"), name = "Sex") +
  ggtitle("") +
  guides(color = guide_legend(order = 1, override.aes = list(size = 10, alpha = 1)), shape = guide_legend(order = 2, override.aes = list(size = 10, fill = "black", alpha = 1))) +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.key.height= unit(1.5, 'cm'),
    legend.key.width= unit(1.5, 'cm'),
    legend.position = "bottom", 
    legend.justification = "left",
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(8, 10.5))

shannons.age.cd8b.plot <- ggplot(div.data[which(div.data$type == "cd8b"), ], aes(x = age, y = shannons, color = site, shape = sex)) +
  geom_point(size = 10, alpha = 0.8) +
  scale_x_continuous(name = "Age (days)") +
  scale_y_continuous(name = "") +
  scale_color_manual(labels = c("Nottingham", "Wirral"), name = "Site", values = c('#d95f02', '#7570b3')) +
  scale_shape_manual(values = c("circle", "triangle"), name = "Sex") +
  ggtitle("") +
  guides(color = guide_legend(order = 1, override.aes = list(size = 10, alpha = 1)), shape = guide_legend(order = 2, override.aes = list(size = 10, fill = "black", alpha = 1))) +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.key.height= unit(1.5, 'cm'),
    legend.key.width= unit(1.5, 'cm'),
    legend.position = "bottom", 
    legend.justification = "left",
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(8, 10.5))

#Save plots
ggsave(here("Results", "Updated_plots", "Final_plots", "Shannons_diversity_age_cd4a.png"), shannons.age.cd4a.plot, bg = "transparent", units = "cm", width = 30, height = 30)
ggsave(here("Results", "Updated_plots", "Final_plots", "Shannons_diversity_age_cd4b.png"), shannons.age.cd4b.plot2, bg = "transparent", units = "cm", width = 30, height = 30)
ggsave(here("Results", "Updated_plots", "Final_plots", "Shannons_diversity_age_cd8a.png"), shannons.age.cd8a.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave(here("Results", "Updated_plots", "Final_plots", "Shannons_diversity_age_cd8b.png"), shannons.age.cd8b.plot, bg = "transparent", units = "cm", width = 30, height = 20)

#patchwork plots
overall.shannon.plot <- (shannons.age.cd4a.plot + shannons.age.cd4b.plot) / (shannons.age.cd8a.plot + shannons.age.cd8b.plot) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 42, face = "bold"), legend.position = "bottom", legend.justification = "left")

#save patchwork
ggsave(here("Results", "Updated_plots", "Final_plots", "Shannons_plot.png"), overall.shannon.plot, bg = "transparent", units = "cm", width = 50, height = 50)

## Exploratory plots for Shannon's diversity -----------------------------------
message("Exploratory plots for Shannon's diversity")

# ##INCLUDING LAB MICE
# shannons.plot.lab <- ggplot(div.data.lab, aes(x = type, y = shannons, color = type)) +
#   geom_boxplot(linewidth = 2, outlier.size = 3) +
#   scale_x_discrete(name = "", labels = c("CD4 alpha", "CD4 beta", "CD8 alpha", "CD8 beta")) +
#   scale_y_continuous(name = "Shannon's") +
#   geom_bracket(xmin = c(1, 2, 3), xmax = c(2, 3, 4), label = "*", y.position = c(10.4, 10.5, 10.6), color = "black", inherit.aes = FALSE,
#                size = 1, label.size = 15, tip.length = 0.01, vjust = 0.5) +
#   ggtitle("") +
#   guides(color = "none") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(8, 10.6))
# 
# ggsave(here("Results", "Updated_plots", "Shannons_diversity_lab.png"), shannons.plot.lab, bg = "transparent", units = "cm", width = 30, height = 20)
# 
# #separate by cd4a/b and cd8a/b
# shannons.age.cd4a.plot.lab <- ggplot(div.data.lab[which(div.data.lab$type == "cd4a"), ], aes(x = age, y = shannons, color = site, shape = sex)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Age (days)") +
#   scale_y_continuous(name = "Shannon's") +
#   scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD4 alpha") +
#   guides(color = guide_legend(override.aes = list(size = 8)), shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(8, 10.5))
# 
# shannons.age.cd4b.plot.lab <- ggplot(div.data.lab[which(div.data.lab$type == "cd4b"), ], aes(x = age, y = shannons, color = site, shape = sex)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Age (days)") +
#   scale_y_continuous(name = "Shannon's") +
#   scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD4 beta") +
#   guides(color = guide_legend(override.aes = list(size = 8)), shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(8, 10.5))
# 
# shannons.age.cd8a.plot.lab <- ggplot(div.data.lab[which(div.data.lab$type == "cd8a"), ], aes(x = age, y = shannons, color = site, shape = sex)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Age (days)") +
#   scale_y_continuous(name = "Shannon's") +
#   scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD8 alpha") +
#   guides(color = guide_legend(override.aes = list(size = 8)), shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(8, 10.5))
# 
# shannons.age.cd8b.plot.lab <- ggplot(div.data.lab[which(div.data.lab$type == "cd8b"), ], aes(x = age, y = shannons, color = site, shape = sex)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Age (days)") +
#   scale_y_continuous(name = "Shannon's") +
#   scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
#   scale_shape_manual(values = c("circle", "square"), name = "Sex") +
#   ggtitle("CD8 beta") +
#   guides(color = guide_legend(override.aes = list(size = 8)), shape = guide_legend(order = 1, override.aes = list(size = 8, fill = "black"))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(8, 10.5))
# 
# #Save plots
# ggsave(here("Results", "Updated_plots", "Shannons_diversity_age_cd4a_lab.png"), shannons.age.cd4a.plot.lab, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "Shannons_diversity_age_cd4b_lab.png"), shannons.age.cd4b.plot.lab, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "Shannons_diversity_age_cd8a_lab.png"), shannons.age.cd8a.plot.lab, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "Shannons_diversity_age_cd8b_lab.png"), shannons.age.cd8b.plot.lab, bg = "transparent", units = "cm", width = 30, height = 20)


## Analysis for TCR sequences shared across mice -------------------------------
message("Analysis for TCR sequences shared across mice")


#function to find how many sequences are shared by X% of mice
sequences_shared <- function(df, n_mice, threshold){
  cutoff <- n_mice * threshold
  total <- sum(df$count >= cutoff)
  percentage <- (total / nrow(df)) * 100
  return(list(total = total, percentage = percentage))
}

#apply function to the data frames
sequences_shared(tcr.cd4a, n_mice = length(cd4.names), threshold = 0.75)
sequences_shared(tcr.cd4b, n_mice = length(cd4.names), threshold = 0.75)
sequences_shared(tcr.cd8a, n_mice = length(cd8.names), threshold = 0.75)
sequences_shared(tcr.cd8b, n_mice = length(cd8.names), threshold = 0.75)


## Plots for TCR sequences shared across mice ----------------------------------
message("Plots for TCR sequences shared across mice")

##CD4 alpha
tcr.cd4a.hist <- ggplot(tcr.cd4a, aes(x = count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "white") +
  scale_x_continuous(name = "Mice sharing a TCR sequence") +
  scale_y_log10(name = "Number of TCR sequences") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
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
  scale_x_continuous(name = "Mice sharing a TCR sequence") +
  scale_y_log10(name = "Number of TCR sequences") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
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
  scale_x_continuous(name = "Mice sharing a TCR sequence") +
  scale_y_log10(name = "Number of TCR sequences") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

##CD8 beta
tcr.cd8b.hist <- ggplot(tcr.cd8b, aes(x = count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "white") +
  scale_x_continuous(name = "Mice sharing a TCR sequence") +
  scale_y_log10(name = "Number of TCR sequences") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'))

#patchwork plots
overall.tcr.freq.plot <- (tcr.cd4a.hist + tcr.cd4b.hist) / (tcr.cd8a.hist + tcr.cd8b.hist) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 42, face = "bold"), legend.position = "bottom", legend.justification = "left",
        plot.margin = unit(c(5, 10, 5, 5), "mm"))

#save patchwork
ggsave(here("Results", "Updated_plots", "Final_plots", "TCR_frequencies.png"), overall.tcr.freq.plot, bg = "transparent", units = "cm", width = 55, height = 40)


## Plot shared TCR sequences ---------------------------------------------------
message("Plot shared TCR sequences")

#set colors
cols.cd4 <- rep('white', nrow(heat.cd4a))
#turn the specified rows the correct color
cols.cd4[grepl("^L", row.names(heat.cd4a))] <- '#1b9e77'
cols.cd4[grepl("^W", row.names(heat.cd4a))] <- '#7570b3'
cols.cd4[grepl("^N", row.names(heat.cd4a))] <- 'white'

#set colors
cols.cd8 <- rep('white', nrow(heat.cd8a))
#turn the specified rows the correct color
cols.cd8[grepl("^L", row.names(heat.cd8a))] <- '#1b9e77'
cols.cd8[grepl("^W", row.names(heat.cd8a))] <- '#7570b3'
cols.cd8[grepl("^N", row.names(heat.cd8a))] <- 'white'

png(here("Results", "Updated_plots", "Final_plots", "CD4_alpha_heatmap.png"),    #create PNG for the heat map        
    width = 6*300,        #5 x 300 pixels
    height = 6*300,
    res = 300,            #300 pixels per inch
    pointsize = 8)
heatmap.2(heat.cd4a, trace = "none", na.color = "white", scale = "none", 
          col = colorRampPalette(c("#ffffd4", "#fed98e", "#fe9929", "#d95f0e", "#993404")) (n = 20),
          breaks = seq(100, 6000, length.out = 21),
          keysize = 2, key.title = "", key.xlab = "", key.ylab = "", main = "", cexRow = 0.85, 
          colRow = cols.cd4, colCol = cols.cd4,
          offsetRow = -0.85, offsetCol = -0.85, ColSideColors = cols.cd4, RowSideColors = cols.cd4,
          labCol=as.expression(lapply(colnames(heat.cd4a), function(a) bquote(bold(.(a))))),
          labRow=as.expression(lapply(rownames(heat.cd4a), function(a) bquote(bold(.(a))))),
          key.xtickfun = function() {
            breaks = c(0, 3000, 6000)
            list(at = parent.frame()$scale01(breaks),
                 labels = breaks)
          })
dev.off()

png(here("Results", "Updated_plots", "Final_plots", "CD4_beta_heatmap.png"),    #create PNG for the heat map        
    width = 6*300,        #5 x 300 pixels
    height = 6*300,
    res = 300,            #300 pixels per inch
    pointsize = 8)
heatmap.2(heat.cd4b, trace = "none", na.color = "white", scale = "none", 
          col = colorRampPalette(c("#ffffd4", "#fed98e", "#fe9929", "#d95f0e", "#993404")) (n = 20),
          breaks = seq(100, 6000, length.out = 21),
          keysize = 2, key.title = "", key.xlab = "", key.ylab = "", main = "", cexRow = 0.85, 
          colRow = cols.cd4, colCol = cols.cd4,
          offsetRow = -0.85, offsetCol = -0.85, ColSideColors = cols.cd4, RowSideColors = cols.cd4,
          labCol=as.expression(lapply(colnames(heat.cd4a), function(a) bquote(bold(.(a))))),
          labRow=as.expression(lapply(rownames(heat.cd4a), function(a) bquote(bold(.(a))))),
          key.xtickfun = function() {
            breaks = c(0, 3000, 6000)
            list(at = parent.frame()$scale01(breaks),
                 labels = breaks)
          })
dev.off()

png(here("Results", "Updated_plots", "Final_plots", "CD8_alpha_heatmap.png"),    #create PNG for the heat map        
    width = 6*300,        #5 x 300 pixels
    height = 6*300,
    res = 300,            #300 pixels per inch
    pointsize = 8)
heatmap.2(heat.cd8a, trace = "none", na.color = "white", scale = "none", 
          col = colorRampPalette(c("#ffffd4", "#fed98e", "#fe9929", "#d95f0e", "#993404")) (n = 20),
          breaks = seq(100, 6000, length.out = 21),
          keysize = 2, key.title = "", key.xlab = "", key.ylab = "", main = "", cexRow = 0.85, 
          colRow = cols.cd8, colCol = cols.cd8,
          offsetRow = -0.85, offsetCol = -0.85, ColSideColors = cols.cd8, RowSideColors = cols.cd8,
          labCol=as.expression(lapply(colnames(heat.cd4a), function(a) bquote(bold(.(a))))),
          labRow=as.expression(lapply(rownames(heat.cd4a), function(a) bquote(bold(.(a))))),
          key.xtickfun = function() {
            breaks = c(0, 3000, 6000)
            list(at = parent.frame()$scale01(breaks),
                 labels = breaks)
          })
dev.off()

png(here("Results", "Updated_plots", "Final_plots", "CD8_beta_heatmap.png"),    #create PNG for the heat map        
    width = 6*300,        #5 x 300 pixels
    height = 6*300,
    res = 300,            #300 pixels per inch
    pointsize = 8)
heatmap.2(heat.cd8b, trace = "none", na.color = "white", scale = "none", 
          col = colorRampPalette(c("#ffffd4", "#fed98e", "#fe9929", "#d95f0e", "#993404")) (n = 20),
          breaks = seq(100, 6000, length.out = 21),
          keysize = 2, key.title = "", key.xlab = "", key.ylab = "", main = "", cexRow = 0.85, 
          key.par=list(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0)), colRow = cols.cd8, colCol = cols.cd8,
          offsetRow = -0.85, offsetCol = -0.85, ColSideColors = cols.cd8, RowSideColors = cols.cd8,
          labCol=as.expression(lapply(colnames(heat.cd4a), function(a) bquote(bold(.(a))))),
          labRow=as.expression(lapply(rownames(heat.cd4a), function(a) bquote(bold(.(a))))),
          key.xtickfun = function() {
            breaks = c(0, 3000, 6000)
            list(at = parent.frame()$scale01(breaks),
                 labels = breaks)
          })
dev.off()


## PCoA for shared TCR sequences -----------------------------------------------
message("PCoA for shared TCR sequences")

#calculate means and modes for sharing data
rep.cd4a.long2 <- rep.cd4a.long[which(!is.na(rep.cd4a.long$Freq)),]
rep.cd4b.long2 <- rep.cd4b.long[which(!is.na(rep.cd4b.long$Freq)),]
rep.cd8a.long2 <- rep.cd8a.long[which(!is.na(rep.cd8a.long$Freq)),]
rep.cd8b.long2 <- rep.cd8b.long[which(!is.na(rep.cd8b.long$Freq)),]

mean(rep.cd4a.long$Freq, na.rm = TRUE)
median(rep.cd4a.long$Freq, na.rm = TRUE)
mlv(rep.cd4a.long2$Freq, method = "mfv")

mean(rep.cd4b.long$Freq, na.rm = TRUE)
median(rep.cd4b.long$Freq, na.rm = TRUE)
mlv(rep.cd4b.long2$Freq, method = "mfv")

mean(rep.cd8a.long$Freq, na.rm = TRUE)
median(rep.cd8a.long$Freq, na.rm = TRUE)
mlv(rep.cd8a.long2$Freq, method = "mfv")

mean(rep.cd8b.long$Freq, na.rm = TRUE)
median(rep.cd8b.long$Freq, na.rm = TRUE)
mlv(rep.cd8b.long2$Freq, method = "mfv")

##rescale to be a dissimilarity matrix between 0 and 1
for(i in 1:62){
  for(j in 2:63){
    if(is.na(share.cd4a[i, j]) == TRUE){
      share.cd4a[i, j] <- NA
      share.cd4b[i, j] <- NA
    }
    else{
      share.cd4a[i, j] <- 1 - (share.cd4a[i, j]/5955)
      share.cd4b[i, j] <- 1 - (share.cd4b[i, j]/4461)
    }
    
  }
}

for(i in 1:58){
  for(j in 2:59){
    if(is.na(share.cd8a[i, j]) == TRUE){
      share.cd8a[i, j] <- NA
      share.cd8b[i, j] <- NA
    }
    else{
      share.cd8a[i, j] <- 1 - (share.cd8a[i, j]/5647)
      share.cd8b[i, j] <- 1 - (share.cd8b[i, j]/4858)
    }
    
  }
}

#make row names
rownames(share.cd4a) <- share.cd4a[, 1]
share.cd4a <- share.cd4a[, -1]

rownames(share.cd4b) <- share.cd4b[, 1]
share.cd4b <- share.cd4b[, -1]

rownames(share.cd8a) <- share.cd8a[, 1]
share.cd8a <- share.cd8a[, -1]

rownames(share.cd8b) <- share.cd8b[, 1]
share.cd8b <- share.cd8b[, -1]

#do pcoa
cd4a_pcoa <- wcmdscale(d = share.cd4a, eig = TRUE)
print(cd4a_pcoa)

cd4b_pcoa <- wcmdscale(d = share.cd4b, eig = TRUE)
print(cd4b_pcoa)

cd8a_pcoa <- wcmdscale(d = share.cd8a, eig = TRUE)
print(cd8a_pcoa)

cd8b_pcoa <- wcmdscale(d = share.cd8b, eig = TRUE)
print(cd8b_pcoa)

#get results ready for plotting
cd4a_pcoa2 <- data.frame(cd4a_pcoa$points)
cd4a_pcoa2$site <- c(rep("Lab", 6), rep("Nottingham", 52), rep("Wirral", 4))

cd4b_pcoa2 <- data.frame(cd4b_pcoa$points)
cd4b_pcoa2$site <- c(rep("Lab", 6), rep("Nottingham", 52), rep("Wirral", 4))

cd8a_pcoa2 <- data.frame(cd8a_pcoa$points)
cd8a_pcoa2$site <- c(rep("Lab", 6), rep("Nottingham", 48), rep("Wirral", 4))

cd8b_pcoa2 <- data.frame(cd8b_pcoa$points)
cd8b_pcoa2$site <- c(rep("Lab", 6), rep("Nottingham", 48), rep("Wirral", 4))

##add mouse age and sex
mouse.data.lab.cd4 <- mouse.data.lab[which(!is.na(mouse.data.lab$CD4.alpha.count)),]

cd4a_pcoa2$sex <- mouse.data.lab.cd4$Sex
cd4a_pcoa2$age <- mouse.data.lab.cd4$Age_days

cd4b_pcoa2$sex <- mouse.data.lab.cd4$Sex
cd4b_pcoa2$age <- mouse.data.lab.cd4$Age_days

mouse.data.lab.cd8 <- mouse.data.lab[which(!is.na(mouse.data.lab$CD8.alpha.count)),]

cd8a_pcoa2$sex <- mouse.data.lab.cd8$Sex
cd8a_pcoa2$age <- mouse.data.lab.cd8$Age_days

cd8b_pcoa2$sex <- mouse.data.lab.cd8$Sex
cd8b_pcoa2$age <- mouse.data.lab.cd8$Age_days


## Plot PCoA for shared TCR sequences ------------------------------------------
message("Plot PCoA for shared TCR sequences")

cd4a.pcoa.plot <- ggplot(data = cd4a_pcoa2, aes(x = Dim1, y = Dim2, color = site)) +
  geom_point(size = 6) +
  scale_x_continuous(name = "PCo1") +
  scale_y_continuous(name = "PCo2") +
  scale_color_manual(labels = c("Lab", "Nottingham", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("") +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 42, family = "Calibri", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

cd4b.pcoa.plot <- ggplot(data = cd4b_pcoa2, aes(x = Dim1, y = Dim2, color = site)) +
  geom_point(size = 6) +
  scale_x_continuous(name = "PCo1") +
  scale_y_continuous(name = "PCo2") +
  scale_color_manual(labels = c("Lab", "Nottingham", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("") +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 42, family = "Calibri", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

cd8a.pcoa.plot <- ggplot(data = cd8a_pcoa2, aes(x = Dim1, y = Dim2, color = site)) +
  geom_point(size = 6) +
  scale_x_continuous(name = "PCo1") +
  scale_y_continuous(name = "PCo2") +
  scale_color_manual(labels = c("Lab", "Nottingham", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("") +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 42, family = "Calibri", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

cd8b.pcoa.plot <- ggplot(data = cd8b_pcoa2, aes(x = Dim1, y = Dim2, color = site)) +
  geom_point(size = 6) +
  scale_x_continuous(name = "PCo1") +
  scale_y_continuous(name = "PCo2") +
  scale_color_manual(labels = c("Lab", "Nottingham", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
  ggtitle("") +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  theme_bw() +
  theme(plot.title = element_text(size = 42, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 42, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 42, family = "Calibri", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

#save plots
ggsave(here("Results", "Updated_plots", "Final_plots", "CD4_alpha_PCOA.png"), cd4a.pcoa.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave(here("Results", "Updated_plots", "Final_plots", "CD4_beta_PCOA.png"), cd4b.pcoa.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave(here("Results", "Updated_plots", "Final_plots", "CD8_alpha_PCOA.png"), cd8a.pcoa.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave(here("Results", "Updated_plots", "Final_plots", "CD8_beta_PCOA.png"), cd8b.pcoa.plot, bg = "transparent", units = "cm", width = 30, height = 20)


## Exploratory plots PCoA for shared TCR sequences -----------------------------
message("Exploratory plots PCoA for shared TCR sequences")

# amin <- 0.2
# 
# highcol <- "black"
# 
# lowcol.hex <- as.hexmode(round(col2rgb(highcol) * amin + 255 * (1 - amin)))
# lowcol <- paste0("#",   sep = "",
#                  paste(format(lowcol.hex, width = 2), collapse = ""))
# 
# cd4a.pcoa.plot2 <- ggplot(data = cd4a_pcoa2[which(!is.na(cd4a_pcoa2$age)),], aes(x = Dim1, y = Dim2, color = site, shape = sex, alpha = age)) +
#   geom_point(aes(color = site), size = 4) +
#   geom_point(aes(fill = age), alpha = 0) +
#   scale_fill_gradient(name = "Age (days)", high = highcol, low = lowcol, limits = c(17, 345)) +
#   scale_x_continuous(name = "PCo1") +
#   scale_y_continuous(name = "PCo2") +
#   scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
#   scale_shape_manual(labels = c("F", "M"), name = "Sex", values = c("circle", "square")) +
#   scale_alpha_continuous(name = "Age (days)", range = c(0.2, 0.9), limits = c(17, 345)) +
#   guides(color = guide_legend(order = 2, override.aes = list(size = 8)),
#          shape = guide_legend(order = 1, override.aes = list(size = 8)),
#          alpha = "none") +
#   ggtitle("CD4 alpha") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm'))
# 
# cd4b.pcoa.plot <- ggplot(data = cd4b_pcoa2[which(!is.na(cd4b_pcoa2$age)),], aes(x = Dim1, y = Dim2, color = site, shape = sex, alpha = age)) +
#   geom_point(aes(color = site), size = 4) +
#   geom_point(aes(fill = age), alpha = 0) +
#   scale_fill_gradient(name = "Age (days)", high = highcol, low = lowcol, limits = c(17, 345)) +
#   scale_x_continuous(name = "PCo1") +
#   scale_y_continuous(name = "PCo2") +
#   scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
#   scale_shape_manual(labels = c("F", "M"), name = "Sex", values = c("circle", "square")) +
#   scale_alpha_continuous(name = "Age (days)", range = c(0.2, 0.9), limits = c(17, 345)) +
#   guides(color = guide_legend(order = 2, override.aes = list(size = 8)),
#          shape = guide_legend(order = 1, override.aes = list(size = 8)),
#          alpha = "none") +
#   ggtitle("CD4 beta") +
#   guides(color = guide_legend(override.aes = list(size = 8))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm'))
# 
# cd8a.pcoa.plot <- ggplot(data = cd8a_pcoa2[which(!is.na(cd8a_pcoa2$age)),], aes(x = Dim1, y = Dim2, color = site, shape = sex, alpha = age)) +
#   geom_point(aes(color = site), size = 4) +
#   geom_point(aes(fill = age), alpha = 0) +
#   scale_fill_gradient(name = "Age (days)", high = highcol, low = lowcol, limits = c(17, 345)) +
#   scale_x_continuous(name = "PCo1") +
#   scale_y_continuous(name = "PCo2") +
#   scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
#   scale_shape_manual(labels = c("F", "M"), name = "Sex", values = c("circle", "square")) +
#   scale_alpha_continuous(name = "Age (days)", range = c(0.2, 0.9), , limits = c(17, 345)) +
#   guides(color = guide_legend(order = 2, override.aes = list(size = 8)),
#          shape = guide_legend(order = 1, override.aes = list(size = 8)),
#          alpha = "none") +
#   ggtitle("CD8 alpha") +
#   guides(color = guide_legend(override.aes = list(size = 8))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm'))
# 
# cd8b.pcoa.plot <- ggplot(data = cd8b_pcoa2[which(!is.na(cd8b_pcoa2$age)),], aes(x = Dim1, y = Dim2, color = site, shape = sex, alpha = age)) +
#   geom_point(aes(color = site), size = 4) +
#   geom_point(aes(fill = age), alpha = 0) +
#   scale_fill_gradient(name = "Age (days)", high = highcol, low = lowcol, limits = c(17, 345)) +
#   scale_x_continuous(name = "PCo1") +
#   scale_y_continuous(name = "PCo2") +
#   scale_color_manual(labels = c("Lab", "Nott.", "Wirral"), name = "", values = c('#1b9e77', '#d95f02', '#7570b3')) +
#   scale_shape_manual(labels = c("F", "M"), name = "Sex", values = c("circle", "square")) +
#   scale_alpha_continuous(name = "Age (days)", range = c(0.2, 0.9), limits = c(17, 345)) +
#   guides(color = guide_legend(order = 2, override.aes = list(size = 8)),
#          shape = guide_legend(order = 1, override.aes = list(size = 8)),
#          alpha = "none") +
#   ggtitle("CD8 beta") +
#   guides(color = guide_legend(override.aes = list(size = 8))) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 42, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm'))
# 
# #save plots
# ggsave(here("Results", "Updated_plots", "CD4_alpha_PCOA2.png"), cd4a.pcoa.plot2, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "CD4_beta_PCOA2.png"), cd4b.pcoa.plot, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "CD8_alpha_PCOA2.png"), cd8a.pcoa.plot, bg = "transparent", units = "cm", width = 30, height = 20)
# ggsave(here("Results", "Updated_plots", "CD8_beta_PCOA2.png"), cd8b.pcoa.plot, bg = "transparent", units = "cm", width = 30, height = 20)


## Data preparation for TCR sequence sharing analyses --------------------------
message("Data preparation for TCR sequence sharing analyses")

#get age and sex data
rep.cd4a.long <- rep.cd4a.long[which(!is.na(rep.cd4a.long$Freq)),]
rep.cd4a.long$mouse_1 <- mouse.data.lab.cd4$ID[rep.cd4a.long$Var1]
rep.cd4a.long$mouse2 <- mouse.data.lab.cd4$ID[rep.cd4a.long$Var2]
rep.cd4a.long$age_1 <- mouse.data.lab.cd4$Age_days[rep.cd4a.long$Var1]
rep.cd4a.long$age_2 <- mouse.data.lab.cd4$Age_days[rep.cd4a.long$Var2]
rep.cd4a.long$sex_1 <- mouse.data.lab.cd4$Sex[rep.cd4a.long$Var1]
rep.cd4a.long$sex_2 <- mouse.data.lab.cd4$Sex[rep.cd4a.long$Var2]
rep.cd4a.long$type <- "cd4a"

rep.cd4b.long <- rep.cd4b.long[which(!is.na(rep.cd4b.long$Freq)),]
rep.cd4b.long$mouse_1 <- mouse.data.lab.cd4$ID[rep.cd4b.long$Var1]
rep.cd4b.long$mouse2 <- mouse.data.lab.cd4$ID[rep.cd4b.long$Var2]
rep.cd4b.long$age_1 <- mouse.data.lab.cd4$Age_days[rep.cd4b.long$Var1]
rep.cd4b.long$age_2 <- mouse.data.lab.cd4$Age_days[rep.cd4b.long$Var2]
rep.cd4b.long$sex_1 <- mouse.data.lab.cd4$Sex[rep.cd4b.long$Var1]
rep.cd4b.long$sex_2 <- mouse.data.lab.cd4$Sex[rep.cd4b.long$Var2]
rep.cd4b.long$type <- "cd4b"

rep.cd8a.long <- rep.cd8a.long[which(!is.na(rep.cd8a.long$Freq)),]
rep.cd8a.long$mouse_1 <- mouse.data.lab.cd8$ID[rep.cd8a.long$Var1]
rep.cd8a.long$mouse2 <- mouse.data.lab.cd8$ID[rep.cd8a.long$Var2]
rep.cd8a.long$age_1 <- mouse.data.lab.cd8$Age_days[rep.cd8a.long$Var1]
rep.cd8a.long$age_2 <- mouse.data.lab.cd8$Age_days[rep.cd8a.long$Var2]
rep.cd8a.long$sex_1 <- mouse.data.lab.cd8$Sex[rep.cd8a.long$Var1]
rep.cd8a.long$sex_2 <- mouse.data.lab.cd8$Sex[rep.cd8a.long$Var2]
rep.cd8a.long$type <- "cd8a"

rep.cd8b.long <- rep.cd8b.long[which(!is.na(rep.cd8b.long$Freq)),]
rep.cd8b.long$mouse_1 <- mouse.data.lab.cd8$ID[rep.cd8b.long$Var1]
rep.cd8b.long$mouse2 <- mouse.data.lab.cd8$ID[rep.cd8b.long$Var2]
rep.cd8b.long$age_1 <- mouse.data.lab.cd8$Age_days[rep.cd8b.long$Var1]
rep.cd8b.long$age_2 <- mouse.data.lab.cd8$Age_days[rep.cd8b.long$Var2]
rep.cd8b.long$sex_1 <- mouse.data.lab.cd8$Sex[rep.cd8b.long$Var1]
rep.cd8b.long$sex_2 <- mouse.data.lab.cd8$Sex[rep.cd8b.long$Var2]
rep.cd8b.long$type <- "cd8b"

#get age differences, sex differences, age summation
rep.cd4a.long$age_dif <- abs(rep.cd4a.long$age_1 - rep.cd4a.long$age_2)
rep.cd4a.long$age_sum <- rep.cd4a.long$age_1 + rep.cd4a.long$age_2
rep.cd4a.long$sex_dif <- ifelse(rep.cd4a.long$sex_1 == rep.cd4a.long$sex_2, "same", "different")
rep.cd4a.long$sex_dif <- as.factor(rep.cd4a.long$sex_dif)
rep.cd4a.long.nott <- rep.cd4a.long[which(grepl("^N", rep.cd4a.long$mouse_1) & grepl("^N", rep.cd4a.long$mouse2)), ]

rep.cd4b.long$age_dif <- abs(rep.cd4b.long$age_1 - rep.cd4b.long$age_2)
rep.cd4b.long$age_sum <- rep.cd4b.long$age_1 + rep.cd4b.long$age_2
rep.cd4b.long$sex_dif <- ifelse(rep.cd4b.long$sex_1 == rep.cd4b.long$sex_2, "same", "different")
rep.cd4b.long$sex_dif <- as.factor(rep.cd4b.long$sex_dif)
rep.cd4b.long.nott <- rep.cd4b.long[which(grepl("^N", rep.cd4b.long$mouse_1) & grepl("^N", rep.cd4b.long$mouse2)), ]

rep.cd8a.long$age_dif <- abs(rep.cd8a.long$age_1 - rep.cd8a.long$age_2)
rep.cd8a.long$age_sum <- rep.cd8a.long$age_1 + rep.cd8a.long$age_2
rep.cd8a.long$sex_dif <- ifelse(rep.cd8a.long$sex_1 == rep.cd8a.long$sex_2, "same", "different")
rep.cd8a.long$sex_dif <- as.factor(rep.cd8a.long$sex_dif)
rep.cd8a.long.nott <- rep.cd8a.long[which(grepl("^N", rep.cd8a.long$mouse_1) & grepl("^N", rep.cd8a.long$mouse2)), ]

rep.cd8b.long$age_dif <- abs(rep.cd8b.long$age_1 - rep.cd8b.long$age_2)
rep.cd8b.long$age_sum <- rep.cd8b.long$age_1 + rep.cd8b.long$age_2
rep.cd8b.long$sex_dif <- ifelse(rep.cd8b.long$sex_1 == rep.cd8b.long$sex_2, "same", "different")
rep.cd8b.long$sex_dif <- as.factor(rep.cd8b.long$sex_dif)
rep.cd8b.long.nott <- rep.cd8b.long[which(grepl("^N", rep.cd8b.long$mouse_1) & grepl("^N", rep.cd8b.long$mouse2)), ]

#combine into one data frame
rep.long <- rbind(rep.cd4a.long, rep.cd4b.long, rep.cd8a.long, rep.cd8b.long)
rep.long$type <- as.factor(rep.long$type)

rep.long.nott <- rep.long[which(grepl("^N", rep.long$mouse_1) & grepl("^N", rep.long$mouse2)), ]


## Statistical tests for TCR sequence sharing analyses -------------------------
message("Statistical tests for TCR sequence sharing analyses")

share.mod <- glm(Freq ~ age_dif + sex_dif + age_sum + type + type*age_dif + type*age_sum, data = rep.long, family = "gaussian")
summary(share.mod)

#test fit
share.output <- simulateResiduals(fittedModel = share.mod)
plot(share.output)
#DHARMa says the fit is okay

#Post hoc test
share.post <- emmeans(share.mod, specs = pairwise ~ type, weights = "proportional")
summary(share.post)
#all significant

#get trends for age difference
emtrends(share.mod, pairwise ~ type, var = "age_dif", adjust = "tukey")

#get trends for age sum
emtrends(share.mod, pairwise ~ type, var = "age_sum", adjust = "tukey")

## Plots for TCR sequence sharing with age differences -------------------------
message("Plots for TCR sequence sharing with age differences")

# CD4.alpha.share.age.plot <- ggplot(data = rep.cd4a.long, aes(x = age_dif, y = Freq)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Difference in pairwise age (days)") +
#   scale_y_continuous(name = "No. shared sequences") +
#   ggtitle("CD4 alpha") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 28, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(0, 6000), xlim = c(0, 330))
# 
# CD4.beta.share.age.plot <- ggplot(data = rep.cd4b.long, aes(x = age_dif, y = Freq)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Difference in pairwise age (days)") +
#   scale_y_continuous(name = "No. shared sequences") +
#   ggtitle("CD4 beta") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 28, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(0, 6000), xlim = c(0, 330))
# 
# CD8.alpha.share.age.plot <- ggplot(data = rep.cd8a.long, aes(x = age_dif, y = Freq)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Difference in pairwise age (days)") +
#   scale_y_continuous(name = "No. shared sequences") +
#   ggtitle("CD8 alpha") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 28, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(0, 6000), xlim = c(0, 330))
# 
# CD8.beta.share.age.plot <- ggplot(data = rep.cd8b.long, aes(x = age_dif, y = Freq)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Difference in pairwise age (days)") +
#   scale_y_continuous(name = "No. shared sequences") +
#   ggtitle("CD8 beta") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 28, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(0, 6000), xlim = c(0, 330))
# 
# ##save plots
# ggsave(here("Results", "Updated_plots", "CD4_alpha_sharing_age.png"), CD4.alpha.share.age.plot, bg = "transparent", units = "cm", width = 30, height = 23)
# ggsave(here("Results", "Updated_plots", "CD4_beta_sharing_age.png"), CD4.beta.share.age.plot, bg = "transparent", units = "cm", width = 30, height = 23)
# ggsave(here("Results", "Updated_plots", "CD8_alpha_sharing_age.png"), CD8.alpha.share.age.plot, bg = "transparent", units = "cm", width = 30, height = 23)
# ggsave(here("Results", "Updated_plots", "CD8_beta_sharing_age.png"), CD8.beta.share.age.plot, bg = "transparent", units = "cm", width = 30, height = 23)


## Plots for TCR sequence sharing with age sums --------------------------------
message("Plots for TCR sequence sharing with age sums")
# 
# CD4.alpha.share.age.plot2 <- ggplot(data = rep.cd4a.long, aes(x = age_sum, y = Freq)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Sum of pairwise age (days)") +
#   scale_y_continuous(name = "No. shared sequences") +
#   ggtitle("CD4 alpha") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 28, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(0, 6000), xlim = c(0, 660))
# 
# CD4.beta.share.age.plot2 <- ggplot(data = rep.cd4b.long, aes(x = age_sum, y = Freq)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Sum of pairwise age (days)") +
#   scale_y_continuous(name = "No. shared sequences") +
#   ggtitle("CD4 beta") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 28, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(0, 6000), xlim = c(0, 660))
# 
# CD8.alpha.share.age.plot2 <- ggplot(data = rep.cd8a.long, aes(x = age_sum, y = Freq)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Sum of pairwise age (days)") +
#   scale_y_continuous(name = "No. shared sequences") +
#   ggtitle("CD8 alpha") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 28, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(0, 6000), xlim = c(0, 660))
# 
# CD8.beta.share.age.plot2 <- ggplot(data = rep.cd8b.long, aes(x = age_sum, y = Freq)) +
#   geom_point(size = 4) +
#   scale_x_continuous(name = "Sum of pairwise age (days)") +
#   scale_y_continuous(name = "No. shared sequences") +
#   ggtitle("CD8 beta") +
#   theme_bw() +
#   theme(plot.title = element_text(size = 28, family = "Calibri"),
#     axis.text = element_text(size = 38, color = "black", family = "Calibri"),
#     axis.title = element_text(size = 42, family = "Calibri"),
#     axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
#     axis.text.x = element_text(size = 38, color = "black"),
#     legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
#     legend.key.height= unit(1, 'cm'),
#     legend.key.width= unit(1, 'cm')) +
#   coord_cartesian(ylim = c(0, 6000), xlim = c(0, 660))
# 
# ##save plots
# ggsave(here("Results", "Updated_plots", "CD4_alpha_sharing_age2.png"), CD4.alpha.share.age.plot2, bg = "transparent", units = "cm", width = 30, height = 23)
# ggsave(here("Results", "Updated_plots", "CD4_beta_sharing_age2.png"), CD4.beta.share.age.plot2, bg = "transparent", units = "cm", width = 30, height = 23)
# ggsave(here("Results", "Updated_plots", "CD8_alpha_sharing_age2.png"), CD8.alpha.share.age.plot2, bg = "transparent", units = "cm", width = 30, height = 23)
# ggsave(here("Results", "Updated_plots", "CD8_beta_sharing_age2.png"), CD8.beta.share.age.plot2, bg = "transparent", units = "cm", width = 30, height = 23)


## Statistical tests for TCR sequence sharing analyses (Nottingham only) -------
message("Statistical tests for TCR sequence sharing analyses (Nottingham only)")

##overall
share.mod.nott <- glm(Freq ~ age_dif + sex_dif + age_sum + type + type*age_dif + type*age_sum, data = rep.long.nott, family = "gaussian")
summary(share.mod.nott)

#test fit
share.output.nott <- simulateResiduals(fittedModel = share.mod.nott)
plot(share.output.nott)
#DHARMa says the fit is okay

#Post hoc test
share.post.nott <- emmeans(share.mod.nott, specs = pairwise ~ type, weights = "proportional")
summary(share.post.nott)
#all significant

#get trends for age difference
emtrends(share.mod.nott, pairwise ~ type, var = "age_dif", adjust = "tukey")

#get trends for age sum
emtrends(share.mod.nott, pairwise ~ type, var = "age_sum", adjust = "tukey")


## Plots for TCR sequence sharing with age differences (Nottingham only) -------
message("Plots for TCR sequence sharing with age differences (Nottingham only)")

CD4.alpha.share.age.plot.nott <- ggplot(data = rep.cd4a.long.nott, aes(x = age_dif, y = Freq)) +
  geom_point(size = 4, color = "#d95f02") +
  scale_x_continuous(name = "Difference in pairwise age (days)") +
  scale_y_continuous(name = "Number shared") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 6000), xlim = c(0, 330))

CD4.beta.share.age.plot.nott <- ggplot(data = rep.cd4b.long.nott, aes(x = age_dif, y = Freq)) +
  geom_point(size = 4, color = "#d95f02") +
  scale_x_continuous(name = "Difference in pairwise age (days)") +
  scale_y_continuous(name = "Number shared") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 6000), xlim = c(0, 330))

CD8.alpha.share.age.plot.nott <- ggplot(data = rep.cd8a.long.nott, aes(x = age_dif, y = Freq)) +
  geom_point(size = 4, color = "#d95f02") +
  scale_x_continuous(name = "Difference in pairwise age (days)") +
  scale_y_continuous(name = "Number shared") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 6000), xlim = c(0, 330))

CD8.beta.share.age.plot.nott <- ggplot(data = rep.cd8b.long.nott, aes(x = age_dif, y = Freq)) +
  geom_point(size = 4, color = "#d95f02") +
  scale_x_continuous(name = "Difference in pairwise age (days)") +
  scale_y_continuous(name = "Number shared") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 6000), xlim = c(0, 330))

##save plots
ggsave(here("Results", "Updated_plots", "Final_plots", "CD4_alpha_sharing_age_nottingham.png"), CD4.alpha.share.age.plot.nott, bg = "transparent", units = "cm", width = 30, height = 23)
ggsave(here("Results", "Updated_plots", "Final_plots", "CD4_beta_sharing_age_nottingham.png"), CD4.beta.share.age.plot.nott, bg = "transparent", units = "cm", width = 30, height = 23)
ggsave(here("Results", "Updated_plots", "Final_plots", "CD8_alpha_sharing_age_nottingham.png"), CD8.alpha.share.age.plot.nott, bg = "transparent", units = "cm", width = 30, height = 23)
ggsave(here("Results", "Updated_plots", "Final_plots", "CD8_beta_sharing_age_nottingham.png"), CD8.beta.share.age.plot.nott, bg = "transparent", units = "cm", width = 30, height = 23)


## Plots for TCR sequence sharing with age sums (Nottingham only) --------------
message("Plots for TCR sequence sharing with age sums (Nottingham only)")

CD4.alpha.share.age.plot2.nott <- ggplot(data = rep.cd4a.long.nott, aes(x = age_sum, y = Freq)) +
  geom_point(size = 4, color = "#d95f02") +
  scale_x_continuous(name = "Sum of pairwise age (days)") +
  scale_y_continuous(name = "Number shared") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 6000), xlim = c(0, 660))

CD4.beta.share.age.plot2.nott <- ggplot(data = rep.cd4b.long.nott, aes(x = age_sum, y = Freq)) +
  geom_point(size = 4, color = "#d95f02") +
  scale_x_continuous(name = "Sum of pairwise age (days)") +
  scale_y_continuous(name = "Number shared") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 6000), xlim = c(0, 660))

CD8.alpha.share.age.plot2.nott <- ggplot(data = rep.cd8a.long.nott, aes(x = age_sum, y = Freq)) +
  geom_point(size = 4, color = "#d95f02") +
  scale_x_continuous(name = "Sum of pairwise age (days)") +
  scale_y_continuous(name = "Number shared") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 6000), xlim = c(0, 660))

CD8.beta.share.age.plot2.nott <- ggplot(data = rep.cd8b.long.nott, aes(x = age_sum, y = Freq)) +
  geom_point(size = 4, color = "#d95f02") +
  scale_x_continuous(name = "Sum of pairwise age (days)") +
  scale_y_continuous(name = "Number shared") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(size = 28, family = "Calibri"),
    axis.text = element_text(size = 38, color = "black", family = "Calibri"),
    axis.title = element_text(size = 42, family = "Calibri"),
    axis.title.y = element_text(size = 42, family = "Calibri", vjust = 1.5),
    axis.text.x = element_text(size = 38, color = "black"),
    legend.text = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.title = element_text(size = 38, family = "Calibri", face = "bold"),
    legend.key.height= unit(1, 'cm'),
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(0, 6000), xlim = c(0, 660))

##save plots
ggsave(here("Results", "Updated_plots", "Final_plots", "CD4_alpha_sharing_age2_nottingham.png"), CD4.alpha.share.age.plot2.nott, bg = "transparent", units = "cm", width = 30, height = 23)
ggsave(here("Results", "Updated_plots", "Final_plots", "CD4_beta_sharing_age2_nottingham.png"), CD4.beta.share.age.plot2.nott, bg = "transparent", units = "cm", width = 30, height = 23)
ggsave(here("Results", "Updated_plots", "Final_plots", "CD8_alpha_sharing_age2_nottingham.png"), CD8.alpha.share.age.plot2.nott, bg = "transparent", units = "cm", width = 30, height = 23)
ggsave(here("Results", "Updated_plots", "Final_plots", "CD8_beta_sharing_age2_nottingham.png"), CD8.beta.share.age.plot2.nott, bg = "transparent", units = "cm", width = 30, height = 23)

#patchwork for cd8a
overall.sharing.plot <- CD8.alpha.share.age.plot.nott / CD8.alpha.share.age.plot2.nott +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 42, face = "bold"))

cd4a.sharing.plot <- CD4.alpha.share.age.plot.nott + CD4.alpha.share.age.plot2.nott +
  plot_annotation(theme = theme(plot.tag = element_text(size = 42, face = "bold")))

cd4b.sharing.plot <- CD4.beta.share.age.plot.nott + CD4.beta.share.age.plot2.nott +
  plot_annotation(theme = theme(plot.tag = element_text(size = 42, face = "bold")))

cd8b.sharing.plot <- CD8.beta.share.age.plot.nott + CD8.beta.share.age.plot2.nott +
  plot_annotation(theme = theme(plot.tag = element_text(size = 42, face = "bold")))

#save patchwork
ggsave(here("Results", "Updated_plots", "Final_plots", "Sharing_plot.png"), overall.sharing.plot, bg = "transparent", units = "cm", width = 30, height = 46)
ggsave(here("Results", "Updated_plots", "Final_plots", "CD4_alpha_sharing_plot.png"), cd4a.sharing.plot, bg = "transparent", units = "cm", width = 60, height = 23)
ggsave(here("Results", "Updated_plots", "Final_plots", "CD4_beta_sharing_plot.png"), cd4b.sharing.plot, bg = "transparent", units = "cm", width = 60, height = 23)
ggsave(here("Results", "Updated_plots", "Final_plots", "CD8_beta_sharing_plot.png"), cd8b.sharing.plot, bg = "transparent", units = "cm", width = 60, height = 23)