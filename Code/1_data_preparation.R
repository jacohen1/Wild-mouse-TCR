# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Project and code information ------------------------------------------------
##
## Project: Wild Mouse TCR
##
## Purpose of script: Prepare data for analysis
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

packages <- c('here', 'ggplot2', 'immunarch', 'stringr', 'igraph', 'showtext',
              'purrr', 'dplyr', 'tidyr', 'vegan', 'poweRlaw')
lapply(packages, library, character.only = TRUE)

#additional code for showtext package
font_add(family = "Calibri", regular = "C:/Windows/Fonts/calibri.ttf")
showtext_auto()
showtext_opts(dpi = 300)

## Load data -------------------------------------------------------------------
message("Load data")

##read in mouse metadata file
mouse.data <- read.csv(here("Data", "additional_data", "Mouse_metadata_tcr_mice_Jacob.csv"), header = TRUE)

##read in all .tsv files
cd4a.files = list.files(path = here("Data", "tcrs_wild_mice", "CD4_alpha"), pattern = "_CD4_1_alpha\\.tsv", full.names = TRUE)
cd4b.files = list.files(path = here("Data", "tcrs_wild_mice", "CD4_beta"), pattern = "_CD4_1_beta\\.tsv", full.names = TRUE)
cd8a.files = list.files(path = here("Data", "tcrs_wild_mice", "CD8_alpha"), pattern = "_CD8_1_alpha\\.tsv", full.names = TRUE)
cd8b.files = list.files(path = here("Data", "tcrs_wild_mice", "CD8_beta"), pattern = "_CD8_1_beta\\.tsv", full.names = TRUE)

#function for reading in files
read_in <- function(file){
  read.table(file, sep = "\t", header = TRUE)
}

cd4a <- lapply(cd4a.files, read_in)
cd4b <- lapply(cd4b.files, read_in)
cd8a <- lapply(cd8a.files, read_in)
cd8b <- lapply(cd8b.files, read_in)


## Add necessary columns to data -----------------------------------------------
message("Add necessary columns to data")

##add column for concatenated v_call, j_call, junction_aa
cd4a <- cd4a %>%
  map(~ .x %>%
        unite(full_seq, 'v_call', 'j_call', 'junction_aa', remove = FALSE, sep=""))

cd4b <- cd4b %>%
  map(~ .x %>%
        unite(full_seq, 'v_call', 'j_call', 'junction_aa', remove = FALSE, sep=""))

cd8a <- cd8a %>%
  map(~ .x %>%
        unite(full_seq, 'v_call', 'j_call', 'junction_aa', remove = FALSE, sep=""))

cd8b <- cd8b %>%
  map(~ .x %>%
        unite(full_seq, 'v_call', 'j_call', 'junction_aa', remove = FALSE, sep=""))

#pull names from files to rename lists
cd4a.names <- str_extract(cd4a.files, "(?<=M_)(.*?)(?=_C)")
names(cd4a) <- cd4a.names

cd4b.names <- str_extract(cd4b.files, "(?<=M_)(.*?)(?=_C)")
names(cd4b) <- cd4b.names

cd8a.names <- str_extract(cd8a.files, "(?<=M_)(.*?)(?=_C)")
names(cd8a) <- cd8a.names

cd8b.names <- str_extract(cd8b.files, "(?<=M_)(.*?)(?=_C)")
names(cd8b) <- cd8b.names

##add mouse id column
cd4a <- mapply(cbind, cd4a, "ID" = names(cd4a), SIMPLIFY = F)
cd4b <- mapply(cbind, cd4b, "ID" = names(cd4b), SIMPLIFY = F)
cd8a <- mapply(cbind, cd8a, "ID" = names(cd8a), SIMPLIFY = F)
cd8b <- mapply(cbind, cd8b, "ID" = names(cd8b), SIMPLIFY = F)


## Get number of rows where the Productive column is False ---------------------
message("Get number of rows where the Productive column is False")

# #count number of productive = FALSE rows in each file
# n.rem <- function(df){
#   n_before <- nrow(df)
#   df_filtered <- df[df$productive == TRUE, ]
#   n_removed <- n_before - nrow(df_filtered)
#   return(n_removed)
# }
# 
# #apply to the data sets
# cd4a_removed <- sapply(cd4a, n.rem)
# cd4b_removed <- sapply(cd4b, n.rem)
# cd8a_removed <- sapply(cd8a, n.rem)
# cd8b_removed <- sapply(cd8b, n.rem)
# 
# #remove the mice that are eventually excluded
# to_remove_cd4 <- c("N63", "N42", "N24", "N01", "N27", "N11", "N12", "N09", "N35")
# to_remove_cd8 <- c("N17", "N42", "N22", "N01", "N21", "N11", "N09", "N27", "N12", "N51", "N10", "N35", "N36")
# 
# #find the positions of these names
# positions_cd4 <- which(cd4a.names %in% to_remove_cd4)
# positions_cd8 <- which(cd8a.names %in% to_remove_cd8)
# 
# #remove the necessary mice
# cd4a_removed_filtered <- cd4a_removed[-positions_cd4]
# cd4b_removed_filtered <- cd4b_removed[-positions_cd4]
# cd8a_removed_filtered <- cd8a_removed[-positions_cd8]
# cd8b_removed_filtered <- cd8b_removed[-positions_cd8]
# 
# #get means and medians for each
# mean(cd4a_removed_filtered)
# median(cd4a_removed_filtered)
# mean(cd4b_removed_filtered)
# median(cd4b_removed_filtered)
# mean(cd8a_removed_filtered)
# median(cd8a_removed_filtered)
# mean(cd8b_removed_filtered)
# median(cd8b_removed_filtered)

## Clean data ------------------------------------------------------------------
message("Clean data")

#remove all rows with FALSE in productive column
f.rem <- function(df){
  df_filtered <- df[df$productive == TRUE, ]
  return(df_filtered)
}

#apply to the data
cd4a <- lapply(cd4a, f.rem)
cd4b <- lapply(cd4b, f.rem)
cd8a <- lapply(cd8a, f.rem)
cd8b <- lapply(cd8b, f.rem)

## Extract repertoire sizes ----------------------------------------------------
message("Extract repertoire sizes")

#CD4
cd4.tcr.count <- c()

for(j in 1:length(cd4a.names)){
  CD4.alpha.count <- sum(cd4a[[j]]$duplicate_count)
  CD4.beta.count <- sum(cd4b[[j]]$duplicate_count)
  results <- cbind(cd4a.names[[j]], CD4.alpha.count, CD4.beta.count)
  cd4.tcr.count <- rbind(cd4.tcr.count, results)
}

cd4.tcr.count <- as.data.frame(cd4.tcr.count)
cd4.tcr.count$CD4.alpha.count <- as.numeric(cd4.tcr.count$CD4.alpha.count)
cd4.tcr.count$CD4.beta.count <- as.numeric(cd4.tcr.count$CD4.beta.count)

#CD8
cd8.tcr.count <- c()

for(j in 1:length(cd8a.names)){
  CD8.alpha.count <- sum(cd8a[[j]]$duplicate_count)
  CD8.beta.count <- sum(cd8b[[j]]$duplicate_count)
  results <- cbind(cd8a.names[[j]], CD8.alpha.count, CD8.beta.count)
  cd8.tcr.count <- rbind(cd8.tcr.count, results)
}

cd8.tcr.count <- as.data.frame(cd8.tcr.count)
cd8.tcr.count$CD8.alpha.count <- as.numeric(cd8.tcr.count$CD8.alpha.count)
cd8.tcr.count$CD8.beta.count <- as.numeric(cd8.tcr.count$CD8.beta.count)

##combine all results into mouse.data data frame
#check IDs are in the same order
cd4.tcr.count$V1 == mouse.data$ID
cd8.tcr.count$V1 == mouse.data$ID

#add repertoire sizes to mouse.data
mouse.data$CD4.alpha.count <- cd4.tcr.count$CD4.alpha.count
mouse.data$CD4.beta.count <- cd4.tcr.count$CD4.beta.count
mouse.data$CD8.alpha.count <- cd8.tcr.count$CD8.alpha.count
mouse.data$CD8.beta.count <- cd8.tcr.count$CD8.beta.count

## Get number of unique amino acid sequences per mouse -------------------------
message("Get number of unique amino acid sequences per mouse")

i <- 1:length(cd4a.files)

#write a function
unique.calc <- function(i, x){
  un <- length(unique(x[[i]]$full_seq))
  return(un)
}

cd4a.unique <- lapply(i, x = cd4a, unique.calc)
cd4b.unique <- lapply(i, x = cd4b, unique.calc)
cd8a.unique <- lapply(i, x = cd8a, unique.calc)
cd8b.unique <- lapply(i, x = cd8b, unique.calc)

#add results to the mouse.data data frame
mouse.data$CD4.alpha.unique.seq <- unlist(cd4a.unique)
mouse.data$CD4.beta.unique.seq <- unlist(cd4b.unique)
mouse.data$CD8.alpha.unique.seq <- unlist(cd8a.unique)
mouse.data$CD8.beta.unique.seq <- unlist(cd8b.unique)

## Get alpha values from distributions of clone sizes --------------------------
message("Get alpha values from distributions of clone sizes")

#write function
alpha.val.calc <- function(i, x){
  power.law <- fit_power_law(x[[i]]$duplicate_count, implementation = "plfit")
  return(power.law$alpha)
}

cd4a.dist <- lapply(i, x = cd4a, alpha.val.calc)
cd4b.dist <- lapply(i, x = cd4b, alpha.val.calc)
cd8a.dist <- lapply(i, x = cd8a, alpha.val.calc)
cd8b.dist <- lapply(i, x = cd8b, alpha.val.calc)

#add results to the mouse.data data frame
mouse.data$CD4.alpha.powerlaw <- unlist(cd4a.dist)
mouse.data$CD4.beta.powerlaw <- unlist(cd4b.dist)
mouse.data$CD8.alpha.powerlaw <- unlist(cd8a.dist)
mouse.data$CD8.beta.powerlaw <- unlist(cd8b.dist)

## Get Simpson's diversity measures for each mouse -----------------------------
message("Get Simpson's diversity measures for each mouse")

#function based on the formula, rather than a package
simpsons.calc <- function(s, i, x, y){
  #set seed to be able to replicate
  set.seed(seed = s)
  #first need to downsample the data set
  #make dataframe with rows based on duplicate count
  accordion <- uncount(x[[i]], x[[i]]$duplicate_count)
  #make duplicate_count column 1 for everything
  accordion$duplicate_count <- 1
  #then randomly sample this data frame
  if(nrow(accordion) >= y){
    row.sample <- sample(1:nrow(accordion), y, replace = FALSE)
    down.sample <- accordion[row.sample, ]
    #then compile back into duplicate counts
    sim <- aggregate(down.sample['duplicate_count'], by = down.sample['full_seq'], sum)
    #then calculate simpson's diversity index
    sim$v2 <- (sim$duplicate_count/sum(sim$duplicate_count))^2
    simpsons <- sum(sim$v2)
  }  else{
    (simpsons <- NA)
  }
  return(simpsons)
}

#apply function
cd4a.simp <- lapply(s = 45, i, x = cd4a, y = 54908, simpsons.calc)
cd4b.simp <- lapply(s = 45, i, x = cd4b, y = 58772, simpsons.calc)
cd8a.simp <- lapply(s = 45, i, x = cd8a, y = 66881, simpsons.calc)
cd8b.simp <- lapply(s = 45, i, x = cd8b, y = 64224, simpsons.calc)

#add results to the mouse.data data frame
mouse.data$CD4.alpha.simpsons <- unlist(cd4a.simp)
mouse.data$CD4.beta.simpsons <- unlist(cd4b.simp)
mouse.data$CD8.alpha.simpsons <- unlist(cd8a.simp)
mouse.data$CD8.beta.simpsons <- unlist(cd8b.simp)


## Get Shannon's diversity measures for each mouse -----------------------------
message("Get Shannon's diversity measures for each mouse")

#function based on the formula, rather than a package
shannons.calc <- function(s, i, x, y){
  #set seed to be able to replicate
  set.seed(seed = s)
  #first need to downsample the data set
  #make dataframe with rows based on duplicate count
  accordion <- uncount(x[[i]], x[[i]]$duplicate_count)
  #make duplicate_count column 1 for everything
  accordion$duplicate_count <- 1
  #then randomly sample this data frame
  if(nrow(accordion) >= y){
    row.sample <- sample(1:nrow(accordion), y, replace = FALSE)
    down.sample <- accordion[row.sample, ]
    #then compile back into duplicate counts
    shan <- aggregate(down.sample['duplicate_count'], by = down.sample['full_seq'], sum)
    #then calculate shannon's diversity index
    shan$v2 <- (shan$duplicate_count/sum(shan$duplicate_count))*(log(shan$duplicate_count/sum(shan$duplicate_count)))
    shannons <- -(sum(shan$v2))
  }  else{
    (shannons <- NA)
  }
  return(shannons)
}

#apply function
cd4a.shan <- lapply(s = 30, i, x = cd4a, y = 54908, shannons.calc)
cd4b.shan <- lapply(s = 30, i, x = cd4b, y = 58772, shannons.calc)
cd8a.shan <- lapply(s = 30, i, x = cd8a, y = 66881, shannons.calc)
cd8b.shan <- lapply(s = 30, i, x = cd8b, y = 64224, shannons.calc)

#add results to the mouse.data data frame
mouse.data$CD4.alpha.shannons <- unlist(cd4a.shan)
mouse.data$CD4.beta.shannons <- unlist(cd4b.shan)
mouse.data$CD8.alpha.shannons <- unlist(cd8a.shan)
mouse.data$CD8.beta.shannons <- unlist(cd8b.shan)

##write csv for mouse.data
write.csv(mouse.data, here("Data" , "Mouse_data.csv"), row.names = FALSE)


## Exact matches of sequences between mice -------------------------------------
message("Exact matches of sequences between mice")

#function for down-sampling
down.sample <- function(s, i, x, y){
  #set seed to be able to replicate
  set.seed(seed = s)
  #first need to downsample the data set
  #make dataframe with rows based on duplicate count
  accordion <- uncount(x[[i]], x[[i]]$duplicate_count)
  #make duplicate_count column 1 for everything
  accordion$duplicate_count <- 1
  #then randomly sample this data frame
  if(nrow(accordion) >= y){
    row.sample <- sample(1:nrow(accordion), y, replace = FALSE)
    down.sample <- accordion[row.sample, ]
    #then compile back into duplicate counts
    res <- aggregate(down.sample['duplicate_count'], by = down.sample['full_seq'], sum)
  }  else{
    res <- as.data.frame(x[[i]]$full_seq)
    colnames(res) <- c("full_seq")
  }
  return(res)
}

#apply function
cd4a.down.sample <- lapply(s = 94, i, x = cd4a, y = 54908, down.sample)
cd4b.down.sample <- lapply(s = 94, i, x = cd4b, y = 58772, down.sample)
cd8a.down.sample <- lapply(s = 94, i, x = cd8a, y = 49379, down.sample)
cd8b.down.sample <- lapply(s = 94, i, x = cd8b, y = 64224, down.sample)

#create blank matrix
rep.cd4a <- matrix(nrow = 71, ncol = 71)
rep.cd4b <- matrix(nrow = 71, ncol = 71)
rep.cd8a <- matrix(nrow = 71, ncol = 71)
rep.cd8b <- matrix(nrow = 71, ncol = 71)

#set list number for rows
for(l in 1:71){
  #then list number for columns
  for(m in 1:71){
    #given that a does not equal b (because then everything would match)
    if(l != m){
      #get a vector which states TRUE when there is a match and FALSE when there isn't
      tr.fa.cd4a <- unique(cd4a.down.sample[[l]]$full_seq) %in% unique(cd4a.down.sample[[m]]$full_seq)
      tr.fa.cd4b <- unique(cd4b.down.sample[[l]]$full_seq) %in% unique(cd4b.down.sample[[m]]$full_seq)
      tr.fa.cd8a <- unique(cd8a.down.sample[[l]]$full_seq) %in% unique(cd8a.down.sample[[m]]$full_seq)
      tr.fa.cd8b <- unique(cd8b.down.sample[[l]]$full_seq) %in% unique(cd8b.down.sample[[m]]$full_seq)
      
      #count the number of TRUEs in the vector
      matches.cd4a <- sum(tr.fa.cd4a)
      matches.cd4b <- sum(tr.fa.cd4b)
      matches.cd8a <- sum(tr.fa.cd8a)
      matches.cd8b <- sum(tr.fa.cd8b)
    } else {
      matches.cd4a <- NA
      matches.cd4b <- NA
      matches.cd8a <- NA
      matches.cd8b <- NA
    }
    #put total count into the final matrix at row l, column m
    rep.cd4a[l, m] <- matches.cd4a
    rep.cd4b[l, m] <- matches.cd4b
    rep.cd8a[l, m] <- matches.cd8a
    rep.cd8b[l, m] <- matches.cd8b
  }
}

#turn matrices into data frames
rep.cd4a <- as.data.frame(rep.cd4a)
rep.cd4b <- as.data.frame(rep.cd4b)
rep.cd8a <- as.data.frame(rep.cd8a)
rep.cd8b <- as.data.frame(rep.cd8b)

#add column and row names
colnames(rep.cd4a) <- cd4a.names
rownames(rep.cd4a) <- cd4a.names
colnames(rep.cd4b) <- cd4a.names
rownames(rep.cd4b) <- cd4a.names
colnames(rep.cd8a) <- cd4a.names
rownames(rep.cd8a) <- cd4a.names
colnames(rep.cd8b) <- cd4a.names
rownames(rep.cd8b) <- cd4a.names

#write to csv files 
write.csv(rep.cd4a, here("Data", "CD4_alpha_Repertoire.csv"), row.names = TRUE)
write.csv(rep.cd4b, here("Data", "CD4_beta_Repertoire.csv"), row.names = TRUE)
write.csv(rep.cd8a, here("Data", "CD8_alpha_Repertoire.csv"), row.names = TRUE)
write.csv(rep.cd8b, here("Data", "CD8_beta_Repertoire.csv"), row.names = TRUE)

## Number of times specific TCR sequences appear -------------------------------
message("Number of times specific TCR sequences appear")

##get list of only unique aa sequences
unique.calc2 <- function(i, x){
  un <- unique(x[[i]]$full_seq)
  un2 <- as.data.frame(cbind(un, x[[i]]$ID[1], 1))
  colnames(un2) <- c("full_seq", "ID", "count")
  return(un2)
}

#apply function
cd4a.unique2 <- lapply(i, x = cd4a, unique.calc2)
cd4b.unique2 <- lapply(i, x = cd4b, unique.calc2)
cd8a.unique2 <- lapply(i, x = cd8a, unique.calc2)
cd8b.unique2 <- lapply(i, x = cd8b, unique.calc2)

#first make large data frame of all of the individual data frames
cd4a.full <- bind_rows(cd4a.unique2, .id = "column_label")
cd4a.full$count <- 1
cd4b.full <- bind_rows(cd4b.unique2, .id = "column_label")
cd4b.full$count <- 1
cd8a.full <- bind_rows(cd8a.unique2, .id = "column_label")
cd8a.full$count <- 1
cd8b.full <- bind_rows(cd8b.unique2, .id = "column_label")
cd8b.full$count <- 1

#then remove all excluded mice
cd4a.full <- subset(cd4a.full, !(cd4a.full$ID %in% c("N63", "N42", "N24", "N01", "N27", "N11", "N12", "N09", "N35")))
cd4b.full <- subset(cd4b.full, !(cd4b.full$ID %in% c("N63", "N42", "N24", "N01", "N27", "N11", "N12", "N09", "N35")))
cd8a.full <- subset(cd8a.full, !(cd8a.full$ID %in% c("N17", "N42", "N22", "N01", "N21", "N11", "N09", "N27", "N12", "N51", "N10", "N35", "N36")))
cd8b.full <- subset(cd8b.full, !(cd8b.full$ID %in% c("N17", "N42", "N22", "N01", "N21", "N11", "N09", "N27", "N12", "N51", "N10", "N35", "N36")))

#then compile back into duplicate counts
tcr.freq.cd4a <- aggregate(cd4a.full['count'], by = cd4a.full['full_seq'], sum)
tcr.freq.cd4b <- aggregate(cd4b.full['count'], by = cd4b.full['full_seq'], sum)
tcr.freq.cd8a <- aggregate(cd8a.full['count'], by = cd8a.full['full_seq'], sum)
tcr.freq.cd8b <- aggregate(cd8b.full['count'], by = cd8b.full['full_seq'], sum)

#write to csv files
write.csv(tcr.freq.cd4a, here("Data", "CD4_alpha_TCR_freq.csv"), row.names = FALSE)
write.csv(tcr.freq.cd4b, here("Data", "CD4_beta_TCR_freq.csv"), row.names = FALSE)
write.csv(tcr.freq.cd8a, here("Data", "CD8_alpha_TCR_freq.csv"), row.names = FALSE)
write.csv(tcr.freq.cd8b, here("Data", "CD8_beta_TCR_freq.csv"), row.names = FALSE)


## Run Kolmogorov-Smirnov test to test for power law distribution --------------
message("Run Kolmogorov-Smirnov test to test for power law distribution")
# 
# #remove the mice that are eventually excluded
# to_remove_cd4 <- c("N63", "N42", "N24", "N01", "N27", "N11", "N12", "N09", "N35")
# to_remove_cd8 <- c("N17", "N42", "N22", "N01", "N21", "N11", "N09", "N27", "N12", "N51", "N10", "N35", "N36")
# 
# #filter the lists
# cd4a_filtered <- cd4a[setdiff(names(cd4a), to_remove_cd4)]
# cd4b_filtered <- cd4b[setdiff(names(cd4b), to_remove_cd4)]
# cd8a_filtered <- cd8a[setdiff(names(cd8a), to_remove_cd8)]
# cd8b_filtered <- cd8b[setdiff(names(cd8b), to_remove_cd8)]
# 
# #function for one-sample KS test against power law
# ks_powerlaw <- function(df){
#   data <- df$duplicate_count
#   
#   #fit discrete power law
#   pl <- displ$new(data)
#   est <- estimate_xmin(pl)
#   pl$setXmin(est)
#   pl$setPars(estimate_pars(pl))
#   
#   #bootstrap KS test for p-value
#   ks <- bootstrap_p(pl, no_of_sims = 500, threads = 8)
#   
#   #classify result based on p-value
#   fit_result <- case_when(ks$p > 0.1 ~ "Fits power law",
#                           ks$p < 0.05 ~ "Does not fit power law",
#                           TRUE ~ "Borderline")
#   #return results
#   tibble(fit_result = fit_result)
# }
# 
# #apply to the lists
# cd4a_ks_test <- lapply(cd4a_filtered, ks_powerlaw)
# cd4b_ks_test <- lapply(cd4b_filtered, ks_powerlaw)
# cd8a_ks_test <- lapply(cd8a_filtered, ks_powerlaw)
# cd8b_ks_test <- lapply(cd8b_filtered, ks_powerlaw)
# 
# #add the list names as dataset identifiers
# cd4a_ks_test_df <- bind_rows(cd4a_ks_test, .id = "dataset_name")
# cd4b_ks_test_df <- bind_rows(cd4b_ks_test, .id = "dataset_name")
# cd8a_ks_test_df <- bind_rows(cd8a_ks_test, .id = "dataset_name")
# cd8b_ks_test_df <- bind_rows(cd8b_ks_test, .id = "dataset_name")
# 
# #get summaries
# cd4a_ks_outcomes <- cd4a_ks_test_df %>%
#   count(fit_result) %>%
#   mutate(percent = 100 * n / sum(n))
# 
# cd4b_ks_outcomes <- cd4b_ks_test_df %>%
#   count(fit_result) %>%
#   mutate(percent = 100 * n / sum(n))
# 
# cd8a_ks_outcomes <- cd8a_ks_test_df %>%
#   count(fit_result) %>%
#   mutate(percent = 100 * n / sum(n))
# 
# cd8b_ks_outcomes <- cd8b_ks_test_df %>%
#   count(fit_result) %>%
#   mutate(percent = 100 * n / sum(n))


## Plot repertoire sizes -------------------------------------------------------
message("Plot repertoire sizes")

##CD4 results
cd4.tcr.plot <- ggplot(data = mouse.data, aes(x = CD4.alpha.count, y = CD4.beta.count, color = Site)) +
  geom_point(size = 2) +
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
  geom_point(size = 2) +
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
  geom_point(size = 2) +
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
  geom_point(size = 2) +
  scale_x_continuous(name = "Log, CD8 alpha repertoire size") +
  scale_y_continuous(name = "Log, CD8 beta repertoire size") +
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
ggsave(here("Results", "Exploratory_plots", "CD4_TCR_count.png"), cd4.tcr.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave(here("Results", "Exploratory_plots", "CD4_TCR_count_log.png"), cd4.tcr.plot.log, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave(here("Results", "Exploratory_plots", "CD8_TCR_count.png"), cd8.tcr.plot, bg = "transparent", units = "cm", width = 30, height = 20)
ggsave(here("Results", "Exploratory_plots", "CD8_TCR_count_log.png"), cd8.tcr.plot.log, bg = "transparent", units = "cm", width = 30, height = 20)

## Plot alpha coefficient against mouse age ------------------------------------
message("Plot alpha coefficient against mouse age")

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

## Plot histogram for two individual mice (N22, N16) ---------------------------
message("Plot histogram for two individual mice (N22, N16)")

cd4a.n22 <- ggplot(data = cd4a[[27]], aes(x = duplicate_count)) +
  geom_histogram(bins = 50, color = "black", fill = "white") +
  scale_y_log10(name = "Frequency (Log)") +
  scale_x_log10(name = "Duplicate count (Log)") +
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
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

cd4b.n22 <- ggplot(data = cd4b[[27]], aes(x = duplicate_count)) +
  geom_histogram(bins = 50, color = "black", fill = "white") +
  scale_y_log10(name = "Frequency (Log)") +
  scale_x_log10(name = "Duplicate count (Log)") +
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
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

cd8a.n16 <- ggplot(data = cd8a[[21]], aes(x = duplicate_count)) +
  geom_histogram(bins = 50, color = "black", fill = "white") +
  scale_y_log10(name = "Frequency (Log)") +
  scale_x_log10(name = "Duplicate count (Log)") +
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
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

cd8b.n16 <- ggplot(data = cd8b[[21]], aes(x = duplicate_count)) +
  geom_histogram(bins = 50, color = "black", fill = "white") +
  scale_y_log10(name = "Frequency (Log)") +
  scale_x_log10(name = "Duplicate count (Log)") +
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
    legend.key.width= unit(1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

#Save plots
ggsave(here("Results", "Updated_plots", "CD4_alpha_N22_duplicatecount_histogram.png"), cd4a.n22, bg = "transparent", units = "cm", width = 30, height = 23)
ggsave(here("Results", "Updated_plots", "CD4_beta_N22_duplicatecount_histogram.png"), cd4b.n22, bg = "transparent", units = "cm", width = 30, height = 23)
ggsave(here("Results", "Updated_plots", "CD8_alpha_N16_duplicatecount_histogram.png"), cd8a.n16, bg = "transparent", units = "cm", width = 30, height = 23)
ggsave(here("Results", "Updated_plots", "CD8_beta_N16_duplicatecount_histogram.png"), cd8b.n16, bg = "transparent", units = "cm", width = 30, height = 23)
