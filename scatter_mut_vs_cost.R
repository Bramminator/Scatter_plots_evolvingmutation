.libPaths('~/R/x86_64-pc-linux-gnu-library/3.6')
library(tidyverse)
library(ggplot2)
library(data.table)
library(stringr)

mutant23 <- list.dirs("/hosts/linuxhome/mutant23/tmp/bramve/", full.names = TRUE, recursive = FALSE)
evolving_mutation <- grep('evolving', mutant23, value = TRUE)

mutfile_em_vector <- c()
costfile_em_vector <- c()

timepointlist <- c()
cycle_points <- c("2500","5000","7500")
for(i in cycle_points){
  timepoints <- seq(paste0(4400,as.numeric(i)),paste0(5300,as.numeric(i)),by = 1000000)
  timepointlist <- c(timepointlist, timepoints)
}


for(i in evolving_mutation){
  for(cyclefraction in timepointlist){
    for(timepoints in cyclefraction){
      file_name <- (paste0(i, paste0("/grid/grid_mut",as.character(timepoints),".txt")))
      mutfile_em_vector <- c(mutfile_em_vector, file_name)
    }
  }
}


for(i in evolving_mutation){
  for(cyclefraction in timepointlist){
    for(timepoints in cyclefraction){
      file_name <- (paste0(i, paste0("/grid/grid_costs",as.character(timepoints),".txt")))
      costfile_em_vector <- c(costfile_em_vector, file_name)
    }
  }
}

mut_name_vector <- c('time', 'seed', 'mut_wt', 'run')
cost_name_vector <- c('time', 'seed','cost','run')

mutframe_em <- data.frame(time=numeric(0),seed=character(0),mut_wt=numeric(0),run=character(0))

for(i in mutfile_em_vector){
  temp <- fread(i)
  colnames(temp) <- c('cost', 'x1', 'y1')
  timepoint <- as.numeric(gsub('.*mut|\\.txt$', '', head(i)))
  mut <- mean(temp$cost[temp$cost!=0])
  run <- as.character(str_match(i, 'bramve//(.*?)/grid'))
  seed <- as.character(str_match(run,'[0-9]+[^r]*$'))
  run_values <- data.frame(timepoint,seed[2],mut,run[2])
  colnames(run_values) <- mut_name_vector
  mutframe_em <- bind_rows(mutframe_em,run_values)
}


mutframe <- mutframe_em

costframe_em <- data.frame(time=numeric(0), seed = character(0),cost=numeric(0),run=character(0))

for(i in costfile_em_vector){
  temp <- fread(i)
  colnames(temp) <- c('cost', 'x1', 'y1')
  timepoint <- as.numeric(gsub('.*costs|\\.txt$', '', head(i)))
  cost <- pmin(0.6,mean(temp$cost[temp$cost!=0]))
  run <- as.character(str_match(i, 'bramve//(.*?)/grid'))
  seed <- as.character(str_match(run,'[0-9]+[^r]*$'))
  run_values <- data.frame(timepoint,seed[2],cost,run[2])
  colnames(run_values) <- cost_name_vector
  costframe_em <- bind_rows(costframe_em,run_values)
}


costframe <- costframe_em

mut_cost_frame <- merge(mutframe,costframe)

print(mut_cost_frame)


mut_cost_frame %>%
  filter(time %in% seq(44002500, 53002500, by = 1000000)) %>%
  ggplot(aes(x = cost, y = mut_wt)) +
  geom_point(shape = 16, size = 3) +
  geom_text(aes(label = seed, hjust = -1)) +
  geom_smooth(method = "lm") +
  ggtitle('mutation-rate vs cost quarter cycle halfway simulation') +
  theme_bw() +
  ylab("mutation-rate") +
  
  ggsave('~/Documents/strepto/scatter_cost_mut_em_nem/scatter_mut_cost_2500halfway.png',
         height = 210,
         width  = 297,
         units= 'mm',
         dpi= 150)

mut_cost_frame %>%
  filter(time %in% seq(44005000, 53005000, by = 1000000)) %>%
  ggplot(aes(x = cost, y = mut_wt)) +
  geom_point(shape = 16, size = 3) +
  geom_text(aes(label=seed, hjust = -1))+
  geom_smooth(method = "lm") +
  ggtitle('mutation-rate vs cost half cycle halfway simulation') +
  theme_bw() +
  ylab("mutation-rate") +

  ggsave('~/Documents/strepto/scatter_cost_mut_em_nem/scatter_mut_cost_5000halfway.png',
         height = 210,
         width  = 297,
         units= 'mm',
         dpi= 150)

mut_cost_frame %>%
  filter(time %in% seq(44007500, 53007500, by = 1000000)) %>%
  ggplot(aes(x = cost, y = mut_wt)) +
  geom_point(shape = 16, size = 3) +
  geom_text(aes(label = seed, hjust = -1)) +
  geom_smooth(method = "lm") +
  ggtitle('mutation-rate vs cost three-quarter cycle halfway simulation') +
  theme_bw() +
  ylab("mutation-rate") +

  
  ggsave('~/Documents/strepto/scatter_cost_mut_em_nem/scatter_mut_cost_7500halfway.png',
         height = 210,
         width  = 297,
         units= 'mm',
         dpi= 150)


