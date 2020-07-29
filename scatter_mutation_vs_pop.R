.libPaths('~/R/x86_64-pc-linux-gnu-library/3.6')
library(tidyverse)
library(ggplot2)
library(data.table)
library(stringr)

mutant23 <- list.dirs("/hosts/linuxhome/mutant23/tmp/bramve/", full.names = TRUE, recursive = FALSE)
evolving_mutation <- grep('evolving', mutant23, value = TRUE)

popfile_em_vector <- c()
mutationfile_em_vector <- c()


timepointlist <- c()
cycle_points <- c("2500","5000","7500")
for(i in cycle_points){
  timepoints <- seq(paste0(100,as.numeric(i)),paste0(900,as.numeric(i)),by = 1000000)
  timepointlist <- c(timepointlist, timepoints)
}

for(i in evolving_mutation){
  for(cyclefraction in timepointlist){
    for(timepoints in cyclefraction){
      file_name <- (paste0(i, paste0("/grid/grid_type",as.character(timepoints),".txt")))
      popfile_em_vector <- c(popfile_em_vector, file_name)
    }
  }
}


for(i in evolving_mutation){
  for(cyclefraction in timepointlist){
    for(timepoints in cyclefraction){
      file_name <- (paste0(i, paste0("/grid/grid_mut",as.character(timepoints),".txt")))
      mutationfile_em_vector <- c(mutationfile_em_vector, file_name)
    }
  }
}

pop_name_vector <- c('time', 'seed','type', 'popsize_wt','popsize_mu', 'run')
mutation_name_vector <- c('time', 'seed', 'type', 'mutation','run')

popframe_em <- data.frame(time=numeric(0),seed=character(0),type=character(0),popsize_wt=numeric(0),popsize_mu=numeric(0),run=character(0))

for(i in popfile_em_vector){
  temp <- fread(i)
  colnames(temp) <- c('type', 'x1', 'y1')
  timepoint <- as.numeric(gsub('.*type|\\.txt$', '', head(i)))
  popsize_wt <- sum(temp$type == "1")
  popsize_mu <- sum(temp$type == "2")
  type <- as.character('evolving mutation')
  run <- as.character(str_match(i, 'bramve//(.*?)/grid'))
  seed <- as.character(str_match(run,'[0-9]+[^r]*$'))
  run_values <- data.frame(timepoint,seed[2],type,popsize_wt,popsize_mu,run[2])
  colnames(run_values) <- pop_name_vector
  popframe_em <- bind_rows(popframe_em,run_values)
}

popframe <- popframe_em

mutationframe_em <- data.frame(time=numeric(0), seed = character(0), type=character(0),mutation=numeric(0),run=character(0))

for(i in mutationfile_em_vector){
  temp <- fread(i)
  colnames(temp) <- c('mutation', 'x1', 'y1')
  timepoint <- as.numeric(gsub('.*mut|\\.txt$', '', head(i)))
  mutation <- mean(temp$mutation[temp$mutation!=0])
  type <- as.character('evolving mutation')
  run <- as.character(str_match(i, 'bramve//(.*?)/grid'))
  seed <- as.character(str_match(run,'[0-9]+[^r]*$'))
  run_values <- data.frame(timepoint,seed[2],type,mutation,run[2])
  colnames(run_values) <- mutation_name_vector
  mutationframe_em <- bind_rows(mutationframe_em,run_values)
}

mutationframe <- mutationframe_em

pop_mutation_frame <- merge(popframe,mutationframe)

long_pop_mutation_frame <- gather(pop_mutation_frame,key = 'wt_mu', value = 'popsize',c(popsize_wt,popsize_mu))

long_pop_mutation_frame %>%
  filter(time %in% seq(1002500, 9002500, by = 1000000)) %>%
  ggplot(aes(x = popsize, y = mutation, color = interaction(wt_mu,type))) +
  geom_point(shape = 16, size = 3) +
  geom_text(aes(label = seed, hjust = -1)) +
  geom_smooth(method = "lm") +
  ggtitle('mutation vs popsize quarter cycle') +
  theme_bw() +
  scale_colour_discrete(name = "Type", labels = c('Mutants evolving mutation', "Wildtypes evolving mutation")) +
  
  ggsave('~/Documents/strepto/scatter_popsize_mutation_em/scatter_2500_mutation_start_simulation.png',
         height = 210,
         width  = 297,
         units= 'mm',
         dpi= 150)

long_pop_mutation_frame %>%
  filter(time %in% seq(1005000, 9005000, by = 1000000)) %>%
  ggplot(aes(x = popsize, y = mutation, color = interaction(wt_mu,type))) +
  geom_point(shape = 16, size = 3) +
  geom_text(aes(label=seed, hjust = -1))+
  geom_smooth(method = "lm") +
  ggtitle('mutation vs popsize half cycle') +
  theme_bw() +
  scale_colour_discrete(name = "Type", labels = c('Mutants evolving mutation', "Wildtypes evolving mutation")) +
  
  ggsave('~/Documents/strepto/scatter_popsize_mutation_em/scatter_5000_mutation_start_simulation.png',
         height = 210,
         width  = 297,
         units= 'mm',
         dpi= 150)

long_pop_mutation_frame %>%
  filter(time %in% seq(1007500, 9007500, by = 1000000)) %>%
  ggplot(aes(x = popsize, y = mutation, color = interaction(wt_mu,type))) +
  geom_point(shape = 16, size = 3) +
  geom_text(aes(label = seed, hjust = -1)) +
  geom_smooth(method = "lm") +
  ggtitle('mutation vs popsize three-quarter cycle') +
  theme_bw() +
  scale_colour_discrete(name = "Type", labels = c('Mutants evolving mutation', "Wildtypes evolving mutation")) +
  
  ggsave('~/Documents/strepto/scatter_popsize_mutation_em/scatter_7500_mutation_start_simulation.png',
         height = 210,
         width  = 297,
         units= 'mm',
         dpi= 150)


