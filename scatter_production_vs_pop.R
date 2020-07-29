.libPaths('~/R/x86_64-pc-linux-gnu-library/3.6')
library(tidyverse)
library(ggplot2)
library(data.table)
library(stringr)

mutant23 <- list.dirs("/hosts/linuxhome/mutant23/tmp/bramve/", full.names = TRUE, recursive = FALSE)
evolving_mutation <- grep('evolving', mutant23, value = TRUE)

mutant7 <- list.dirs("/hosts/linuxhome/mutant7/tmp/bramve/", full.names = TRUE, recursive = FALSE)
evolving_production <- grep('new_seeding', mutant7, value = TRUE)

popfile_em_vector <- c()
popfile_nem_vector <- c()
productionfile_em_vector <- c()
productionfile_nem_vector <- c()

timepointlist <- c()
cycle_points <- c("2500","5000","7500")
for(i in cycle_points){
  timepoints <- seq(paste0(8600,as.numeric(i)),paste0(9500,as.numeric(i)),by = 1000000)
  timepointlist <- c(timepointlist, timepoints)
}

for(i in evolving_production){
  for(cyclefraction in timepointlist){
    for(timepoints in cyclefraction){
      file_name <- (paste0(i, paste0("/grid/grid_type",as.character(timepoints),".txt")))
      popfile_nem_vector <- c(popfile_nem_vector,file_name)
    }
  }
}

for(i in evolving_production){
  for(cyclefraction in timepointlist){
    for(timepoints in cyclefraction){
      file_name <- (paste0(i, paste0("/grid/grid_production",as.character(timepoints),".txt")))
      productionfile_nem_vector <- c(productionfile_nem_vector, file_name)
    }
  }
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
      file_name <- (paste0(i, paste0("/grid/grid_production",as.character(timepoints),".txt")))
      productionfile_em_vector <- c(productionfile_em_vector, file_name)
    }
  }
}

pop_name_vector <- c('time', 'seed','type', 'popsize_wt','popsize_mu', 'run')
production_name_vector <- c('time', 'seed', 'type', 'production','run')

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

popframe_nem <- data.frame(time=numeric(0),seed=character(0),type=character(0),popsize_wt=numeric(0),popsize_mu=numeric(0),run=character(0))

for(i in popfile_nem_vector){
  temp <- fread(i)
  colnames(temp) <- c('type', 'x1', 'y1')
  timepoint <- as.numeric(gsub('.*type|\\.txt$', '', head(i)))
  popsize_wt <- sum(temp$type == "1")
  popsize_mu <- sum(temp$type == "2")
  type <- as.character("no evolving mutation")
  run <- as.character(str_match(i, 'bramve//(.*?)/grid'))
  seed <- as.character(str_match(run,'[0-9]+[^r]*$'))
  run_values <- data.frame(timepoint,seed[2],type,popsize_wt,popsize_mu,run[2])
  colnames(run_values) <- pop_name_vector
  popframe_nem <- bind_rows(popframe_nem,run_values)
}

popframe <- bind_rows(popframe_em, popframe_nem)

productionframe_em <- data.frame(time=numeric(0), seed = character(0), type=character(0),production=numeric(0),run=character(0))

for(i in productionfile_em_vector){
  temp <- fread(i)
  colnames(temp) <- c('production', 'x1', 'y1')
  timepoint <- as.numeric(gsub('.*production|\\.txt$', '', head(i)))
  production <- mean(temp$production[temp$production!=0])
  type <- as.character('evolving mutation')
  run <- as.character(str_match(i, 'bramve//(.*?)/grid'))
  seed <- as.character(str_match(run,'[0-9]+[^r]*$'))
  run_values <- data.frame(timepoint,seed[2],type,production,run[2])
  colnames(run_values) <- production_name_vector
  productionframe_em <- bind_rows(productionframe_em,run_values)
}

productionframe_nem <- data.frame(time=numeric(0),seed=character(0), type=character(0),production=numeric(0),run=character(0))

for(i in productionfile_nem_vector){
  temp <- fread(i)
  colnames(temp) <- c('production', 'x1', 'y1')
  timepoint <- as.numeric(gsub('.*production|\\.txt$', '', head(i)))
  production <- mean(temp$production[temp$production!=0])
  type <- as.character('no evolving mutation')
  run <- as.character(str_match(i, 'bramve//(.*?)/grid'))
  seed <- as.character(str_match(run,'[0-9]+[^r]*$'))
  run_values <- data.frame(timepoint, seed[2], type,production,run[2])
  colnames(run_values) <- production_name_vector
  productionframe_nem <- bind_rows(productionframe_nem,run_values)
}

productionframe <- bind_rows(productionframe_em,productionframe_nem)

pop_production_frame <- merge(popframe,productionframe)

long_pop_production_frame <- gather(pop_production_frame,key = 'wt_or_mu', value = 'popsize',c(popsize_wt,popsize_mu))

long_pop_production_frame %>%
  filter(time %in% seq(86002500, 95002500, by = 1000000)) %>%
  ggplot(aes(x = popsize, y = production, color = interaction(wt_or_mu,type))) +
  geom_point(shape = 16, size = 3) +
  geom_text(aes(label = seed, hjust = -1)) +
  geom_smooth(method = "lm") +
  ggtitle('production vs popsize quarter cycle') +
  theme_bw() +
  scale_colour_discrete(name = "Type", labels = c('Mutants evolving mutation', "Wildtypes evolving mutation", "Mutants constant mutation", "Wildtypes constant mutation")) +
  
  
  ggsave('~/Documents/strepto/scatter_popsize_production_em_nem/scatter_2500_production.png',
         height = 210,
         width  = 297,
         units= 'mm',
         dpi= 150)

long_pop_production_frame %>%
  filter(time %in% seq(86005000, 95005000, by = 1000000)) %>%
  ggplot(aes(x = popsize, y = production, color = interaction(wt_or_mu,type))) +
  geom_point(shape = 16, size = 3) +
  geom_text(aes(label=seed, hjust = -1))+
  geom_smooth(method = "lm") +
  ggtitle('production vs popsize half cycle') +
  theme_bw() +
  scale_colour_discrete(name = "Type", labels = c('Mutants evolving mutation', "Wildtypes evolving mutation", "Mutants constant mutation", "Wildtypes constant mutation")) +
  
  ggsave('~/Documents/strepto/scatter_popsize_production_em_nem/scatter_5000_production.png',
         height = 210,
         width  = 297,
         units= 'mm',
         dpi= 150)

long_pop_production_frame %>%
  filter(time %in% seq(86007500, 95007500, by = 1000000)) %>%
  ggplot(aes(x = popsize, y = production, color = interaction(wt_or_mu,type))) +
  geom_point(shape = 16, size = 3) +
  geom_text(aes(label = seed, hjust = -1)) +
  geom_smooth(method = "lm") +
  ggtitle('production vs popsize three-quarter cycle') +
  theme_bw() +
  scale_colour_discrete(name = "Type", labels = c('Mutants evolving mutation', "Wildtypes evolving mutation", "Mutants constant mutation", "Wildtypes constant mutation")) +
  
  
  ggsave('~/Documents/strepto/scatter_popsize_production_em_nem/scatter_7500_production.png',
         height = 210,
         width  = 297,
         units= 'mm',
         dpi= 150)


