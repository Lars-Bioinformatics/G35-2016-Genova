library(tidyverse)

setwd("sduvarcall/G35-2016-Genova/MarkDuplicates/MRDetect_result/variants_final/mrdetect/")

patients = paste0("G35-P",1:10,"-")
for (patient in patients){
  files = list.files(pattern="_RESULT", recursive=T)
  data = map_dfr(.x = files, .f = function(file){read.csv(file, header = T, nrows = 1)})
  data = setNames(data,c(names(data)[2:6]))
  data = data[,1:5]
  case = data %>% filter(str_detect(row.names(data),patient))
  controls = data %>% filter(!str_detect(row.names(data),patient))
  z_score = case$sites.detected - mean(controls$sites.detected) / sd(controls$sites.detected)
  print(paste0(patient,": ",z_score, " - sites detected: ", case$sites.detected))
}
