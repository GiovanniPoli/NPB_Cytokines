### 24/10/2023
##### Set directory  -----

library(rstudioapi)

directory_source = setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory_source = setwd(dirname(rstudioapi::getSourceEditorContext()$path))

setwd(directory_source)

Folder_cases = c("v3", 
  "v4_1", "v4_2", 
  "v3_equal", "v3_quarter","v3_doble","v3_diffuse",
  "hierarchical",
  "v3_M_1_1_4", "v3_M_10_2", "v3_M_50_10", 
  "v3_ns_1", "v3_ns_15", "v3_ns_30",         
  "v3_pi_5_15", "v3_pi_1_1", "v3_pi_15_5")       

dir.create("Chains_real_data")
for(s in Folder_cases){
  dir.create(paste("Chains_real_data",s,sep = "/"))
}

