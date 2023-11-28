### 24/10/2023
##### Set directory  -----

library(rstudioapi)

directory_source = setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory_source = setwd(dirname(rstudioapi::getSourceEditorContext()$path))

setwd(directory_source)

##### Create Folders -----

for(s in c("Scenario_I","Scenario_II","Scenario_III","Scenario_IV")){
  dir.create(s)
  for(folder in c( "Mod_hierarchical","Mod_hierarchical_ExFDR","Mod_v3","Mod_v3_ExFDR")){
    dir.create(paste(s,folder,sep = "/"))
    for(subfolder in c("high","low","mid")){
      dir.create(paste(s,folder,subfolder,sep = "/"))
    }
  }
}