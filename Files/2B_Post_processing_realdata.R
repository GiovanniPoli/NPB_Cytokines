##### Library and Functions         -----

library(ROCR)
library(devtools)
library(mcclust.ext)
library(BiocManager)
library(limma)
library(statmod)
library(coda)
library(kableExtra)
library(formattable)
library(knitr)
library(tidyverse)
library(rstudioapi)

directory_source = setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory_source = setwd(dirname(rstudioapi::getSourceEditorContext()$path))

setwd(directory_source)

source("Scripts/functions_v3.0.R")
source("Scripts/My_GG_Scripts.R")  

##### Data (required for log-lik.)  -----

cytokines_data <- read_csv("Data/Cytokine_2020.csv", 
                   col_types = cols(state = readr::col_factor(levels = c("healthy", "diseased")),
                                 location = readr::col_factor(levels = c("mucosa", "submucosa","serum")),
                                       id = readr::col_factor(levels = c("IBD1","IBD2","BD5","BD6","BD7","BD8",
                                                                         "BD9","BD10","BD11","BD12","BD20","BD21"))))

lower_limit = c(10,1246,160	,305 ,0.55,	12	,0.73,2.01,	2.26,	6.86,	3.27,2.33,	12	,6.57,	9.69,	20  ,17	,20	,8.84,
                7.42,9.23,	2.4	,7.57,1.89,	3.61,	5.1 ,	7.2)
upper_limit = c(41000,5104000,653700,312450,2250,12675,3000,	8250,9250,	28100,13400,9550,	49500,26900,
                39700,82500,	68500,	82000,	36200,30400,	37800,	9850,31000,	7750,	14800,5225,	29500)

cbind(lower_limit,upper_limit)

labels_variable = c("GM-CSF",
"sP-Selectin",
"ICAM-1",
"sE-Selectin",
"IFN-$\\alpha$",
"IFN-$\\gamma$",
"IL-1$\\alpha$",
"IL-1$\\beta$",
"IL-10",
"IL-12p70",
"IL-13",
"IL-17A",
"IL-18",
"IL-2",
"IL-21",
"IL-22",
"IL-23",
"IL-27",
"IL-4",
"IL-5",
"IL-6",
"IL-8",
"IL-9",
"IP-10",
"MCP-1",
"MIP-1$\\alpha$",
"TNF-$\\alpha$")
  
  
names(upper_limit) = names(lower_limit) = labels_variable

labels_location    = levels(cytokines_data$location)
labels_id          = levels(cytokines_data$id)
labels_state       = levels(cytokines_data$state)

cytokines_data = reshape2::melt(cytokines_data,  id.vars = c("id", "state", "location"))
cytokines_data = na.omit(cytokines_data)

cytokines_data$censored =   cytokines_data$value == 0

cytokines_data$value = log(cytokines_data$value)
cytokines_data$value[cytokines_data$censored] = log(lower_limit)[cytokines_data$variable][cytokines_data$censored]


##### CHAINS POST                   ----- 

CHAINS  = list.files("Chains_real_data/v3/")
chain.1 = readRDS(paste0("Chains_real_data/v3/", CHAINS[1]))
chain.2 = readRDS(paste0("Chains_real_data/v3/", CHAINS[2]))
chain.3 = readRDS(paste0("Chains_real_data/v3/", CHAINS[3]))
chain.4 = readRDS(paste0("Chains_real_data/v3/", CHAINS[4]))


inference_1 = post.inference(chain.1$chain.MCMC, plot = FALSE)
inference_2 = post.inference(chain.2$chain.MCMC, plot = FALSE)
inference_3 = post.inference(chain.3$chain.MCMC, plot = FALSE)
inference_4 = post.inference(chain.4$chain.MCMC, plot = FALSE)

chain.1 = Mass_point_group(chain.1$chain)
chain.2 = Mass_point_group(chain.2$chain)
chain.3 = Mass_point_group(chain.3$chain)
chain.4 = Mass_point_group(chain.4$chain)

Coocorence. = (Coocorence(chain.1) + Coocorence(chain.2) + Coocorence(chain.3) + Coocorence(chain.4) ) /4
colnames(Coocorence.) = rownames(Coocorence.) = labels_variable

cl_binder   = minbinder.ext(Coocorence.)
cl_VI       = minVI(Coocorence.)

g_Matrix = matrix(c(cl_binder$cl, 
                    cl_VI$cl), ncol = 27, byrow = TRUE)


##### FIGURE 3                      -----

## Mod v3 

rownames(g_Matrix) = c("BL","VIC")

p2_1   = My_gg.Heatmap_group(Coocorence., from = 0, to = 1, g_Matrix = g_Matrix) + ggtitle("Similarity matrix with no restrictions")
legend = get_legend(p2_1)
p2_1   = p2_1 + theme(legend.position = "none")

## Mod v4 0.5

CHAINS  = list.files("Chains_real_data/v4_1/")
chain.5 = readRDS(paste0("Chains_real_data/v4_1/", CHAINS[1]))
chain.6 = readRDS(paste0("Chains_real_data/v4_1/", CHAINS[2]))
chain.7 = readRDS(paste0("Chains_real_data/v4_1/", CHAINS[3]))
chain.8 = readRDS(paste0("Chains_real_data/v4_1/", CHAINS[4]))

chain.5 = Mass_point_group(chain.5$chain)
chain.6 = Mass_point_group(chain.6$chain)
chain.7 = Mass_point_group(chain.7$chain)
chain.8 = Mass_point_group(chain.8$chain)

Coocorence.2 = (Coocorence(chain.5) + 
                Coocorence(chain.6) +
                Coocorence(chain.7) + 
                Coocorence(chain.8) ) /4

colnames(Coocorence.2) = rownames(Coocorence.2) = labels_variable


cl_binder_2   = minbinder.ext(Coocorence.2)
cl_VI_2       = minVI(Coocorence.2)

g_Matrix_2 = matrix(c(cl_binder_2$cl, 
                       cl_VI_2$cl), ncol = 27, byrow = TRUE)

rownames(g_Matrix_2) = c("BL","VIC")

p2_2 = My_gg.Heatmap_group(Coocorence.2, from = 0, to = 1, g_Matrix = g_Matrix_2)  +
       ggtitle(TeX("$\\sigma^2_j \\approx 0.5$")) + theme(legend.position = "none")

CHAINS  = list.files("Chains_real_data/v4_2/")
chain.9  = readRDS(paste0("Chains_real_data/v4_2/", CHAINS[1]))
chain.10 = readRDS(paste0("Chains_real_data/v4_2/", CHAINS[2]))
chain.11 = readRDS(paste0("Chains_real_data/v4_2/", CHAINS[3]))
chain.12 = readRDS(paste0("Chains_real_data/v4_2/", CHAINS[4]))

chain.9  = Mass_point_group(chain.9$chain)
chain.10 = Mass_point_group(chain.10$chain)
chain.11 = Mass_point_group(chain.11$chain)
chain.12 = Mass_point_group(chain.12$chain)

Coocorence.3 = ( Coocorence(chain.9) + 
                 Coocorence(chain.10) +
                 Coocorence(chain.11) + 
                 Coocorence(chain.12) ) /4

colnames(Coocorence.3) = rownames(Coocorence.3) = labels_variable


cl_binder_3   = minbinder.ext(Coocorence.3)
cl_VI_3       = minVI(Coocorence.3)

g_Matrix_3 = matrix(c(cl_binder_3$cl, 
                      cl_VI_3$cl), ncol = 27, byrow = TRUE)

rownames(g_Matrix_3) = c("BL","VIC")

p2_3 = My_gg.Heatmap_group(Coocorence.3, from = 0, to = 1, g_Matrix = g_Matrix_3) +
       ggtitle(TeX("$\\sigma^2_j \\approx 0.1$")) + theme(legend.position = "none")

plot.export2 =  grid.arrange( p2_1, legend,
                              p2_2, p2_3, ncol = 2, nrow = 2)

ggsave("FIG2.PDF", grid.arrange( p2_1, legend,
                                 p2_2, p2_3, ncol = 2, nrow = 2),
       dpi = 800, width = 26, height = 26, units = "cm")




##### TABLE 3 and APPENDIX TABLE    -----

Ex_tab      = (inference_1$EX +
               inference_2$EX +
               inference_3$EX +
               inference_4$EX)/4

rownames(Ex_tab) = labels_variable
# C.3

C.3 = round(Ex_tab,4)

Ex_slab_tab = (inference_1$EX.Slab + 
               inference_2$EX.Slab + 
               inference_3$EX.Slab +
               inference_4$EX.Slab)/4
rownames(Ex_slab_tab) = labels_variable

# C.4
C.4 = round(Ex_slab_tab,4)



Ex_prob     = (inference_1$zero.probs +
               inference_2$zero.probs +
               inference_3$zero.probs + 
               inference_4$zero.probs)/4

TAB_3 = C.5 = Ex_prob
rownames(TAB_3) = labels_variable

# C.5 and TAB.3
 round(TAB_3,4)

#### log-likelihood chain 1
CHAINS  = list.files("Chains_real_data/v3/")
chain.1 = readRDS(paste0("Chains_real_data/v3/", CHAINS[1]))
chain.2 = readRDS(paste0("Chains_real_data/v3/", CHAINS[2]))
chain.3 = readRDS(paste0("Chains_real_data/v3/", CHAINS[3]))
chain.4 = readRDS(paste0("Chains_real_data/v3/", CHAINS[4]))

##### PLOT APPENDIX C               -----
dir.create("plot_Sup_C")   

gg_list = list()
  
labels_variable = levels(cytokines_data$variable)

for(cyto in labels_variable){
  for(loc in labels_location){
    tij_1 = sapply(chain.1$chain.MCMC, function(it) it$theta.star[[ it$rho[ cyto] ]][ loc])
    tij_2 = sapply(chain.2$chain.MCMC, function(it) it$theta.star[[ it$rho[ cyto] ]][ loc])
    tij_3 = sapply(chain.3$chain.MCMC, function(it) it$theta.star[[ it$rho[ cyto] ]][ loc])
    tij_4 = sapply(chain.4$chain.MCMC, function(it) it$theta.star[[ it$rho[ cyto] ]][ loc])
    
    gelman_th = coda::gelman.diag(mcmc.list( as.mcmc(tij_1),
                                             as.mcmc(tij_2), 
                                             as.mcmc(tij_3),
                                             as.mcmc(tij_4)))
    
    df.plot = data.frame(y   = c(tij_1,tij_2,tij_3,tij_4),
                         x   = rep(1:2500, 4),
                         col = factor(rep(1:4, each = 2500, levels = paste0("chain ",1:4))))
    
    a = ggplot(data = df.plot, aes(x = x, y = y, col = col)) + 
      geom_line() + theme_bw() +
      ylab(paste(cyto,"\n" ,loc)) + xlab(TeX(paste0("$\\rho_{Gelman}$=", round(gelman_th$psrf[1,],5)))) +
      theme(text = element_text( family = "serif"), legend.position = "none") 
    
    ggsave2(paste0("plot_Sup_C/",cyto,"_",loc,".pdf"),a)
    
    gg_list[[paste0(cyto,"/",loc)]] = a
  }
}
    
    
Cp1 = grid.arrange(grobs = gg_list[ 1 : 27], ncol = 3, norw = 9)
Cp2 = grid.arrange(grobs = gg_list[28 : 54], ncol = 3, norw = 9)
Cp3 = grid.arrange(grobs = gg_list[55 : 81], ncol = 3, norw = 9)

ggsave2("FIG_C1.pdf", Cp1,
        dpi = 800,
        height = 30, 
        width  = 20,
        units  = "cm")

ggsave2("FIG_C2.pdf", Cp2,
        dpi = 800,
        height = 30, 
        width  = 20,
        units  = "cm")

ggsave2("FIG_C3.pdf", Cp3,
        dpi = 800,
        height = 30, 
        width  = 20,
        units  = "cm")
