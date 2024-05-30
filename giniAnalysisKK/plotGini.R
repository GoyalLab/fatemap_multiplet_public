#set root directory where you will access raw files and output tables and plots
root <- "/Users/karunkiani/Desktop/doublet/"
setwd(root)

#point to the other script where the functions are kept
source(paste0(root, "extractionScripts/Gini_utils.R"))

#Load necessary librarires, set themes, set seed
loadLibs()
theme_set(theme_classic())
set.seed(123)

#set data and plot directories
dataDir <- here('data')
plotDir <- here('plots')

###read in Gini values from .tsv file created by extractGini.R

Gini_tib <- read_tsv(here('extractedData', "Gini_Plot_data.tsv"),
                     col_names = FALSE)

####Make Distribution Plot for Gini tibble####

p <- gghistogram(Gini_tib,  x = "Gini_Coeff",
                 add='mean', rug=TRUE,
                 color="resolution", fill ="resolution",
                 palette = pal[1:3], facet.by ="resolution",alpha = .9) +
  xlim(c(0,1))

ggsave(p, file = paste0(plotDir, "/giniDistributions.svg"),
       dpi = 300, height = 6, width = 12, units = "in")


### Create histograms of example distributions and their respective ####
## Gini distributiosn


#create empty distribution tibble
Dist_tib <- tibble(
  Gini_Coeff = numeric(), 
  sample = character(),
  distribution = character()
)  

#Function to find Gini distribution for multiple iterations of different dists
Gini_dist <- function(Dist_tib){
  for (i in 1:23){
    distribution <- "uniform"
    nums <- runif(n = 20, min = 0, max = 1)
    
    sample = paste0(distribution, "_", i)
    
    Dist_tib <- Dist_tib %>%
      add_row(Gini_Coeff = Gini(nums), 
              sample = paste0(distribution, "_", i), 
              distribution = distribution) 
  }
  
  for (i in 1:23){
    distribution <- "exponential"
    nums <- rexp(n = 20)
    
    sample = paste0(distribution, "_", i)
    
    Dist_tib <- Dist_tib %>%
      add_row(Gini_Coeff = Gini(nums), 
              sample = paste0(distribution, "_", i), 
              distribution = distribution) 
  }
  
  for (i in 1:23){
    distribution <- "power"
    nums <- poweRlaw::rplcon(20,xmin = 1, 1.5)
    
    sample = paste0(distribution, "_", i)
    
    Dist_tib <- Dist_tib %>%
      add_row(Gini_Coeff = Gini(nums), 
              sample = paste0(distribution, "_", i), 
              distribution = distribution) 
  }
  
  return(Dist_tib)
  
}

Dist_tib <- Gini_dist(Dist_tib)

Dist_tib <- Dist_tib %>% 
  mutate(distribution = fct_relevel(distribution,
                                    "uniform",
                                    "exponential",
                                    "power"))

##plot Gini Distributions of these simulated distributions
p2 <- gghistogram(Dist_tib,  x = "Gini_Coeff",
                  add='mean', rug=TRUE,
                  color="distribution", fill ="distribution",
                  facet.by ="distribution",alpha = .9) + xlim(c(0,1))

ggsave(p2, file = paste0(plotDir, "/exampleDistributions.svg"),
       dpi = 300, height = 6, width = 12, units = "in")

dat <- as.data.frame(runif(n = 10000, min = 0, max = 1))
colnames(dat) <- "val"

##plot histogram of the three example distrbutions used above for inset
p3 <- gghistogram(dat, x = "val",
                  add = "mean", rug = TRUE, 
                  color = "dodgerblue", fill = "dodgerblue",
                  alpha=0.9)

ggsave(p3, file = paste0(plotDir, "/uniform.svg"),
       dpi = 300)

dat <- as.data.frame(rexp(n =100000))
colnames(dat) <- "val"

p4 <- gghistogram(dat, x = "val",
                  add = "mean", rug = TRUE, 
                  color = "dodgerblue", fill = "dodgerblue",
                  alpha=0.9)

ggsave(p4, file = paste0(plotDir, "/exponential.svg"),
       dpi = 300)

dat <- as.data.frame(poweRlaw::rplcon(10000,xmin = 1, 1.5))
colnames(dat) <- "val"

p5 <- gghistogram(dat, x = "val",
                  add = "mean", rug = TRUE, 
                  color = "dodgerblue", fill = "dodgerblue",
                  alpha=0.9) + xlim(c(0, 100))
ggsave(p5, file = paste0(plotDir, "/power.svg"),
       dpi = 300)
