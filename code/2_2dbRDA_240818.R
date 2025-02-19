
set.seed(999)

library(vegan)
library(tidyverse)
library(dendextend)

# Figure 5 =============================================================

sampleinfo <- read.csv("./input/sample_info.csv")
df_org <- read.csv("./input/dataset.csv", row.names = 1)

df_pa <- replace(df_org, df_org > 0, 1)
df_relative <- as.data.frame(apply(df_org, 2, function(x) x/sum(x)))
df_relative[,is.na(colSums(df_relative))] <- rep(0, nrow(df_relative))

# case ================================================================

remove_OB <- T #or F

# dbrda ================================================================

# month_main <- c("(a) Nov", "(b) Jan")
month_list <- c("N", "J")
main <- c("Bray-Curtis", "", "Sorensen", "")
sub <- c("(a) November", "(b) January", "(c) November", "(d) January")
cols <- c("salmon", "royalblue")

dev.new()
device1 <- max(dev.list())
par(mfrow = c(2,2))

dev.new()
par(mfrow = c(2,2))
device2 <- max(dev.list())

dev.new()
par(mfrow = c(2,2))
device3 <- max(dev.list())

k <- 0
a<- NULL
for (j in 1:2) {
  
  ifelse(j == 1, df <- df_relative, df <- df_pa)
  ifelse(j == 1, distance <- "Bray-Curtis", distance <- "Sorensen")

  for (i in 1:2) {
    k <- k + 1
    
    month <- month_list[i]
    
    if (remove_OB == F) {
      symbols <- sampleinfo[sampleinfo$Month == month,  "Symbol"]
    } else if (remove_OB == T) {
      symbols <- sampleinfo[
        sampleinfo$Site2 != "S8" &
        sampleinfo$Month == month,  "Symbol"]
    }
    
    d <- df[rowSums(df) != 0, symbols]
    d <- d[, colSums(d) != 0]
    
    envs <- sampleinfo[sampleinfo$Symbol %in% colnames(d), ]
    envs$qPCR_Carp <- log10(envs$qPCR_Carp + 1)
    
    if (month == "N") {
      env <- select(envs, PSU, WT, -Chl_a, Turb., qPCR_Carp, Depth) # -Chl_a
    } else if (month == "J") {
      env <- select(envs, PSU, WT, Chl_a, Turb., qPCR_Carp, Depth)
    }
    
    env <- as.data.frame(scale(env))
    
    a <- append(a, print(colnames(d)))
    if (j == 1) {
      dist_eco <- vegdist(t(d), distance = "bray")
    } else if (j == 2) {
      dist_eco <- vegdist(t(d), distance = "bray")
      # dist_eco <- vegdist(t(d), distance = "jacca
    }
    
    dev.set(device1)
    hclust <- hclust(dist_eco, method = "average")

    dend <- as.dendrogram(hclust)
    labels(dend)
    labels_colors(dend) <- sapply(labels(dend), function (x) {
      ifelse( sampleinfo[sampleinfo$Symbol == x, "Group1"] == "R",
              "salmon", "royalblue")
    })
    dev.set(device1)
    plot(dend)
    # rect.dendrogram(dend, k = 0) ##
    title(main = main[k], adj = 0)
    mtext(sub[k], side = 3, line = 0, adj = 0)
        
    res_0 <- dbrda(dist_eco ~ 1,
                   data = env, distance="bray")
    res_full <- dbrda(dist_eco ~ .,
                      data = env, distance="bray")
    res_best <- ordistep(res_full, scope = formula(res_0), 
                         direction = "backward", permutations = 999)
    # = plot =
    summary <- summary(res_best)
    rda <- round(summary$cont$importance["Proportion Explained", 1:2],3)
    
    dev.set(device2)
    
    plot(res_best, main = "", display ="wa", type = "p",
         xlab = paste0(names(rda)[1], "(", rda[1]*100, "%)"),
         ylab = paste0(names(rda)[2], "(", rda[2]*100, "%)"))    
    points(summary$sites[,1], summary$sites[,2], 
           bg = cols[factor(envs$Group3)], cex = 2,
           pch = c(21,24)[factor(envs$Tide)])
    text(res_best, main = "", display ="bp", col = "black")  
    
    title(main = main[k], adj = 0)
    mtext(sub[k], side = 3, line = 0, adj = 0)
    
    
    # dev.set(device3)
    # plot(res_best, type = "t") # check
  }

}

## VIFによる説明変数の選択 =====================================================

sample <- sampleinfo[sampleinfo$Month == "J" ,  "Symbol"]
# sample <- sampleinfo[,  "Symbol"]


envs <- sampleinfo %>%
  filter(Symbol %in% sample) %>%
  select(Symbol, PSU, WT, Chl_a,  Turb., qPCR_Carp, Depth)

envs$qPCR_Carp <- log10( envs$qPCR_Carp + 1)

env <- select(envs, -Symbol)

vif <- round(diag(solve(cor(env))), 3)
vif
# OK


sample <- sampleinfo[sampleinfo$Month == "N" &
                       sampleinfo$Site2 != "S8",  "Symbol"]

envs <- sampleinfo %>%
  filter(Symbol %in% sample) %>%
  select(Symbol, PSU, WT, -Chl_a,  Turb., qPCR_Carp, Depth)

envs$qPCR_Carp <- log10( envs$qPCR_Carp + 1)

env <- select(envs, -Symbol)

vif <- round(diag(solve(cor(env))), 3)
vif

# end






