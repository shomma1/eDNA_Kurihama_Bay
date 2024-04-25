
set.seed(999)

library(vegan)
library(tidyverse)

# Figure 5 =============================================================
# dbrda ================================================================

sampleinfo <- read.csv("./input/sample_info.csv")
df_org <- read.csv("./input/dataset.csv", row.names = 1)
df_pa <- replace(df_org, df_org > 0, 1)
df_relative <- as.data.frame(apply(df_org, 2, function(x) x/sum(x)))
df_relative[,is.na(colSums(df_relative))] <- rep(0, nrow(df_relative))

df <- df_relative

month_main <- c("(a) Nov", "(b) Jan")
month_list <- c("N", "J")
cols <- c("salmon", "royalblue")
dev.new()
par(mfrow = c(1,2))
for (i in 1:2) {
  month <- month_list[i]
  
  symbols <- sampleinfo[sampleinfo$Month == month &
                        sampleinfo$remove != T,  "Symbol"]
  envs <- sampleinfo[sampleinfo$Symbol %in% symbols, ]
  
  d <- df[rowSums(df) != 0, symbols]
  dist_eco <- vegdist(t(d), distance = "bray")
  hclust <- hclust(dist_eco, method = "average")
  group <- cutree(hclust, k=3)
  plot(hclust)
  rect.hclust(hclust, 3)
  
  res_0 <- dbrda(dist_eco ~ 1,
                 data = envs, distance="bray")
  res_full <- dbrda(dist_eco ~ PSU + WT + Chl_a +  Turb. + log10(qPCR_Carp),
                    data = envs, distance="bray")
  res_best <- ordistep(res_full, scope = formula(res_0), 
                       direction = "backward", permutations = 999)
  # ===== plot ====
  summary <- summary(res_best)
  rda1 <- round(summary$cont$importance["Proportion Explained", "dbRDA1"],3)
  rda2 <- round(summary$cont$importance["Proportion Explained", "dbRDA2"],3)
  
  plot(res_best, main = "", type ="p",
       xlab = paste0("dbRDA1 (", rda1*100, "%)"),
       ylab = paste0("dbRDA2 (", rda2*100, "%)"))
  points(summary$sites[,1], summary$sites[,2], 
         col = cols[factor(envs$Group3)], cex = 2,
         pch = c(16,17)[factor(envs$Tide)])
  title(month_main[i], line = 2, adj = 0)
  # plot(res_best, type = "t") # check
  
}

# legend ------------------------------------------------------------------

# plot(0, type = "n", ann = F, axes = F)
# legend("center", c("River High", "River Low", "Sea High", "Sea Low"),
#        pch = c(16,17,16,17),
#        col = c("salmon", "salmon", "royalblue", "royalblue"),box.lty = 0,
#        title = "Site")

