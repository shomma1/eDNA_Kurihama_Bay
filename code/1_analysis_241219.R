
library(vegan)
library(tidyverse)

# ==============================================================================
# 共通
sampleinfo <- read.csv("./input/sample_info.csv")

df_org <- read.csv("./input/dataset.csv", row.names = 1)
df_relative <- apply(df_org, 2, function(x) {x/sum(x)})
df_relative <- replace(df_relative, is.na(df_relative), 0)
df_pa <- replace(df_relative, df_relative > 0, 1)

habitat <- read.csv("./input/habitat.csv", row.names = 1)

#===============================================================================
# Figure 2 ---------------------------------------------------------------------
#===============================================================================

tides <- c("High", "Low")
cols <- c("salmon", "royalblue") 
pchs <- c(21, 24)
titles <- c("(a) November", "(b) January")
dev.new()
par(mfrow = c(1,2))
for (month in c("Nov", "Jan")) {
  for (i in 1:2) {
    
    tide <- tides[[i]]
    mon <- substr(month, 1, 1)
    sample <- sampleinfo[sampleinfo$Tide == tide & sampleinfo$Month == mon,]
    
    col <- ifelse(sample$Group3 == "R", col <- "salmon", col <- "royalblue")
    
    plot(sample$distance_from_upstream_cumulative, log10(sample$qPCR_Carp + 1),  
         bg = col, pch = pchs[i], cex = 2,
         xlim = c(0, 5500), ylim = c(0,6),
         xlab = "Stream distance from RA [m]",
         ylab = "Log 10 (eDNA concentration + 1) [copies / L]")
    axis(side=3, sample$distance_from_upstream_cumulative, sample$Site2)
    
    ifelse (month == "Nov", title <- titles[[1]], title <- titles[[2]])
    mtext(title, side = 3, line = 2, cex = 1.6, adj = 0)

    # linear regression in River ===
    y  <- sample %>% filter(Group1 == "R") %>%
      select(Site, qPCR_Carp, distance_from_upstream_cumulative)
    x <- y$distance_from_upstream_cumulative
    model <- lm(log10(y$qPCR_Carp + 1) ~ x)
    p_value <- summary(model)$coefficients["x","Pr(>|t|)"]
    print(paste("p_value", mon, tide, round(p_value, 3)))
    
    if (p_value < 0.05) { # significant
      abline(model, lty = 2)
      y_pred <- predict(model, data.frame(x))
      lines(x, y=y_pred, lwd =2)
      
    } else { # not significant -> null model was accepted
      model <- lm(log10(y$qPCR_Carp + 1) ~ 1) # null model
      abline(model, lty = 2)
      y_pred <- predict(model, data.frame(x))
      lines(x, y=y_pred, lwd = 2)
    }
    
    # plot
    y <- sample %>% filter(Group1 == "IB") %>%
      select(Site, qPCR_Carp, distance_from_upstream_cumulative)
    mean_y_bay <- mean(log10(y$qPCR_Carp +1))
    lines(c(-100, mean(y$distance_from_upstream_cumulative)), # S2 to S7
          c(mean_y_bay,mean_y_bay),
          lty = 3)
    mean_x_bay <- mean(y$distance_from_upstream_cumulative)
    abline(v=mean_x_bay,
           col="black", lty = 4)
    print(paste("mean_y_bay",mon, tide, round(mean_y_bay, 2)))
    pred <- predict(model, data.frame(x = mean_x_bay))
    print(paste("predicted", mon, tide, round(pred,2)))
    par(new = T)
  }
  par(new = F)
}
par(new = F)
# ==============================================================================
# legend

dev.new()
cols <- c("salmon", "salmon", "royalblue","royalblue")
pchs <- c(21, 24, 21, 24)
plot.new()
legend("center", 
       c("River High", "River Low", "Marine High", "Marine Low"),
       pch = pchs,col = "black", title = "Sample", pt.bg = cols)

#===============================================================================
# Figure 3 ---------------------------------------------------------------------
#===============================================================================
# Heatmap

# processing ===================================================================
ttl_reads <- data.frame(Species = rownames(df_org),
                        ttl_reads = rowSums(df_org))

# habitatの並び替え．
habitat_ordered <- ttl_reads %>%
  left_join(habitat, by = "Species") %>%
  #   filter(Species %in% rownames(df_org)) %>%
  arrange(habitat_class, ttl_reads)

df_with_fishtype <- df_relative[habitat_ordered$Species, ]
df_with_fishtype <- cbind(df_with_fishtype, habitat_ordered)

# df_with_fishtype <- df_relative %>%
#   data.frame(Species = rownames(.)) %>%
#   left_join(habitat, by = "Species") %>%
#   select(-X)

# summarize
df_aggre <- df_with_fishtype %>%
  group_by(habitat_class) %>%
  summarise(across(-Species, sum))

df_aggre_sp <- df_with_fishtype %>%
  group_by(habitat_class) %>%
  summarise(across(-Species, function(x) {sum(x>0)} ))

# order
habitat <- ttl_reads %>%
  left_join(habitat, by = "Species") %>%
  arrange(habitat_class, ttl_reads)
g.list <- list()
d_long_2 <- data.frame()
for (i in 1:2) {
  
  month <- ifelse(i == 1, month <- "N", month <- "J")
  month_plot <- ifelse(i == 1, "November", "January")
  sample <- sampleinfo[sampleinfo$Month == month, "Symbol"]
  sample_2 <- gsub("(_J|_N)$", "", sample)

  d <- df_org[sample]
  d <- d[habitat$Species, ]
  
  d_long <- d %>%
    mutate(Species = rownames(.)) %>%
    pivot_longer(
      cols = -Species,
      names_to = "Sample",
      values_to = "n_reads") %>%
    mutate(color_group = case_when(
      n_reads == 0 ~ "0 ~ 9",
      n_reads > 0 & n_reads < 100 ~ "10 ~ 99",
      n_reads >= 100 & n_reads < 1000 ~ "100 ~ 999",
      n_reads >= 1000 & n_reads < 10000 ~ "1000 ~ 9999",
      n_reads >= 10000 ~ "10000 ~")) %>%
    mutate(
      Sample_2 = sample_2[match(Sample, sample)]
    )
  
  d_long <- d_long %>%
    left_join(habitat, by = "Species") %>%
    mutate(Species = factor(Species, levels = habitat$Species) )
  
  d_long_2 <- d_long %>%
    cbind(month = factor(month_plot, levels = c("November", "January"))) %>%
    rbind(d_long_2)

  # dev.new()
  g.list[[i]] <- ggplot(d_long, aes(x = Sample, y = Species, fill = color_group)) +
    geom_tile() +
    scale_fill_manual(
      values = c("0 ~ 9" = "white", 
                 "10 ~ 99" = "yellow",
                 "100 ~ 999" = "orange",
                 "1000 ~ 9999" = "red",
                 "10000 ~" = "purple")) +
    theme_minimal() +
    labs(title = "Heatmap of Species by Sample", x = "Sample", y = "Species",
         fill = "Relative reads") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5), 
      axis.text.y = element_text(size = 3),
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides(fill = guide_legend(
      override.aes = list(color = "black")
    ))
  dev.new()
  print(g.list[[i]])
}

# =====
d_long_2 %>%
  ggplot(aes(x = Sample, y = Species, fill = color_group, group = month)) +
  geom_tile() +
  scale_fill_manual(
    values = c("0 ~ 9" = "white", 
               "10 ~ 99" = "yellow",
               "100 ~ 999" = "orange",
               "1000 ~ 9999" = "red",
               "10000 ~" = "purple")) +
  theme_minimal() +
  # scale_x_discrete(
  #   labels = c( unique(d_long_2$Sample) ) # カスタムラベル
  # ) +
  labs(title = "(a) Metabarcoding results", x = "Sample", y = "Species",
       fill = "Relative reads") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5), 
    axis.text.y = element_text(size = 3),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  guides(fill = guide_legend(
    override.aes = list(color = "black")
  )) +
  facet_wrap(~month, scales = "free_x")



# ==============================================================================

## arrange order
sampleinfo <- rbind(sampleinfo[sampleinfo$Group1 == "R",],
                    sampleinfo[sampleinfo$Group1 == "IB",],
                    sampleinfo[sampleinfo$Group1 == "OB",])
sites <- unique(sampleinfo$Site)
symbol_order <- vector()
for (i in 1: length(sites)) {
  symbol_order <- 
    append(symbol_order, sampleinfo[sampleinfo$Site %in% sites[i], "Symbol"])
}
sampleinfo <- sampleinfo[order(sampleinfo$Symbol, symbol_order), ]

removed  <- sampleinfo[sampleinfo$remove == T, "Symbol"]

##
sampleinfo <- sampleinfo %>%
  mutate(Tide_1 = substring(Tide,1,1)) %>%
  mutate(Loc = Site) %>%
  mutate(Site = paste0(Site, Tide_1))
##

dev.new()
par(mfrow = c(4,2), mar = c(5.1,4.1,4.1,2.1))
# F, F-B, B, F-B-M, B-M, M
cols = c("red", "salmon","green3", "lightgreen","skyblue","royalblue")

# Number of Species ===
# November
sample <- sampleinfo[sampleinfo$Month == "N" , 
                     c("Symbol", "Site", "Loc", "Tide_1")]
plt <- df_aggre_sp[, sample$Symbol]
bp <- barplot(as.matrix(plt), xaxt = "n",
              las = 2, ylim = c(0, 70), col = cols)
title("(b) Number of species", adj = 0)
# January
sample <- sampleinfo[sampleinfo$Month == "J" , 
                     c("Symbol", "Site", "Loc", "Tide_1")]
plt <- df_aggre_sp[, sample$Symbol]
barplot(as.matrix(plt), xaxt = "n",
        las = 2, ylim = c(0, 70), col = cols)

# Effective number of species ===
# November
sample <- sampleinfo[sampleinfo$Month == "N", c("Symbol", "Site")]
plt <- sampleinfo[match(sample$Symbol, sampleinfo$Symbol), 
                  c("Symbol", "en_sp", "Loc", "Tide_1")]
bp <- barplot(plt$en_sp, xaxt = "n", 
              las = 2, ylim = c(0, 14))
title("(c) Effective number of species", adj = 0)

# January
sample <- sampleinfo[sampleinfo$Month == "J", c("Symbol", "Site")]
plt <- sampleinfo[match(sample$Symbol, sampleinfo$Symbol), 
                  c("Symbol", "en_sp", "Loc", "Tide_1")]
barplot(plt$en_sp,  xaxt = "n", las = 2, ylim = c(0, 14))

# relative read counts ===
# November
sample <- sampleinfo[sampleinfo$Month == "N" , 
                     c("Symbol", "Site", "Loc", "Tide_1")]
plt <- df_aggre[, sample$Symbol]
bp <- barplot(as.matrix(plt),  xaxt = "n", las = 2, col = cols)
title(main = "(d) Relative reads", adj = 0)

# January
sample <- sampleinfo[sampleinfo$Month == "J" ,
                     c("Symbol", "Site", "Loc", "Tide_1")]
plt <- df_aggre[, sample$Symbol]
barplot(as.matrix(plt),  xaxt = "n", las = 2, col = cols)

# total reads ===
# November
sample <- sampleinfo[sampleinfo$Month == "N" , 
                     c("Symbol", "Site", "Loc", "Tide_1")]
plt <- df_org[, sample$Symbol]
bp <- barplot(colSums(plt), las = 2, ylim = c(0, 120000))
title(main = "(e) Total reads", adj = 0)
# January
sample <- sampleinfo[sampleinfo$Month == "J" ,
                     c("Symbol", "Site", "Loc", "Tide_1")]
plt <- df_org[, sample$Symbol]
barplot(colSums(plt), las = 2, ylim = c(0, 120000))

# legend
dev.new()
plot.new()
legend("center", rev(c("F", "F-B","F-B-M", "B", "B-M", "M")), fill = rev(cols))

#==================================================================
# Figure 4 --------------------------------------------------------
#==================================================================
# dominant
sampleinfo <- read.csv("./input/sample_info.csv")

sampleinfo <- sampleinfo %>%
  mutate(Tide_1 = substring(Tide,1,1)) %>%
  mutate(Loc = Site) %>%
  mutate(Site = paste0(Site, Tide_1))
plot_bar <- function(month, sp, title, text) {
  sample <- sampleinfo[sampleinfo$Month == month, c("Symbol")]
  plt <- sampleinfo[match(sample, sampleinfo$Symbol), 
                    c("Symbol", sp, "Site", "Tide_1", "Loc")]
  bp <- barplot(t(cbind(plt[sp], 1 - plt[sp])), names.arg = plt$Symbol,
                las = 2,
                col = c("salmon", "white"),
                ylab = "Relative reads")
  title(bquote(italic(.(title))), line = 2, adj = 0, cex.main = 1.5)
  mtext(text, side = 3, line = 0.5, cex = 1.0, adj = 0)
  
  # detect <- plt[sp] != 0
  # points(bp[detect], rep(0.9, sum(detect)), pch = 8, col = "salmon")
  }

dev.new()
par(mfrow = c(2,2))
# month, sp, title, text
plot_bar("N", "Cyprinus.carpio",       "C. carpio",     "(a) November")
plot_bar("J", "Cyprinus.carpio",       "C. carpio",     "(b) January")
plot_bar("N", "Mugil.cephalus",        "M. cephalus",   "(c) November")
plot_bar("J", "Mugil.cephalus",        "M. cephalus",   "(d) January")
dev.new()
par(mfrow = c(1,2))
plot_bar("N", "Girella.punctata",      "G. punctata",   "(c) November")
plot_bar("J", "Lateolabrax.japonicus", "L. japonicus",  "(d) January")

#==================================================================
# Table 1 --------------------------------------------------------
#==================================================================

# for table1
df_pa <- as.data.frame(df_pa)

# aggre_loc <- function (month, group3) {
#   symbols <- sampleinfo %>%
#     filter(remove == FALSE, Month == month, Group1 == group3) %>% pull(Symbol)
#   return(rowSums(df_pa[symbols]))
# }
# aggre_loc_all <- function (month) {
#   symbols <- sampleinfo %>%
#     filter(remove == FALSE,  Month == month) %>% pull(Symbol)
#   return(rowSums(df_pa[symbols]))
# }
aggre_loc <- function (month, group3) {
  symbols <- sampleinfo %>%
    filter(Month == month, Group1 == group3) %>% pull(Symbol)
  return(rowSums(df_pa[symbols]))
}
aggre_loc_all <- function (month) {
  symbols <- sampleinfo %>%
    filter( Month == month) %>% pull(Symbol)
  return(rowSums(df_pa[symbols]))
}

df_aggre_loc <- cbind(
  "N_R" = aggre_loc("N", "R"),
  "N_IB" = aggre_loc("N", "IB"),
  "N_OB" = aggre_loc("N", "OB"),
  "N" = aggre_loc_all("N"),
  "J_R" = aggre_loc("J", "R"),
  "J_IB" = aggre_loc("J", "IB"),
  "J_OB" = aggre_loc("J", "OB"),
  "J" = aggre_loc_all("J"))
df_aggre_loc <- replace(df_aggre_loc, df_aggre_loc > 0, 1)
# df_aggre_loc <- cbind(df_aggre_loc, "Fish.type" = df_with_fishtype$habitat_class)
tab <- df_aggre_loc %>%
  as.data.frame() %>%
  mutate(Species = rownames(.)) %>%
  left_join(select(habitat, Species, habitat_class), by = "Species") %>%
  arrange(habitat_class) %>%
  group_by(habitat_class) %>%
  summarise(across(everything(), 
                   function(x) {sum(x>0)}),
            n_samp = n())
tab
write.csv(tab, "./output/table.csv")

# ==============================================================================
# Check results
# ==============================================================================

# wilcox test
# qPCR inside the bay (high vs low)

# library(dplyr)

sampleinfo <- read.csv("./input/sample_info.csv")
high <- sampleinfo %>%
  filter(Group1 ==  "IB") %>%
  select(Site, Month, Tide, qPCR_Carp) %>%
  filter(Tide == "High")

low <- sampleinfo %>%
  filter(Group1 ==  "IB") %>%
  select(Site, Month, Tide, qPCR_Carp) %>%
  filter(Tide == "Low")

d <- left_join(high, low, by = c("Site", "Month"))
colnames(d)[colnames(d) == "qPCR_Carp.x"] <- "High"
colnames(d)[colnames(d) == "qPCR_Carp.y"] <- "Low"

wilcox.test(d$High, d$Low, paired = T)

# ======================
# qPCR S7 U/B

upper <- sampleinfo %>%
  filter(Site2 == "S7", UB == "U") %>%
  select(Month, Symbol, UB, qPCR_Carp)

bottom <- sampleinfo %>%
  filter(Site2 == "S7", UB == "B") %>%
  select(Month, Symbol, UB, qPCR_Carp)

c(mean(upper$qPCR_Carp), mean(bottom$qPCR_Carp))


# == end ==

