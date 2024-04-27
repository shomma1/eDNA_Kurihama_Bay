
library(vegan)
library(tidyverse)
library(geosphere)

#==================================================================
# Figure 2 --------------------------------------------------------
#==================================================================

# geographical distance - C carpio eDNA
sampleinfo <- read.csv("./input/sample_info.csv")

tides <- c("High", "Low")
cols <- c("salmon", "royalblue") 
pchs <- c(16, 16)
dev.new()
par(mfrow = c(1,2))
for (month in c("Nov", "Jan")) {
  for (i in 1:2) {
    
    tide <- tides[[i]]
    mon <- substr(month, 1, 1)
    sample <- sampleinfo[sampleinfo$Tide == tide & sampleinfo$Month == mon,]
    dist_geo <- distm(sample[c("long", "lat")])
    rownames(dist_geo) <- colnames(dist_geo) <- sample$Site
    
    plot(dist_geo["RA",], log10(sample$qPCR_Carp),  
         col = cols[i], pch = pchs[i], cex = 2,
         xlim = c(0, 6000), ylim = c(0,6),
         xlab = "Geographical distance from RA [m]",
         ylab = "Log 10 (eDNA concentration + 1) [copies / L]")
    axis(side=3, dist_geo["RA",], rownames(dist_geo))
    mtext(month, side = 3, line = 2, cex = 1.6, adj = 0)
    
    # linear regression in River ===
    y  <- sample %>% filter(Group1 == "R") %>%
        select(Site, qPCR_Carp)
    x <- dist_geo["RA", y$Site]
    model <- lm(log10(y$qPCR_Carp + 1) ~ x)
    p_value <- summary(model)$coefficients["x","Pr(>|t|)"]
    
    if (p_value < 0.05) { # significant
      abline(model, col=cols[i], lty = 2)
      y_pred <- predict(model, data.frame(x))
      lines(x, y=y_pred, col = cols[i],lwd =2)
    
    } else { # not significant -> null model was accepted
      model <- lm(log10(y$qPCR_Carp + 1) ~ 1) # null model
      abline(model, col=cols[i], lty = 2)
      y_pred <- predict(model, data.frame(x))
      lines(x, y=y_pred, col = cols[i], lwd = 2)
    }
    
    # plot
    y <- sample %>% filter(Group1 == "IB") %>%
      select(Site, qPCR_Carp)
    mean_y_bay <- mean(log10(y$qPCR_Carp +1))
    lines(c(-100,dist_geo["RA","S7B"]), # S2 to S7
          c(mean_y_bay,mean_y_bay),
          col=cols[i], lty = 3)
    mean_x_bay <- mean(dist_geo["RA", y$Site])
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

dev.new()
plot.new()
legend("center", c("High", "Low"), pch = pchs, col = cols, title = "Tide")

#==================================================================
# Figure 3 --------------------------------------------------------
#==================================================================
# need to run processing and get object "df_aggre_sp", "df_aggre" before

sampleinfo <- read.csv("./input/sample_info.csv")

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

#dev.new()
par(mfrow = c(3,2))
# F, F-B, B, F-B-M, B-M, M
cols = c("red", "salmon","green3", "lightgreen","skyblue","royalblue")

# Number of Species
# November
sample <- sampleinfo[sampleinfo$Month == "N" , 
                     c("Symbol", "Site", "Loc", "Tide_1")]
plt <- df_aggre_sp[, sample$Symbol]
plt[, removed] <- 0
bp <- barplot(as.matrix(plt), names.arg = sample$Tide_1,
              las = 2, ylim = c(0, 70), col = cols)
title(main = "Number of species", line = 2, adj = 0, cex.main = 1.5)
mtext("(a) November", side = 3, line = 0.5, adj = 0)
remove <- match(removed, colnames(plt))
text(bp[remove], 0, " excluded", srt = 90, adj = 0)
new_x <- (bp[seq(1,length(bp),2)] + bp[1 + seq(1,length(bp),2)]) / 2
axis(side = 1, at = new_x, labels = unique(sample$Loc), tick = F,
     las = 2, line = 0.9, cex.axis = 1.3)
axis(side = 1, at = bp, tick = T, labels = F)


# January
sample <- sampleinfo[sampleinfo$Month == "J" , 
                     c("Symbol", "Site", "Loc", "Tide_1")]
plt <- df_aggre_sp[, sample$Symbol]
barplot(as.matrix(plt), names.arg = sample$Tide_1, 
        las = 2, ylim = c(0, 70), col = cols)
mtext("(b) January", side = 3, line = 0.5, adj = 0)
new_x <- (bp[seq(1,length(bp),2)] + bp[1 + seq(1,length(bp),2)]) / 2
axis(side = 1, at = new_x, labels = unique(sample$Loc), tick = F,
     las = 2, line = 0.9, cex.axis = 1.3)
axis(side = 1, at = bp, tick = T, labels = F)

# Effective number of species ===
# November
sample <- sampleinfo[sampleinfo$Month == "N", c("Symbol", "Site")]
plt <- sampleinfo[match(sample$Symbol, sampleinfo$Symbol), 
                  c("Symbol", "en_sp", "Loc", "Tide_1")]
plt[plt$Symbol %in% removed, "en_sp"] <- 0
bp <- barplot(plt$en_sp, names.arg = plt$Tide_1, las = 2, ylim = c(0, 14))
title(main = "Effective number of species", line = 2, adj = 0, cex.main = 1.5)
mtext("(c) November", side = 3, line = 0.5, adj = 0)
remove <- match(removed, plt$Symbol)
text(bp[remove], 0, " excluded", srt = 90, adj = 0)
new_x <- (bp[seq(1,length(bp),2)] + bp[1 + seq(1,length(bp),2)]) / 2
axis(side = 1, at = new_x, labels = unique(plt$Loc), tick = F,
     las = 2, line = 0.9, cex.axis = 1.3)
axis(side = 1, at = bp, tick = T, labels = F)

# January
sample <- sampleinfo[sampleinfo$Month == "J", c("Symbol", "Site")]
plt <- sampleinfo[match(sample$Symbol, sampleinfo$Symbol), 
                  c("Symbol", "en_sp", "Loc", "Tide_1")]
barplot(plt$en_sp, names.arg = plt$Tide_1, las = 2, ylim = c(0, 14))
mtext("(d) January", side = 3, line = 0.5, adj = 0)
new_x <- (bp[seq(1,length(bp),2)] + bp[1 + seq(1,length(bp),2)]) / 2
axis(side = 1, at = new_x, labels = unique(plt$Loc), tick = F,
     las = 2, line = 0.9, cex.axis = 1.3)
axis(side = 1, at = bp, tick = T, labels = F)

# relative read counts ===
# November
sample <- sampleinfo[sampleinfo$Month == "N" , 
                     c("Symbol", "Site", "Loc", "Tide_1")]
plt <- df_aggre[, sample$Symbol]
plt[, removed] <- 0
bp <- barplot(as.matrix(plt),names.arg = sample$Tide_1, las = 2, col = cols)
title(main = "Relative read counts", line = 2, adj = 0, cex.main = 1.5)
mtext("(e) November", side = 3, line = 0.5, adj = 0)
remove <- match(removed, colnames(plt))
text(bp[remove], 0, " excluded", srt = 90, adj = 0)
new_x <- (bp[seq(1,length(bp),2)] + bp[1 + seq(1,length(bp),2)]) / 2
axis(side = 1, at = new_x, labels = unique(sample$Loc), tick = F,
     las = 2, line = 0.9, cex.axis = 1.3)
axis(side = 1, at = bp, tick = T, labels = F)

# January
sample <- sampleinfo[sampleinfo$Month == "J" ,
                     c("Symbol", "Site", "Loc", "Tide_1")]
plt <- df_aggre[, sample$Symbol]
barplot(as.matrix(plt), names.arg = sample$Tide_1, las = 2, col = cols)
mtext("(f) January", side = 3, line = 0.5, adj = 0)
new_x <- (bp[seq(1,length(bp),2)] + bp[1 + seq(1,length(bp),2)]) / 2
axis(side = 1, at = new_x, labels = unique(sample$Loc), tick = F,
     las = 2, line = 0.9, cex.axis = 1.3)
axis(side = 1, at = bp, tick = T, labels = F)

# legend
dev.new()
plot.new()
legend("center", rev(c("F", "F-B","F-B-M", "B", "B-M", "M")), fill = rev(cols))

#==================================================================
# Figure 4 --------------------------------------------------------
#==================================================================
# dominant
sampleinfo <- read.csv("./input/sample_info.csv")
remove_list <- c("S3.L_N", "S4.L_N", "S5.H_N", "S7U.H_N", "S7B.H_N")
sampleinfo <- sampleinfo %>%
  mutate(Tide_1 = substring(Tide,1,1)) %>%
  mutate(Loc = Site) %>%
  mutate(Site = paste0(Site, Tide_1))


plot_bar <- function(month, sp, title, text) {
  sample <- sampleinfo[sampleinfo$Month == month, c("Symbol")]
  plt <- sampleinfo[match(sample, sampleinfo$Symbol), 
                    c("Symbol", sp, "Site", "Tide_1", "Loc")]
  bp <- barplot(t(cbind(plt[sp], 1 - plt[sp])), names.arg = plt$Tide_1, las = 2,
                col = c("salmon", "white"))
  title(bquote(italic(.(title))), line = 2, adj = 0, cex.main = 1.5)
  mtext(text, side = 3, line = 0.5, cex = 0.8, adj = 0)
  
  detect <- plt[sp] != 0
  points(bp[detect], rep(0.9, sum(detect)), pch = 8, col = "salmon")
  
  remove <- match(remove_list, plt$Symbol)
  rect(bp[remove] - 0.5, 0, bp[remove] + 0.5, 1, col = "black", border = NA)
  
  axis(side = 1, at = new_x, labels = unique(plt$Loc), tick = F,
       las = 2, line = 0.9, cex.axis = 1.3)
  axis(side = 1, at = bp, tick = T, labels = F)
}

#dev.new()
par(mfrow = c(3,2))
# month, sp, title, text
plot_bar("N", "Mugil.cephalus",        "M. cephalus",   "(a) November")
plot_bar("J", "Mugil.cephalus",        "M. cephalus",   "(b) January")
plot_bar("N", "Cyprinus.carpio",       "C. carpio",     "(c) November")
plot_bar("J", "Cyprinus.carpio",       "C. carpio",     "(d) January")
plot_bar("N", "Girella.punctata",      "G. punctata",   "(e) November")
plot_bar("J", "Lateolabrax.japonicus", "L. japonicus",  "(f) January")

#==================================================================
# Table 1 --------------------------------------------------------
#==================================================================
# need to run processing and get object "df_pa" "df_with_fishtype" before

# for table1
df_pa <- as.data.frame(df_pa)

aggre_loc <- function (month, group3) {
  symbols <- sampleinfo %>%
    filter(remove == FALSE, Month == month, Group3 == group3) %>% pull(Symbol)
  return(rowSums(df_pa[symbols]))
}
aggre_loc_all <- function (month) {
  symbols <- sampleinfo %>%
    filter(remove == FALSE,  Month == month) %>% pull(Symbol)
  return(rowSums(df_pa[symbols]))
}

df_aggre_loc <- cbind(
  "N_R" = aggre_loc("N", "R"),
  "N_S" = aggre_loc("N", "S"),
  "N" = aggre_loc_all("N"),
  "J_R" = aggre_loc("J", "R"),
  "J_S" = aggre_loc("J", "S"),
  "J" = aggre_loc_all("J"))
df_aggre_loc <- replace(df_aggre_loc, df_aggre_loc > 0, 1)
df_aggre_loc <- cbind(df_aggre_loc, "Fish.type" = df_with_fishtype$Fish.type)
tab <- df_aggre_loc %>%
  as.data.frame() %>%
  group_by(Fish.type) %>%
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

