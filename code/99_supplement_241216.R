
set.seed(777)

library(vegan)
library(tidyverse)
library(geosphere)

#===============================================================================
# Figure S1
# ==============================================================================

# Salinity - # qPCR 
sampleinfo <- read.csv("./input/sample_info.csv")

sampleinfo %>%
  mutate(
    Group3 = case_when(
      Group3 == "R" ~ "River",
      Group3 == "S" ~ "Marine",
      TRUE ~ Group3),
    Month = case_when(
      Month == "N" ~ "November",
      Month == "J" ~ "January",
      TRUE ~ Month),
    Group3 = factor(Group3, levels = c("River", "Marine")),
    Month = factor(Month, levels = c("November", "January"))) %>%
  ggplot(aes(x = PSU, y = qPCR_Carp + 1, col = Group3, pch = Tide)) +
  geom_point() +
  facet_wrap(~Month) +
  geom_text(aes(label = Symbol), hjust=1.1) +
  scale_y_continuous(trans = "log10") +
  labs(color = "Location",
       x = "Salinity [psu]",
       y = "eDNA concentration [copies / L]")


# Figure S3 ================================================================
df_org <- read.csv("./input/dataset.csv", row.names = 1)
# box plot of total read counts
boxp <- boxplot(log10(colSums(df_org + 1)), horizontal = T, 
        range = 2,
        main = "Log10-transformed total read counts")
text(min(boxp$out), 0.75, paste(names(boxp$out),collapse = ", "),
     pos = 4, cex = 0.7)


#===============================================================================
# Figure S2 
# ==============================================================================

# prepare for spatial analysis --------------------------------------------

library(ggmap)
register_stadiamaps("d70e9693-b38d-4d02-a8f2-0598f7f733c0", write = FALSE)
map <- get_stadiamap(bbox = c(left = 139.712,
                              bottom = 35.219,
                              right =139.725,
                              top = 35.23),
                     maptype = "stamen_toner_lite",
                     zoom = 17, color =  "bw")

sampleinfo <- read.csv("./input/sample_info.csv")
d <- sampleinfo %>%
  mutate(
    Group3 = case_when(
      Group3 == "R" ~ "River",
      Group3 == "S" ~ "Marine",
      TRUE ~ Group3),
    Month = case_when(
      Month == "N" ~ "November",
      Month == "J" ~ "January",
      TRUE ~ Month),
    Group3 = factor(Group3, levels = c("River", "Marine")),
    Month = factor(Month, levels = c("November", "January"))) %>%
  filter(Group1 == "IB") %>%
  mutate(qPCR_class = 
           cut(qPCR_Carp, breaks = c(-100, 10, 100, 1000, Inf), 
               labels = c("0-9", "10-99", "100-999", "1000+"))) %>%
  mutate(UB_size = case_when(
    UB == "U" ~ 1,
    UB == "B" ~ 2))

ggmap(map) +
  geom_point(aes(x = long, y = lat, fill = qPCR_class, size = UB_size),
             shape = 21, data = d) +
  facet_wrap(~ Month + Tide) +
  scale_fill_manual(values = c("0-9" = "blue", "10-99" = "purple", 
                               "100-999" = "#FA0101", "1000+" = "#FABA00")) +
  ggtitle("eDNA concentrations of C. carpio [copies / L]") +
  labs(fill = "copies / L") +
  guides(size = FALSE) +
  scale_size_continuous(range = c(3,6)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# spatial dispersion ===========================================================
# library(geosphere)
dist_geo <- as.matrix(distm(sampleinfo[c("long", "lat")]))
rownames(dist_geo) <- colnames(dist_geo) <- sampleinfo$Symbol


data <- sampleinfo %>%
  mutate(dist_from_S3 = dist_geo[,"S3.L_N"]) %>%
  mutate( Month = case_when(
    Month == "N" ~ "November",
    Month == "J" ~ "January",
    TRUE ~ Month)) %>%
  filter(Group1 == "IB")

col <- matrix(c("red", "skyblue", "darkgreen", "purple"), ncol= 2)
p <- NULL
for (i in 1:2) {
  ifelse(i == 1, month <- "November", month <- "January")
  for (j in 1:2) {
    ifelse(j == 1, tide <- "High", tide <- "Low")
    
    d <- data[data$Month == month & data$Tide == tide, ]
    x <- d$dist_from_S3
    y <- log10(d$qPCR_Carp+1)
    model <- lm(y ~ x)
    
    if(i == 1 & j == 1) {
      dev.new()
      plot(x, y, xlim = c(0, 650), ylim = c(0, 5),
           col = col[i,j], pch = 20, cex = 3,
           ylab = "Log10(Carp's eDNA + 1) [copies/L]",
           xlab = "Distance from S3 [m]")
    } else {
      points(x, y, xlim = c(0, 650), ylim = c(0, 5),
             col = col[j,i], pch = 20, cex = 3)
    }
    abline(model, col = col[j,i], lwd = 3)
    sum <- summary(model)
    p <- append(p, sum$coefficients["x", "Pr(>|t|)"])

  }
}
p <- round(p, 2)
dev.new()
plot.new()
legend("center", paste(c("Nov-High (p=", "Nov-Low (p=", "Jan-High (p=", "Jan-Low (p="),p, ")"),
       pch = 20, col = as.vector(col), cex = 3)


# ==============================================================================

#===============================================================================
# Figure S3
# ==============================================================================

# total read counts
df_org <- read.csv("./input/dataset.csv", row.names = 1)
sampleinfo <- read.csv("./input/sample_info.csv")
dev.new()
# box plot of total read counts
bp <- boxplot(log10(colSums(df_org + 1)), horizontal = F, 
        range = 1.5, ylab = "log10(read counts + 1)")
text(names(bp$out), x = 1.1, y = bp$out, adj = 0)

###
bp <- boxplot(log10(colSums(df_org)+1), horizontal = F, 
              range = 1.5, ylab = "log10(read counts + 1)")
text(names(bp$out), x = 1.1, y = bp$out, adj = 0)
###

#===============================================================================
# Figure S4
# ==============================================================================

# dominant
# November
dev.new()
par(mfrow = c(2,1))

removed <- sampleinfo[sampleinfo$remove == T, "Symbol"]

sample <- sampleinfo[sampleinfo$Month == "N" , "Symbol"]

plt <- sampleinfo[match(sample, sampleinfo$Symbol), 
                  c("Symbol", "dominant_sp", "dominant_sp_reads")]
plt$others <- 1 - plt[3]
plt[plt$Symbol %in% removed, "dominant_sp"] <- ""
bp <- barplot(t(plt[,3:4]), names.arg = plt$Symbol, las = 2, 
              col = c("darkgray", "white"))
title("Dominant species and its relative read counts", line = 2, adj = 0)
mtext("(a) November", side = 3, line = 0.5, cex = 0.8, adj = 0)
remove_list <- c("S3.L_N", "S4.L_N", "S5.H_N", "S7U.H_N", "S7B.H_N")
remove <- match(remove_list, plt$Symbol)
rect(bp[remove] - 0.5, 0, bp[remove] + 0.5, 1, col = "black", border = NA)
text(bp+0.5, 1, labels = plt$dominant_sp, col = "red",
     cex = 0.8, srt = 45, font.main = 2, adj = 0)

#January
sample <- sampleinfo[sampleinfo$Month == "J" , "Symbol"]
plt <- sampleinfo[match(sample, sampleinfo$Symbol), 
                  c("Symbol", "dominant_sp", "dominant_sp_reads")]
plt$others <- 1 - plt[3]
plt[plt$Symbol %in% removed, "dominant_sp"] <- ""
bp <- barplot(t(plt[,3:4]), names.arg = plt$Symbol, las = 2,
              col = c("darkgray", "white"))
mtext("(b) January", side = 3, line = 0.5, cex = 0.8, adj = 0)
remove_list <- c("S3.L_N", "S4.L_N", "S5.H_N", "S7U.H_N", "S7B.H_N")
remove <- match(remove_list, plt$Symbol)
rect(bp[remove] - 0.5, 0, bp[remove] + 0.5, 1, col = "black", border = NA)
text(bp, 1, labels = plt$dominant_sp, col = "red",
     cex = 0.8, srt = 45, font.main = 2, adj = 0)

#===============================================================================
# Figure S6 
# ==============================================================================

# PSU vertical

# (a)
# Salinity
sampleinfo <- read.csv("./input/sample_info.csv")
sampleinfo %>%
  group_by(Month) %>%
  mutate(Group1 = factor(Group1, levels = c("R", "IB", "OB")),
         Month = factor(Month, levels = c("N", "J"))) %>%
  ggplot(aes(x = Group1, y = PSU)) +
  geom_boxplot() +
  facet_wrap(~Month) +
  xlab("") +
  ylab("Salinity [psu]")


# Water temp
sampleinfo %>%
  group_by(Month) %>%
  mutate(Group1 = factor(Group1, levels = c("R", "IB", "OB")),
         Month = factor(Month, levels = c("N", "J"))) %>%
  ggplot(aes(x = Group1, y = WT)) +
  geom_boxplot() +
  facet_wrap(~Month) +
  xlab("") +
  ylab("Water temperature [°C]")


# (b)
tab_psu <- read.csv("./input/salinity/PSU_vertical.csv")

tab_psu %>%
  pivot_longer(-depth, names_to = c("survey", "location"), names_sep = "_") %>%
  filter(!is.na(value)) %>%
  group_by(survey, location) %>%
  mutate(survey = factor(survey, levels = c("NL", "NH", "JL", "JH"))) %>%
  arrange(depth) %>%
  ggplot(aes(x = value, y = depth, col = location)) +
  geom_point() +
  geom_path() +
  facet_wrap(~survey) +
  xlab("Salinity [psu] ") +
  ylab("Depth [m]")

#===============================================================================
# Heatmap
#===============================================================================
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

# order
habitat <- ttl_reads %>%
  left_join(habitat, by = "Species") %>%
  arrange(habitat_class, ttl_reads)

df_filtered <- df_org[rowSums(df_org>=1000)>=1, ]

g.list <- list()
for (i in 1:2) {
  
  month <- ifelse(i == 1, month <- "N", month <- "J")
  sample <- sampleinfo[sampleinfo$Month == month, "Symbol"]
  
  d <- df_filtered[sample]
  d <- d[rownames(d) %in%  habitat$Species, ]
  
  d_long <- d %>%
    mutate(Species = rownames(.)) %>%
    pivot_longer(
      cols = -Species,
      names_to = "Sample",
      values_to = "n_reads") %>%
#     filter(n_reads >= 1000) %>%
    mutate(color_group = case_when(
      n_reads == 0 ~ "0 ~ 9",
      n_reads > 0 & n_reads < 100 ~ "10 ~ 99",
      n_reads >= 100 & n_reads < 1000 ~ "100 ~ 999",
      n_reads >= 1000 & n_reads < 10000 ~ "1000 ~ 9999",
      n_reads >= 10000 ~ "10000 ~"))
  
  d_long <- d_long %>%
    left_join(habitat, by = "Species") %>%
    mutate(Species = factor(Species, levels = habitat$Species))
  
  # dev.new()
  g.list[[i]] <- ggplot(d_long, aes(x = Sample, y = Species, fill = color_group, labels = habitat_class)) +
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
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides(fill = guide_legend(
      override.aes = list(color = "black")
    ))
  
  dev.new()
  print(g.list[[i]])
}
# ==============================================================================
# Fresh water 
# ==============================================================================

ttl_reads <- data.frame(Species = rownames(df_org),
                        ttl_reads = rowSums(df_org))

capture_hirasaku_river <- c(
  "Cyprinus carpio",
  "Tridentiger obscurus",
  "Tridentiger sp.1", # 置き換え
  "Anguilla japonica",
  "Favonigobius gymnauchen / Acanthogobius lactipes",
  "Eleotris oxycephala",
  "Gymnogobius urotaenia", # < under 10
  "Rhinogobius giurinus",
  "Mugil cephalus", # Mugil cephalus cephalus
  "Mugilogobius abei",
  "Carassius sp.",
  "Carassius spp."
)


# habitatの並び替え．
habitat_ordered <- ttl_reads %>%
  left_join(habitat, by = "Species") %>%
  #   filter(Species %in% rownames(df_org)) %>%
  arrange(habitat_class, ttl_reads)

df_with_fishtype <- df_relative[habitat_ordered$Species, ]
df_with_fishtype <- cbind(df_with_fishtype, habitat_ordered)

# order
habitat <- ttl_reads %>%
  left_join(habitat, by = "Species") %>%
  arrange(habitat_class, ttl_reads)

df_filtered <- df_org %>%
  filter(rownames(.) %in% capture_hirasaku_river)

g.list <- list()
for (i in 1:2) {
  
  month <- ifelse(i == 1, month <- "N", month <- "J")
  month_plot <- ifelse(i == 1, tmp <- "(a) November", tmp <- "(b) January")
  sample <- sampleinfo[sampleinfo$Month == month, "Symbol"]
  
  d <- df_filtered[,sample]
  d <- d[rownames(d) %in%  habitat$Species, ]
  
  d_long <- d %>%
    mutate(Species = rownames(.)) %>%
    pivot_longer(
      cols = -Species,
      names_to = "Sample",
      values_to = "n_reads") %>%
    #     filter(n_reads >= 1000) %>%
    mutate(color_group = case_when(
      n_reads == 0 ~ "0 ~ 9",
      n_reads > 0 & n_reads < 100 ~ "10 ~ 99",
      n_reads >= 100 & n_reads < 1000 ~ "100 ~ 999",
      n_reads >= 1000 & n_reads < 10000 ~ "1000 ~ 9999",
      n_reads >= 10000 ~ "10000 ~"))
  
  d_long <- d_long %>%
    left_join(habitat, by = "Species") %>%
    mutate(Species = factor(Species, levels = habitat$Species))
  
  # dev.new()
  g.list[[i]] <- ggplot(d_long, aes(x = Sample, y = Species, fill = color_group, labels = habitat_class)) +
    geom_tile() +
    scale_fill_manual(
      values = c("0 ~ 9" = "white", 
                 "10 ~ 99" = "yellow",
                 "100 ~ 999" = "orange",
                 "1000 ~ 9999" = "red",
                 "10000 ~" = "purple")) +
    theme_minimal() +
    labs(title = month_plot, x = "Sample", y = "Species",
         fill = "Relative reads") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides(fill = guide_legend(
      override.aes = list(color = "black")
    ))
  
  dev.new()
  print(g.list[[i]])
}

#===============================================================================
# Supplementary info B
# ==============================================================================
# Positive Field Blank 

df_org <- read.csv("./input/dataset.csv", row.names = 1)
sampleinfo <- read.csv("./input/sample_info.csv")

bl <- read.csv("./input/blank/positive_field_blank.csv")

sample <- sampleinfo %>%
  filter(Month == "N", remove != T) %>%
  pull(Symbol)

df <- df_org[sample] %>%
  mutate(Species = rownames(.)) %>%
  left_join(bl, by = "Species")

rownames(df) <- df$Species
df <- df[colnames(df) !=  "Species"]

df[is.na(df)] <- 0
df <- df[rowSums(df) > 0, ]

df_relative <- apply(df, 2, function(x) { x/ sum(x) })

##
df_pa <- replace(df_relative, df_relative > 0, 1)
##

dist_sor <- vegdist(t(df_pa), distance ="bray")
dist_bray <- vegdist(t(df_relative), distance ="bray")


res <- capscale(t(df_relative) ~ 1, distance = "bray",)

plot(res, type = "t", display = "site")

# == end ==

