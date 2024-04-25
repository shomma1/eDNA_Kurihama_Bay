
set.seed(777)

library(vegan)
library(tidyverse)

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
  ggplot(aes(x = PSU, y = qPCR_Carp, col = Group3)) +
  geom_point() +
  facet_wrap(~Month) +
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

# 注意
sampleinfo <- read.csv("./input/sample_info.csv")
sampleinfo$qPCR_Carp <- sampleinfo$qPCR_Carp -1
is.na(sampleinfo$qPCR_Carp) <- 0
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

# spatial dispersion ======================================================
# library(geosphere)
dist_geo <- as.matrix(distm(sampleinfo[c("long", "lat")]))
rownames(dist_geo) <- colnames(dist_geo) <- sampleinfo$Symbol

sampleinfo %>%
  mutate(dist_from_S3 = dist_geo[,"S3.L_N"]) %>%
  mutate( Month = case_when(
    Month == "N" ~ "November",
    Month == "J" ~ "January",
    TRUE ~ Month)) %>%
  filter(Group1 == "IB") %>%
  ggplot() +
  geom_point(aes(x = dist_from_S3, y = log10(qPCR_Carp+1), 
                 col = paste(Month, Tide))) +
  geom_smooth(formula = y ~ x, method = "lm", 
              aes(x = dist_from_S3, y = log10(qPCR_Carp + 1), 
                  group = paste(Month, Tide), col = paste(Month, Tide)),
              se = T, alpha = 0.2) +
  labs(x = "Geographical distance from S3 (m)",
       y = "Log10(eDNA concentrations + 1) [copies / L]",
       col = "Month, Tide") +
  theme_minimal()

#===============================================================================
# Figure S3
# ==============================================================================

df_org <- read.csv("./input/dataset.csv", row.names = 1)
sampleinfo <- read.csv("./input/sample_info.csv")
dev.new()
# box plot of total read counts
bp <- boxplot(log10(colSums(df_org + 1)), horizontal = F, 
        range = 2, ylab = "log10(read counts + 1)")
text(names(bp$out), x = 1.1, y = bp$out, adj = 0)

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
tab_psu <- read.csv("./input/PSU/PSU_vertical.csv")

tab_psu %>%
  pivot_longer(-depth, names_to = c("survey", "location"), names_sep = "_") %>%
  filter(!is.na(value)) %>%
  group_by(survey, location) %>%
  mutate(survey = factor(survey, levels = c("NL", "NH", "JL", "JH"))) %>%
  ggplot(aes(x = value, y = depth, col = location)) +
  geom_point() +
  geom_line() +
  facet_wrap(~survey) +
  xlab("Salinity [psu] ") +
  ylab("Depth [m]")

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

res <- capscale(t(df_relative) ~ 1, distance = "bray",)

plot(res, type = "t", display = "site")

# == end ==

