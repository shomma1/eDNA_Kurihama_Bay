
library(tidyverse)

# Preprocessing ================================================================

df_org <- read.csv("./input/dataset.csv", row.names = 1)
sampleinfo <- read.csv("./input/sample_info.csv")

# basic info
print(paste("total read counts", sum(colSums(df_org))))

# removed sample
print("removed samples")
print(sampleinfo[sampleinfo$remove == "TRUE", "Symbol"])

# number of fresh in S8
print("number of fresh in S8")
print(sampleinfo[sampleinfo$Site2 %in% "S8", "n_sp_ff"])

# number of marine sites carp detected
print("number of marine sites carp detected")
sample <- 
  sampleinfo[sampleinfo$Group3 == "S" & sampleinfo$remove != TRUE, "Symbol"]
d <- df_org["Cyprinus carpio", sample]
paste(length(colnames(d)[d>0]), "/",  ncol(d))


# ==============================================================================
# create df_aggre 

# read
df_org <- read.csv("./input/dataset.csv", row.names = 1)
df_relative <- apply(df_org, 2, function(x) {x/sum(x)})
df_relative[,is.na(colSums(df_relative))] <- rep(0, nrow(df_relative))
df_pa <- replace(df_relative, df_relative > 0, 1)

sampleinfo <- read.csv("./input/sample_info.csv")
fish_type <- read.csv("./input/habitat.csv")

# order fishes
fish_type <- fish_type %>%
  arrange(match(Species, rownames(df_org))) %>%
  rename(Fish.type = habitat_class)
  
#check
ifelse(all(rownames(df_relative) %in% fish_type$Species), "OK", "NG")
####fish_list_fresh <- fish_type[fish_type$Fish.type == 1, "Species"]

df_with_fishtype <- df_relative %>%
  data.frame(Species = rownames(.)) %>%
  left_join(fish_type, by = "Species") %>%
  select(-X)
  
# summarize
df_aggre <- df_with_fishtype %>%
  group_by(Fish.type) %>%
  summarise(across(-Species, sum))

df_aggre_sp <- df_with_fishtype %>%
  group_by(Fish.type) %>%
  summarise(across(-Species, function(x) {sum(x>0)} ))

#===============================================================================
# check

symbols <- sampleinfo %>%
  filter(remove == FALSE,
         Group3 == "R") %>%
  pull(Symbol)
print("number of F and F-B at River sites")
round(mean(rowMeans(df_aggre[symbols])[1:2]),3) # F, F&B

print("number of B-M and M at River sites")
round(mean(rowMeans(df_aggre[symbols])[5:6]),3) # M, M&B


symbols <- sampleinfo %>%
  filter(remove == FALSE,
         Group3 == "S") %>%
  pull(Symbol)
print("number of F and F-B at Marine sites")
round(mean(rowMeans(df_aggre[symbols])[1:2]),3) # F, F&B
print("number of B-M and M at Marine sites")
round(mean(rowMeans(df_aggre[symbols])[5:6]),3) # M, M&B

# ==============================================================================

sampleinfo_org <- read.csv("./input/original/sample_info.csv")

# doimnant species
dominant <- as.data.frame(cbind(apply(df_relative, 2, which.max),
                                apply(df_relative, 2, max)))
colnames(dominant) <- c("row_index", "reads")
dominant["dominant_sp"] <- rownames(df_relative)[dominant$row_index]

data <- data.frame(
  Symbol = colnames(df_org),
  ttl_reads = colSums(df_org),
  n_sp = colSums(df_org>0),
  en_sp = apply(df_org, 2, function (x) { 1 / sum( ( x/sum(x) )^2 ) }),
  Cyprinus.carpio = df_relative["Cyprinus carpio", ],
  Mugil.cephalus = df_relative["Mugil cephalus", ],
  Lateolabrax.japonicus = df_relative["Lateolabrax japonicus", ],
  Girella.punctata = df_relative["Girella punctata", ],
  dominant_sp = dominant$dominant_sp,
  dominant_sp_reads = dominant$reads)

sampleinfo <- left_join(sampleinfo_org, data, by = "Symbol")

write.csv(sampleinfo, "./output/sample_info.csv")

# finish =======================================================================




