library('ggplot2')
library(dplyr)
library(plyr)
library(ggpubr)


INV <- read.table('INV_vep.csv', header = TRUE, sep = '\t', dec = ".")
INS <- read.table('INS_vep.csv', header = TRUE, sep = '\t', dec = ".")
DEL <- read.table('DEL_vep.csv', header = TRUE, sep = '\t', dec = ".")
DUP <- read.table('DUP_vep.csv', header = TRUE, sep = '\t', dec = ".")

# analysis of the distribution of all SVs in the COVID-19 classes through box plots

INV$class = ifelse(INV$WHOCPS == "1" | INV$WHOCPS == "2" | INV$WHOCPS == "3", "1-3",
                   ifelse(INV$WHOCPS == "4" | INV$WHOCPS == "5", "4-5",
                          ifelse(INV$WHOCPS == "10", "10", "6-9")))
INV <- INV[complete.cases(INV$WHOCPS), ]

df_inv = INV %>%
  dplyr::group_by(sample_name, SVTYPE, class) %>%
  dplyr::summarise(conting =n()
  )

class_order <- c("1-3", "4-5", "6-9", "10")
custom_colors <- c("green", "lightyellow", "gold", "red")

df_inv$class <- factor(df_inv$class, levels = class_order, ordered = TRUE)
plot <- ggplot(df_inv, aes(x = class, y = conting, fill = class)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Boxplot of Inversions per COVID-19 WHO phenotypic classes") +
  theme_classic() 
  
print(plot)



# --- INS
INS$class = ifelse(INS$WHOCPS == "1" | INS$WHOCPS == "2" | INS$WHOCPS == "3", "1-3",
                   ifelse(INS$WHOCPS == "4" | INS$WHOCPS == "5", "4-5",
                          ifelse(INS$WHOCPS == "10", "10", "6-9")))
INS <- INS[complete.cases(INS$WHOCPS), ]

df_ins = INS %>%
  dplyr::group_by(sample_name, SVTYPE, class) %>%
  dplyr::summarise(conting =n()
  )

class_order <- c("1-3", "4-5", "6-9", "10")
custom_colors <- c("green", "lightyellow", "gold", "red")

df_ins$class <- factor(df_ins$class, levels = class_order, ordered = TRUE)
plot_ins <- ggplot(df_ins, aes(x = class, y = conting, fill = class)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Boxplot of Insertions per COVID-19 WHO phenotypic classes") +
  theme_classic() 

print(plot_ins)



# --- DELS
DEL$class = ifelse(DEL$WHOCPS == "1" | DEL$WHOCPS == "2" | DEL$WHOCPS == "3", "1-3",
                   ifelse(DEL$WHOCPS == "4" | DEL$WHOCPS == "5", "4-5",
                          ifelse(DEL$WHOCPS == "10", "10", "6-9")))
DEL <- DEL[complete.cases(DEL$WHOCPS), ]

df_DEL = DEL %>%
  dplyr::group_by(sample_name, SVTYPE, class) %>%
  dplyr::summarise(conting =n()
  )

class_order <- c("1-3", "4-5", "6-9", "10")
custom_colors <- c("green", "lightyellow", "gold", "red")

df_DEL$class <- factor(df_DEL$class, levels = class_order, ordered = TRUE)
plot <- ggplot(df_DEL, aes(x = class, y = conting, fill = class)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Boxplot of Deletions per COVID-19 WHO phenotypic classes") +
  theme_classic() 

print(plot)

# --- DUP
DUP$class = ifelse(DUP$WHOCPS == "1" | DUP$WHOCPS == "2" | DUP$WHOCPS == "3", "1-3",
                   ifelse(DUP$WHOCPS == "4" | DUP$WHOCPS == "5", "4-5",
                          ifelse(DUP$WHOCPS == "10", "10", "6-9")))
DUP <- DUP[complete.cases(DUP$WHOCPS), ]

df_DUP = DUP %>%
  dplyr::group_by(sample_name, SVTYPE, class) %>%
  dplyr::summarise(conting =n()
  )

class_order <- c("1-3", "4-5", "6-9", "10")
custom_colors <- c("green", "lightyellow", "gold", "red")

df_DUP$class <- factor(df_DUP$class, levels = class_order, ordered = TRUE)
plot <- ggplot(df_DUP, aes(x = class, y = conting, fill = class)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Boxplot of Duplications per COVID-19 WHO phenotypic classes") +
  theme_classic() 

print(plot)




# Split the rows per "," in column Consequence
split_string = function(lst){
  lst %>% 
    mutate(SYMBOL = strsplit(as.character(SYMBOL), ",", fixed = TRUE)) %>% 
    mutate(Consequence = strsplit(as.character(Consequence), ",", fixed = TRUE)) %>%
    tidyr::unnest(SYMBOL, Consequence)
}


inv_splitted = split_string(INV)
dup_splitted = split_string(DUP)
del_splitted = split_string(DEL)
ins_splitted = split_string(INS)


# take the coding SV
specific_words = list('coding_sequence_variant', 'inframe_deletion', 
                   'frameshift_variant', 'inframe_insertion', 'transcript_ablation', 
                   'feature_truncation', 'feature_elongation', 'stop_lost', 
                   'transcript_amplification', 'protein_altering_variant', 
                   'stop_gained', 'start_lost')


INVFiltered_coding <- inv_splitted[grepl(paste(specific_words, collapse = "|"), inv_splitted$Consequence), ]
DELFiltered_coding <- del_splitted[grepl(paste(specific_words, collapse = "|"), del_splitted$Consequence), ]
DUPFiltered_coding <- dup_splitted[grepl(paste(specific_words, collapse = "|"), dup_splitted$Consequence), ]
INSFiltered_coding <- ins_splitted[grepl(paste(specific_words, collapse = "|"), ins_splitted$Consequence), ]









# analysis of the distribution of coding SVs in the COVID-19 classes through box plots
# --- INV

INVFiltered_coding$class = ifelse(INVFiltered_coding$WHOCPS == "1" | INVFiltered_coding$WHOCPS == "2" | INVFiltered_coding$WHOCPS == "3", "1-3",
                           ifelse(INVFiltered_coding$WHOCPS == "4" | INVFiltered_coding$WHOCPS == "5", "4-5",
                                  ifelse(INVFiltered_coding$WHOCPS == "10", "10", "6-9")))
INVFiltered_coding <- INVFiltered_coding[complete.cases(INVFiltered_coding$WHOCPS), ]

df_inv = INVFiltered_coding %>%
  dplyr::group_by(sample_name, SVTYPE, class) %>%
  dplyr::summarise(conting =n()
  )

class_order <- c("1-3", "4-5", "6-9", "10")
custom_colors <- c("green", "lightyellow", "gold", "red")

df_inv$class <- factor(df_inv$class, levels = class_order, ordered = TRUE)
plot <- ggplot(df_inv, aes(x = class, y = conting, fill = class)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Boxplot of Inversions within coding regions per COVID-19 WHO phenotypic classes") +
  theme_classic() 

print(plot)



# --- INS
INSFiltered_coding$class = ifelse(INSFiltered_coding$WHOCPS == "1" | INSFiltered_coding$WHOCPS == "2" | INSFiltered_coding$WHOCPS == "3", "1-3",
                           ifelse(INSFiltered_coding$WHOCPS == "4" | INSFiltered_coding$WHOCPS == "5", "4-5",
                                  ifelse(INSFiltered_coding$WHOCPS == "10", "10", "6-9")))
INSFiltered_coding <- INSFiltered_coding[complete.cases(INSFiltered_coding$WHOCPS), ]

df_ins = INSFiltered_coding %>%
  dplyr::group_by(sample_name, SVTYPE, class) %>%
  dplyr::summarise(conting =n()
  )

class_order <- c("1-3", "4-5", "6-9", "10")
custom_colors <- c("green", "lightyellow", "gold", "red")

df_ins$class <- factor(df_ins$class, levels = class_order, ordered = TRUE)
plot_ins <- ggplot(df_ins, aes(x = class, y = conting, fill = class)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Boxplot of Insertions eithin coding regions per COVID-19 WHO phenotypic classes") +
  theme_classic() 

print(plot_ins)



# --- DELS
DELFiltered_coding$class = ifelse(DELFiltered_coding$WHOCPS == "1" | DELFiltered_coding$WHOCPS == "2" | DELFiltered_coding$WHOCPS == "3", "1-3",
                           ifelse(DELFiltered_coding$WHOCPS == "4" | DELFiltered_coding$WHOCPS == "5", "4-5",
                                  ifelse(DELFiltered_coding$WHOCPS == "10", "10", "6-9")))
DELFiltered_coding <- DELFiltered_coding[complete.cases(DELFiltered_coding$WHOCPS), ]

df_DEL = DELFiltered_coding %>%
  dplyr::group_by(sample_name, SVTYPE, class) %>%
  dplyr::summarise(conting =n()
  )

class_order <- c("1-3", "4-5", "6-9", "10")
custom_colors <- c("green", "lightyellow", "gold", "red")

df_DEL$class <- factor(df_DEL$class, levels = class_order, ordered = TRUE)
plot <- ggplot(df_DEL, aes(x = class, y = conting, fill = class)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  stat_compare_means(comparisons = class)+
  stat_compare_means(label.y = 50) +
  labs(title = "Boxplot of Deletions within coding regions per COVID-19 WHO phenotypic classes") +
  theme_classic() 
print(plot)



# --- DUP
DUP <- read.table('/Users/gaiacorona/Downloads/DUP_vep.csv', header = TRUE, sep = '\t', dec = ".")
DUPFiltered_coding$class = ifelse(DUPFiltered_coding$WHOCPS == "1" | DUPFiltered_coding$WHOCPS == "2" | DUPFiltered_coding$WHOCPS == "3", "1-3",
                           ifelse(DUPFiltered_coding$WHOCPS == "4" | DUPFiltered_coding$WHOCPS == "5", "4-5",
                                  ifelse(DUPFiltered_coding$WHOCPS == "10", "10", "6-9")))
DUPFiltered_coding <- DUPFiltered_coding[complete.cases(DUPFiltered_coding$WHOCPS), ]

df_DUP = DUPFiltered_coding %>%
  dplyr::group_by(sample_name, SVTYPE, class) %>%
  dplyr::summarise(conting =n()
  )

class_order <- c("1-3", "4-5", "6-9", "10")
custom_colors <- c("green", "lightyellow", "gold", "red")

df_DUP$class <- factor(df_DUP$class, levels = class_order, ordered = TRUE)
plot <- ggplot(df_DUP, aes(x = class, y = conting, fill = class)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Boxplot of Duplications within coding regions per COVID-19 WHO phenotypic classes") +
  theme_classic() 

print(plot)
ggsave("DUP_total.png", plot, width = 10, height = 5, dpi = 600)





# ANALYSIS ON CLASS 10: see the genomic coordinates of the variants involved in genes selected by the collapse

class <- "10"

class10_del_coding = DELFiltered_coding %>%
  filter(WHOCPS == class)
class10_dup_coding = DUPFiltered_coding %>%
  filter(WHOCPS == class)
class10_ins_coding = INSFiltered_coding %>%
  filter(WHOCPS == class)
class10_inv_coding = INVFiltered_coding %>%
  filter(WHOCPS == class)


# take the excel tables from the collapse analysis 
del_collapse = read_excel("DEL_coding.xlsx")
dup_collapse = read_excel("DUP_coding.xlsx")
ins_collapse = read_excel("INS_coding.xlsx")
inv_collapse = read_excel("INV_coding.xlsx")

# remove all the rows for which the last col is empty
del_collapse <- del_collapse[!is.na(del_collapse[, ncol(del_collapse)]), ]
dup_collapse <- dup_collapse[!is.na(dup_collapse[, ncol(dup_collapse)]), ]
ins_collapse <- ins_collapse[!is.na(ins_collapse[, ncol(ins_collapse)]), ]
inv_collapse <- inv_collapse[!is.na(inv_collapse[, ncol(inv_collapse)]), ]

# select from the collapse analysis table all the genes that are present in the class 10 tables

del_collapse_10 = del_collapse[del_collapse$SYMBOL %in% class10_del_coding$SYMBOL, ]
dup_collapse_10 = dup_collapse[dup_collapse$SYMBOL %in% class10_dup_coding$SYMBOL, ]
ins_collapse_10 = ins_collapse[ins_collapse$SYMBOL %in% class10_ins_coding$SYMBOL, ]
inv_collapse_10 = inv_collapse[inv_collapse$SYMBOL %in% class10_inv_coding$SYMBOL, ]

# select from the class 10 dataframe all the rows that have the genes as in the dataframe now created
# so that we know the positions of the variants in those genes
class10_del_genes = class10_del_coding[class10_del_coding$SYMBOL %in% del_collapse_10$SYMBOL, ]
class10_dup_genes = class10_dup_coding[class10_dup_coding$SYMBOL %in% dup_collapse_10$SYMBOL, ]
class10_ins_genes = class10_ins_coding[class10_ins_coding$SYMBOL %in% ins_collapse_10$SYMBOL, ]













