####
#### R script for Ohigashi et al (2024)
#### Check effect of extraction buffer on microbial community
#### 2025.06.13 written by Ohigashi
#### R 4.3.3
####


### load packages and functions
source("Function/F2_HelperFunctions_for_Visualization.R")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(tibble); packageVersion("tibble")


### load data
# prokaryotes
b_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_16S.txt", header = T)
b_ASV <- b_ASV.table[,1:(ncol(b_ASV.table)-7)]
b_ASV.t <- t(b_ASV)

# fungi
f_ASV.table <- read.table("01_DADA2_out/rarefied_ASV_table_ITS.txt", header = T)
f_ASV <- f_ASV.table[,1:(ncol(f_ASV.table)-7)]
f_ASV.t <- t(f_ASV)

# environment data
env_data <- read.table("Data/soil_metadata.txt", header = T)

## PCoA data
# prokaryotes
b_data.scores <- readRDS("04_PCoA_PERMANOVA_out/pcoa_scores_prok_df.rds")
# fungi
f_data.scores <- readRDS("04_PCoA_PERMANOVA_out/pcoa_scores_fungi_df.rds")


### format data
# enter buffer data
sl2_samples <- c("K05", "K08", "K17", "K27", "K41")

# add buffer data to envdata
env_data <- env_data |>
  mutate(buffer = case_when(
    Sample %in% sl2_samples ~ "SL2",
    TRUE ~ "SL1"
  ))

# get buffer vector
buffer <- env_data$buffer


### PERMANOVA
# set the method to generate random values
set.seed(123)

# prokaryotes
b_perm <- adonis2(b_ASV.t~buffer, method = "bray", permutations = 100000) # permanova
b_permsummary <- as.matrix(b_perm[1,5]) # extract p-values
rownames(b_permsummary) <- c("Buffer")
colnames(b_permsummary) <- "p-value"
b_perm_result <- changep(b_permsummary) # convert the p-values to asterisks

# fungi
f_perm <- adonis2(f_ASV.t~buffer, method = "bray", permutations = 100000) # permanova
f_permsummary <- as.matrix(f_perm[1,5]) # extract p-values
rownames(f_permsummary) <- c("Buffer")
colnames(f_permsummary) <- "p-value"
f_perm_result <- changep(f_permsummary) # convert the p-values to asterisks



### plot
# set maximum values to annotate text (significance table)
b_yRoof <- (max(b_data.scores$Dim2))
f_yRoof <- (max(f_data.scores$Dim2))


# add buffer to the data frames
b_data.scores$Buffer <- buffer
f_data.scores$Buffer <- buffer

# plot
# prokaryotes
gb_buffer <- ggplot() +
  geom_point(data=b_data.scores, aes(x=Dim1, y=Dim2, color=Buffer), size=4) +
  theme_bw() +
  theme(panel.border = element_rect(size = 1.2, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 15, face = "bold", colour = "black"),
        legend.text = element_text(size = 15, face = "bold", colour = "black"),
        plot.title = element_text(size=15, face = "bold", colour = "black"),
        axis.text=element_text(size=10, face = "bold", colour = "black"),
        axis.title=element_text(size=10,face="bold", colour = "black"))+
  labs(title = "Prokaryotes", color = "Buffer",
       x = "Axis 1 (19.7%)", y = "Axis 2 (10.7%)") +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN)+
  annotate("text", x=min(b_data.scores$Dim1), y=b_yRoof,
           label=b_perm_result,
           hjust="left",
           size = 4)
plot(gb_buffer)

# fungi
gf_buffer <- ggplot() +
  geom_point(data=f_data.scores, aes(x=Dim1, y=Dim2, color=Buffer), size=4) +
  theme_bw() +
  theme(panel.border = element_rect(size = 1.2, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 15, face = "bold", colour = "black"),
        legend.text = element_text(size = 15, face = "bold", colour = "black"),
        plot.title = element_text(size=15, face = "bold", colour = "black"),
        axis.text=element_text(size=10, face = "bold", colour = "black"),
        axis.title=element_text(size=10,face="bold", colour = "black"))+
  labs(title = "Fungi", color = "Buffer",
       x = "Axis 1 (10.8%)", y = "Axis 2 (10%)") +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN)+
  annotate("text", x=max(f_data.scores$Dim1), y=f_yRoof,
           label=f_perm_result,
           hjust="right",
           size = 4)
plot(gf_buffer)


### save data
dir.create("FigCode/FigS_PCoA_check.buffer.effect_out")
# permanova result
write.csv(b_perm, file = "FigCode/FigS_PCoA_check.buffer.effect_out/permanova.buffer_result_prok.csv", quote = F, row.names = T)
write.csv(f_perm, file = "FigCode/FigS_PCoA_check.buffer.effect_out/permanova.buffer_result_fungi.csv", quote = F, row.names = T)

# plot objects
saveRDS(gb_buffer, file = "FigCode/FigS_PCoA_check.buffer.effect_out/PCoA_plot_prok.buffer.rds")
saveRDS(gf_buffer, file = "FigCode/FigS_PCoA_check.buffer.effect_out/PCoA_plot_fungi.buffer.rds")

# plot images
ggsave(filename = "FigCode/FigS_PCoA_check.buffer.effect_out/PCoA_plot_prok.buffer.png", plot = gb_buffer)
ggsave(filename = "FigCode/FigS_PCoA_check.buffer.effect_out/PCoA_plot_fungi.buffer.png", plot = gf_buffer)

