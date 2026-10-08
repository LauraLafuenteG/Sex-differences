library(readxl)
library(tidyverse)
library(ggtext)

# ============================================================
# Read data
# ============================================================

df <- read_excel("heatmap_total_all_N8192.xlsx")

# ============================================================
# Convert to long format
# ============================================================

data_long <- df %>%
  pivot_longer(
    cols = -Output,
    names_to = "Parameter",
    values_to = "ST"
  )

# Preserve original ordering from Excel

data_long$Output <- factor(
  data_long$Output,
  levels = rev(df$Output)
)

data_long$Parameter <- factor(
  data_long$Parameter,
  levels = names(df)[-1]
)

# ============================================================
# Parameter labels
# Sex-specific parameters are shown in bold
# ============================================================

output_labels <- c(
  "M1 peak"           = "M<sub>1</sub> peak",
  "M2 peak"           = "M<sub>2</sub> peak",
  "c1 peak"           = "c<sub>1</sub> peak",
  "c2 peak"           = "c<sub>2</sub> peak",
  "Cs peak"           = "C<sub>s</sub> peak",
  "Cb peak"           = "C<sub>b</sub> peak",
  "mf peak"           = "m<sub>f</sub> peak",
  "mb peak"           = "m<sub>b</sub> peak",
  
  "M1 time-to-peak"   = "M<sub>1</sub> TTP",
  "M2 time-to-peak"   = "M<sub>2</sub> TTP",
  "c1 time-to-peak"   = "c<sub>1</sub> TTP",
  "c2 time-to-peak"   = "c<sub>2</sub> TTP",
  "Cs time-to-peak"   = "C<sub>s</sub> TTP",
  "Cb time-to-peak"   = "C<sub>b</sub> TTP",
  "mf time-to-peak"   = "m<sub>f</sub> TTP"
)

param_labels <- c(
  
  # sex-specific
  "kps"  = "<b>k<sub>ps</sub></b>",
  "Mmax" = "<b>M<sub>max</sub></b>",
  "ds"   = "<b>d<sub>s</sub></b>",
  "k01"  = "<b>k<sub>01</sub></b>",
  "k02"  = "<b>k<sub>02</sub></b>",
  "ke1"  = "<b>k<sub>e1</sub></b>",
  "ke2"  = "<b>k<sub>e2</sub></b>",
  "k1"   = "<b>k<sub>1</sub></b>",
  "k2"   = "<b>k<sub>2</sub></b>",
  "k3"   = "<b>k<sub>3</sub></b>",
  "kls"  = "<b>k<sub>ls</sub></b>",
  "klb"  = "<b>k<sub>lb</sub></b>",
  "kpb"  = "<b>k<sub>pb</sub></b>",
  
  # non-sex-specific
  "dc1"  = "d<sub>c1</sub>",
  "k0"   = "k<sub>0</sub>",
  "aps"  = "a<sub>ps</sub>",
  "dc2"  = "d<sub>c2</sub>",
  "asb1" = "a<sub>sb1</sub>",
  "a12"  = "a<sub>12</sub>",
  "aps1" = "a<sub>ps1</sub>",
  "db"   = "d<sub>b</sub>",
  "d2"   = "d<sub>2</sub>",
  "qbd"  = "q<sub>bd</sub>",
  "pbs"  = "p<sub>bs</sub>",
  "pcs"  = "p<sub>cs</sub>",
  "d1"   = "d<sub>1</sub>",
  "qcd1" = "q<sub>cd1</sub>",
  "k12"  = "k<sub>12</sub>",
  "k21"  = "k<sub>21</sub>",
  "a22"  = "a<sub>22</sub>",
  "aed"  = "a<sub>ed</sub>",
  "a02"  = "a<sub>02</sub>",
  "d0"   = "d<sub>0</sub>",
  "qcd2" = "q<sub>cd2</sub>",
  "a01"  = "a<sub>01</sub>",
  "apb"  = "a<sub>pb</sub>"
)

# ============================================================
# Heatmap
# ============================================================

ggplot(
  data_long,
  aes(
    x = Parameter,
    y = Output,
    fill = ST
  )
) +
  
  # geom_tile(color = "grey95", linewidth = 0.15) +
  geom_tile() +
  
  scale_fill_gradientn(
     colours = c("#f7fbff","#c6dbef","#6baed6","#2171b5","#08306b"), # white to dark blue
    # colours = c("#f7fcfd","#ccece6","#66c2a4","#238b45","#00441b"), #teal-blue
    # colours = c("#f7fcfd", "#d0d1e6","#a6bddb","#3690c0","#034e7b"), #blue-purple
    name = NULL
  ) +
  
  scale_x_discrete(labels = param_labels ) +
  
  scale_y_discrete(labels = output_labels) + 
  
  labs(title = "Total-effect sensitivity indices (N = 8192)",
       x = NULL,
       y = NULL) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    panel.grid = element_blank(),
    
    plot.title = element_text(hjust = 0.5, face = "bold", size=14),
    
    legend.key.height = unit(1.8, "cm"),
    
    axis.text.x = element_markdown(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10,
      colour = "black"
    ),
    
    axis.text.y = element_markdown(size = 10, colour = "black"),
    
    axis.title = element_text( face = "bold"),
    
    legend.title = element_text( face = "bold"),
    
    plot.background = element_rect(fill = "white",colour = NA),
    
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave(
  "heatmap2.pdf",
  width = 28,
  height = 12,
  units = "cm"
)