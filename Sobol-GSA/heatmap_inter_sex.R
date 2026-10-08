library(readxl)
library(tidyverse)
library(ggtext)

# ============================================================
# Read data
# ============================================================

df <- read_excel("heatmap_inter_sex_N8192.xlsx")

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
# Improved output labels
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

# ============================================================
# Improved parameter labels
# ============================================================

param_labels <- c(
  "kps"  = "k<sub>ps</sub>",
  "Mmax" = "M<sub>max</sub>",
  "ds"   = "d<sub>s</sub>",
  "k01"  = "k<sub>01</sub>",
  "kpb"  = "k<sub>pb</sub>",
  "ke1"  = "k<sub>e1</sub>",
  "k02"  = "k<sub>02</sub>",
  "kls"  = "k<sub>ls</sub>",
  "klb"  = "k<sub>lb</sub>",
  "ke2"  = "k<sub>e2</sub>",
  "k3"   = "k<sub>3</sub>",
  "k1"   = "k<sub>1</sub>",
  "k2"   = "k<sub>2</sub>"
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
  
  geom_tile() +
  
  scale_fill_gradientn(
    colours = c(
      "#f7fbff",
      "#c6dbef",
      "#6baed6",
      "#2171b5",
      "#08306b"
    ),
    name = NULL
  ) +
  
  scale_x_discrete(
    labels = param_labels
  ) +
  
  scale_y_discrete(
    labels = output_labels
  ) +
  
  labs(
    title = "Parameter interaction effects (N = 8192)",
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    panel.grid = element_blank(),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    
    legend.key.height = unit(1.8, "cm"),
    
    axis.text.x = element_markdown(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10,
      colour = "black"
    ),
    
    axis.text.y = element_markdown(
      size = 10,
      colour = "black"
    ),
    
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    panel.background = element_rect(
      fill = "white",
      colour = NA
    )
  )

# ============================================================
# Save figure
# ============================================================

ggsave(
  "heatmap_interactions.pdf",
  width = 13,
  height = 12,
  units = "cm"
)
