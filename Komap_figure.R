# =========================================
# Circular Bar Plot (AUC 0.8 to 1.0) 
# Includes per-bar custom colors
# =========================================

# --- Libraries ---
library(ggplot2)
library(dplyr)
library(tidyverse)

# --- Input data ---
# AUROC values for two sites (UPMC, MGB) and two feature sets (Codified, Both)
df <- tribble(
  ~Method,    ~Site,    ~AUROC,
  "Codified", "UPMC",   0.912,
  "Both",     "UPMC",   0.922,
  "Codified", "MGB",    0.974,
  "Both",     "MGB",    0.994
)

# --- Custom colors ---
custom_colors <- c("#a5dff9", "#4ea1d3", "#D1B6E1", "#9B59B6")

# --- Prepare plotting fields ---
df <- df %>%
  mutate(
    label = paste(Site, Method, sep = " - "),  # label shown at the inner radius
    id = row_number(),                         # bar id around the circle
    value_trans = AUROC - 0.8,                 # re-center AUROC to start at 0.8 -> 0 to 0.2
    mycolor = custom_colors                    # assign per-bar color
  )

# --- Radial axis: show ticks from 0.8 to 1.0 while plotting 0.0 to 0.2 internally ---
radial_breaks  <- seq(0, 0.2, by = 0.01)
radial_labels  <- seq(0.8, 1.0, by = 0.01)

# --- Plot ---
p <- ggplot(df, aes(x = factor(id), y = value_trans, fill = mycolor)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y") +
  scale_y_continuous(
    limits = c(0, 0.2),
    breaks = radial_breaks,
    labels = radial_labels
  ) +
  scale_fill_identity() +  # use the hex color values directly
  theme_minimal() +
  theme(
    axis.text.y   = element_blank(),
    axis.title    = element_blank(),
    panel.grid    = element_blank(),
    legend.position = "none"
  ) +
  geom_text(
    aes(x = factor(id), y = 0, label = label),
    size = 1.5,
    hjust = 1,
    color = "black"
  ) +
  ggtitle("Circular Bar Plot: AUC 0.8 to 1.0")

# --- Save and print ---
ggsave("auroc_circular_0.8_to_1.svg", plot = p, width = 3, height = 7, dpi = 600)
p


# =========================================
# Circular Bar Plot with Legend
# AUROC values for multiple methods and sites
# =========================================

# --- 1. Input AUROC data ---
# Each row contains the method, site, and corresponding AUROC
df <- tribble(
  ~Method,        ~Site,   ~AUROC,
  "Main PheCode", "UPMC",  0.854,
  "Main CUI",     "UPMC",  0.798,
  "Codified",     "UPMC",  0.912,
  "Both",         "UPMC",  0.922,
  "Main PheCode", "MGB",   0.942,
  "Main CUI",     "MGB",   0.983,
  "Codified",     "MGB",   0.974,
  "Both",         "MGB",   0.994
)

# --- 2. Define custom colors ---
# Blues for UPMC; pink-red gradient for MGB
custom_colors <- c(
  "#B4D8F9", "#7FB3D5", "#4EA1D3", "#1E91D6",  # UPMC
  "#FADADD", "#F7A8A8", "#F46A6A", "#D62F2F"   # MGB
)

# --- 3. Add plotting labels and transformed AUROC values ---
df <- df %>%
  mutate(
    label = paste(Site, Method, sep = " - "),  # used in legend and text labels
    id = row_number(),                         # order around the circle
    value_trans = AUROC - 0.75                  # shift values so that 0.75 is baseline
  )

# --- 4. Create a named color vector for manual fill mapping ---
named_colors <- setNames(custom_colors, df$label)

# --- 5. Plot ---
# Radial scale: plotting values 0.0–0.25, but labeling as 0.75–1.0
radial_breaks <- seq(0, 0.25, by = 0.01)
radial_labels <- seq(0.75, 1.00, by = 0.01)

p <- ggplot(df, aes(x = factor(id), y = value_trans, fill = label)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y") +
  scale_y_continuous(
    limits = c(0, 0.25),
    breaks = radial_breaks,
    labels = radial_labels
  ) +
  scale_fill_manual(values = named_colors) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Circular Bar Plot: AUROC")

# --- 6. Save and print ---
ggsave("auroc_circular_July_with_legend.svg", plot = p, width = 7, height = 6, dpi = 800)

# =========================================
# Radar Plots (5 rings) for KOMAP Performance at Specificity Thresholds
# Outputs one PDF per threshold (0.90, 0.95, 0.97)
# =========================================

library(fmsb)

# --- Input metric tables per threshold ---
# Rows will be coerced to sites (UPMC, MGB); columns are metrics
# NOTE: Values are in [0.5, 1.0] to match the 5-ring scale and axis labels.

data_list <- list(
  "0.90" = data.frame(
    AUROC       = c(0.854, 0.942),
    AUPRC       = c(0.941, 0.831),
    Sensitivity = c(0.662, 0.716),
    PPV         = c(0.951, 0.732),
    NPV         = c(0.518, 0.907)
  ),
  "0.95" = data.frame(
    AUROC       = c(0.854, 0.942),
    AUPRC       = c(0.941, 0.831),
    Sensitivity = c(0.545, 0.637),
    PPV         = c(0.967, 0.876),
    NPV         = c(0.454, 0.888)
  ),
  "0.97" = data.frame(
    AUROC       = c(0.854, 0.942),
    AUPRC       = c(0.941, 0.831),
    Sensitivity = c(0.502, 0.605),
    PPV         = c(0.976, 0.876),
    NPV         = c(0.436, 0.881)
  )
)

# --- Colors: UPMC (blue), MGB (red) ---
colors <- c("#1f77b4", "#d62728")

# --- Loop over thresholds and export one PDF each ---
for (threshold in names(data_list)) {
  
  df <- data_list[[threshold]]
  
  # Reorder columns for a consistent clockwise metric sequence on the radar
  df <- df[, c("AUROC", "AUPRC", "NPV", "PPV", "Sensitivity")]
  
  # Site names as row labels
  rownames(df) <- c("UPMC", "MGB")
  
  # fmsb::radarchart expects first two rows to be the max and min bounds
  df_radar <- rbind(
    rep(1, 5),     # max values across all axes
    rep(0.5, 5),   # min values across all axes
    df             # actual site values
  )
  
  # File name per threshold (PDF)
  svg_filename <- paste0("radar_KOMAP_specificity_", threshold, ".pdf")
  
  # Open device (7x7 inches)
  pdf(svg_filename, width = 7, height = 7)
  
  # --- Radar chart ---
  radarchart(
    df_radar,
    seg        = 5,                 # five rings between min and max
    axistype   = 1,
    pcol       = colors,            # line colors per site
    pfcol      = NA,                # no fill (transparent)
    plwd       = 2,                 # line width
    plty       = 1,                 # line type
    cglcol     = "grey",            # grid line color
    cglty      = 1,                 # grid line type
    cglwd      = 0.8,               # grid line width
    axislabcol = "black",           # axis label color
    caxislabels = seq(0.5, 1.0, 0.1),  # tick labels shown on the radial axes
    vlcex      = 1,                 # vertex label size
    title      = paste("KOMAP Performance @ Specificity", threshold)
  )
  
  # Legend
  legend(
    "bottomright",
    legend = c("UPMC", "MGB"),
    col    = colors,
    lty    = 1,
    lwd    = 2,
    bty    = "n",
    cex    = 0.9
  )
  
  dev.off()
}


# =========================================
# Radar Plots (Supplementary Figure)
# 4 series: UPMC Codified, UPMC Both, MGB Codified, MGB Both
# Exports one SVG per specificity threshold (0.970, 0.950, 0.900)
# =========================================

library(fmsb)

# --- Metrics for three specificity thresholds ---
# Rows (to be labeled later): UPMC Codified, UPMC Both, MGB Codified, MGB Both
# Columns (metrics): AUROC, AUPRC, PPV, NPV, Sensitivity
data_list <- list(
  "0.970" = data.frame(
    AUROC       = c(0.912, 0.922, 0.974, 0.994),
    AUPRC       = c(0.963, 0.966, 0.879, 0.940),
    PPV         = c(0.981, 0.982, 0.868, 0.912),
    NPV         = c(0.510, 0.541, 0.884, 0.986),
    Sensitivity = c(0.631, 0.674, 0.609, 0.957)
  ),
  "0.950" = data.frame(
    AUROC       = c(0.912, 0.922, 0.974, 0.994),
    AUPRC       = c(0.963, 0.966, 0.879, 0.940),
    PPV         = c(0.975, 0.974, 0.808, 0.861),
    NPV         = c(0.567, 0.557, 0.892, 0.985),
    Sensitivity = c(0.711, 0.698, 0.645, 0.957)
  ),
  "0.900" = data.frame(
    AUROC       = c(0.912, 0.922, 0.974, 0.994),
    AUPRC       = c(0.963, 0.966, 0.879, 0.940),
    PPV         = c(0.952, 0.953, 0.764, 0.764),
    NPV         = c(0.630, 0.634, 1.000, 1.000),
    Sensitivity = c(0.791, 0.794, 1.000, 1.000)
  )
)

# --- Series colors and labels ---
# UPMC (light/dark blue), MGB (light/dark red)
colors <- c("#a5dff9", "#4ea1d3", "#e85a71", "#a6172d")
labels <- c("UPMC Codified", "UPMC Both", "MGB Codified", "MGB Both")

# --- Loop over thresholds and export SVGs ---
for (threshold in names(data_list)) {
  
  df <- data_list[[threshold]]
  rownames(df) <- labels
  
  # Ensure a consistent clockwise order of metrics on the radar
  df <- df[, c("AUROC", "AUPRC", "PPV", "NPV", "Sensitivity")]
  
  # fmsb::radarchart requires the first two rows as max and min bounds
  # Here we cap the scale between 0.6 and 1.0 to emphasize high-performance range
  df_radar <- rbind(
    rep(1.0, 5),   # max values across all axes
    rep(0.6, 5),   # min values across all axes
    df             # actual series values
  )
  
  # Open SVG device (7x7 inches)
  svg_filename <- paste0("radar_KOMAP_4lines_", threshold, ".svg")
  svg(svg_filename, width = 7, height = 7)
  
  # Draw radar chart
  radarchart(
    df_radar,
    axistype    = 1,
    pcol        = colors,           # line colors for the 4 series
    pfcol       = NA,               # transparent fill
    plwd        = 2,                # line width
    plty        = 1,                # line type
    cglcol      = "grey",           # grid color
    cglty       = 1,                # grid line type
    cglwd       = 0.8,              # grid line width
    axislabcol  = "black",          # axis tick label color
    caxislabels = seq(0.6, 1.0, 0.1),
    vlcex       = 1,                # vertex label size
    title       = paste("KOMAP Performance @ Specificity", threshold)
  )
  
  # Legend
  legend(
    "bottomright",
    legend = labels,
    col    = colors,
    lty    = 1,
    lwd    = 2,
    bty    = "n",
    cex    = 0.8
  )
  
  # Close device
  dev.off()
}

