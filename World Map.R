# ==============================================================================
# FIGURE 2B: GLOBAL DISTRIBUTION MAP OF CLINICAL PRACTICE GUIDELINES
# Script prepared for Nature Communications code availability
# ==============================================================================

# 1. Environment Setup and Package Loading
if(!require(pacman)) install.packages("pacman")
pacman::p_load(
  ggplot2, 
  rnaturalearth, 
  rnaturalearthdata, 
  sf, 
  dplyr, 
  countrycode, 
  readxl,
  extrafont
)

# 2. Data Acquisition
# Note: Update the path to the relative location of your Source Data file
file_path <- "Source_Data.xlsx" 
data <- read_excel(file_path, sheet = "Figure 2b")

# 3. Standardize Country Identifiers
# Convert country names to ISO 3166-1 alpha-3 codes for reliable mapping
data <- data %>%
  mutate(iso_a3 = countrycode(Country, origin = "country.name", destination = "iso3c"))

# 4. Prepare Geospatial Base Map
# Using Robinson projection for better aesthetic balance in academic publishing
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(continent != "Antarctica") # Exclude Antarctica for cleaner visualization

# 5. Data Integration
# Merge guideline counts with geospatial data based on ISO alpha-3 codes
world_data <- left_join(world, data, by = c("iso_a3" = "iso_a3"))

# 6. Data Categorization (Binning)
# Define intervals for guideline counts for discrete color mapping
world_data$Count_group <- cut(
  world_data$Count,
  breaks = c(-Inf, 100, 200, 300, 400, 500, 600, Inf),
  labels = c("0-100", "101-200", "201-300", "301-400", "401-500", "501-600", ">600"),
  include.lowest = TRUE, 
  right = TRUE
)

# 7. Visualization Construction (NC-style Red Gradient Palette)
p <- ggplot(data = world_data) +
  geom_sf(aes(fill = Count_group), color = "white", size = 0.1) + 
  scale_fill_manual(
    # Sequential red palette for quantitative gradient
    values = c(
      "0-100"   = "#fee5d9",
      "101-200" = "#fcbba1",
      "201-300" = "#fc9272",
      "301-400" = "#fb6a4a",
      "401-500" = "#ef3b2c",
      "501-600" = "#cb181d",
      "600"    = "#67000d" 
    ),
    na.value = "#f7f7f7", # Light grey for missing data/non-represented regions
    name = "Number of Guidelines",
    guide = guide_legend(
      direction = "horizontal",
      keyheight = unit(3, units = "mm"),
      keywidth = unit(12, units = "mm"),
      title.position = 'top',
      title.hjust = 0.5,
      label.hjust = 0.5,
      nrow = 1,
      label.position = "bottom"
    )
  ) +
  coord_sf(crs = "+proj=robin") + # Robinson projection
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Global Distribution of Clinical Practice Guidelines")

# 8. Data Validation and Export
# Check for unmatched country names
unmatched_countries <- data %>% filter(is.na(iso_a3))
if(nrow(unmatched_countries) > 0){
  warning("Unrecognized country names found:")
  print(unmatched_countries$Country)
}

# Exporting as PDF (vector format recommended for high-impact journals)
ggsave(
  "Figure2b_Global_Map.pdf", 
  p, 
  width = 10, 
  height = 6, 
  units = "in", 
  device = cairo_pdf
)

# Execution complete
cat("Figure 2b exported successfully to:", getwd())