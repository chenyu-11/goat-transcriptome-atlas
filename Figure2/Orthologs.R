library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)

df <- read_excel("Orthologs.xlsx")
df <- df %>% filter(Species != "Goat")

# Order species by evolutionary distance (from closest to farthest)
species_order <- c("Sheep", "Cattle", "Dog", "Horse", "Human", "Mouse", "Chicken", "Fly", "Worm")

df_long <- df %>%
  select(Species, ends_with("(%)")) %>%
  pivot_longer(cols = -Species, names_to = "Category", values_to = "Percentage") %>%
  mutate(Category = gsub("_orthologs \\(%)", "", Category)) %>%
  mutate(Category = factor(Category, levels = c("HKG", "TSG", "HUBG", "allgenes"))) %>%
  mutate(Species = factor(Species, levels = species_order))

#=========================================================Dot plot

# Shapes and colors
shape_map <- c(HKG = 21, TSG = 22, HUBG = 24, allgenes = 25)
#shape_map <- c(HKG = 21, TSG = 21, HUBG = 21, allgenes = 21)
fill_colors <- c("HKG" = "#E41A1C", "TSG" = "#377EB8", 
                 "HUBG" = "#4DAF4A", "allgenes" = "#FF7F00")

# Plot
p <- ggplot(df_long, aes(x = Species, y = Percentage, fill = Category, shape = Category)) +
  geom_point(size = 3.5, color = "black", stroke = 0.8) + # Ensure the fill attribute is used
  scale_shape_manual(values = shape_map, name = "Category") + # Add name parameter to unify legend title
  scale_fill_manual(name = "Category", values = fill_colors) + # Ensure legend colors match point colors
  labs(
    title = "Proportion of Goat Orthologous Genes in Different Species (Ordered by Evolutionary Distance)",
    x = "Species (from closest to farthest)",
    y = "Proportion of orthologous genes (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    legend.position = "right"
  )

print(p)


#=========================================================Line plot
# Color definitions (used for both fill and color)
colors <- c("HKG" = "#E41A1C", "TSG" = "#377EB8", 
            "HUBG" = "#4DAF4A", "allgenes" = "#FF7F00")
# Shapes: 21-25 are fillable
shapes <- c(HKG = 21, TSG = 22, HUBG = 24, allgenes = 25)

# Plot
p <- ggplot(df_long, aes(x = Species, y = Percentage)) +
  # Draw lines first (below points)
  geom_line(aes(group = Category, color = Category), linewidth = 1) +
  # Then draw points (with fill)
  geom_point(aes(fill = Category, shape = Category), 
             size = 3.5, color = "black", stroke = 0.8) +
  # Manually set colors, shapes, and fills
  scale_color_manual(name = "Category", values = colors) +
  scale_fill_manual(name = "Category", values = colors) +
  scale_shape_manual(name = "Category", values = shapes) +
  labs(
    title = "Proportion of Goat Orthologous Genes in Different Species (Ordered by Evolutionary Distance)",
    x = "Species (from closest to farthest)",
    y = "Proportion of orthologous genes (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    legend.position = "right"
  )

print(p)