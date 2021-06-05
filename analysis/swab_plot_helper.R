# plotting
library(ggplot2)
library(ggbeeswarm) # <- geom_quasirandom

# stats
library(MASS) # <- glm.nb
library(speedglm)

# tidyverse
library(furrr) # <- parallel map (future_map, plan) (devtools for walk)
library(readxl) # <- read_xlsx
library(magrittr)
library(tidyverse)

select = dplyr::select #, MASS::select masks dplyr...

# ------------------------------------------------------------------------------------
# style plots

theme_pub <- function(base_size = 11, base_family = "") {
  # based on https://github.com/noamross/noamtools/blob/master/R/theme_nr.R
  # start with theme_bw and modify from there!
  theme_bw(base_size = base_size, base_family = base_family) +# %+replace%
    theme(
      # grid lines
      panel.grid.major.x = element_line(colour="#ECECEC", size=0.5, linetype=1),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(colour="#ECECEC", size=0.5, linetype=1),
      panel.background   = element_blank(),
      
      # axis options
      axis.ticks.y   = element_blank(),
      axis.title.x   = element_text(size=rel(2), vjust=0.25),
      axis.title.y   = element_text(size=rel(2), vjust=0.35),
      axis.text      = element_text(color="black", size=rel(1)),
      
      # legend options
      legend.title    = element_text(size=rel(1.5)),
      legend.key      = element_rect(fill="white"),
      legend.key.size = unit(1, "cm"),
      legend.text     = element_text(size=rel(1.5)),
      
      # facet options
      strip.text = element_text(size=rel(2)),
      strip.background = element_blank(),
      
      # title options
      plot.title = element_text(size=rel(2.25), vjust=0.25, hjust=0.5)
    )
}
theme_set(theme_pub(base_size=8))

