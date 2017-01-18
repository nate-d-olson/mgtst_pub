## 16S PCR amplicon diagram

 
 
library(tidyverse)

# ```{r rra_diag, fig.height = 2, fig.cap = "Diagram of forward and reverse reads relative to 16S rRNA variable region.  Diagonal in forward and reverse reads represent overlap region. Forward (341F), and reverse (806R) primers indicates by vertical bars in variable regions V3 and V5 respecively.", echo = FALSE}
amp_diagram <- function(){
      variable_regions <- data_frame(region_id = rep(paste0("V",c(1:9)),2),
                                     pos = c(8, 95, 305, 486, 745, 884, 
                                             1028, 1179, 1371, 
                                             96, 306, 487, 746, 885, 
                                             1029, 1180, 1372, 1468))
      variable_region_id <- variable_regions %>% 
            group_by(region_id) %>% summarise(pos = mean(pos)) %>% 
            mutate(y = 0.5)
      
      
      primer_df <- data_frame(region_type = "amplicon",
                              region_id = c("341F","806R"),
                              pos_str = c(341, 785),
                              pos_end = c(357, 805)) %>%
            gather("str_end","pos", -region_id, -region_type)
      
      read_df <- frame_data(
            ~region_id, ~str_end, ~pos, ~ymin, ~ymax,
            "F-Read",    "str",  341,     1,     2,
            "F-Read",    "ovr",  505,     1,     2,
            "F-Read",    "end",  641,     1,     1,
            "R-Read",    "str",  505,     2,     2,
            "R-Read",    "ovr",  641,     1,     2,
            "R-Read",    "end",  805,     1,     2
      )
      
      read_id <- read_df %>% filter(str_end != "ovr") %>% 
            group_by(region_id) %>% summarise(pos = mean(pos))
      # ```
      
      
      
      ggplot(variable_regions) + 
            geom_area(data = primer_df, 
                      aes(x = pos, y = 1, fill = region_id)) +
            geom_area(aes(x = pos, y = 1, group = region_id),
                      color = "black", alpha = 0.10) +
            geom_text(data = variable_region_id, 
                      aes(x = pos, y= y, label = region_id)) +
            geom_ribbon(data = read_df, 
                        aes(x = pos, ymin = ymin, ymax = ymax, 
                            fill = region_id)) +
            geom_text(data = read_id, 
                      aes(x = pos, y= 1.5, label = region_id)) +
            theme_void()
}