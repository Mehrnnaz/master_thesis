library(dplyr)
library(stringr)
library(ggplot2)

pig <- readxl::read_excel("~/OneDrive - University of Helsinki/pigm_int.xlsx")
pig <- pig %>% 
        mutate(`cells.per.mm2` = `Cell Count`/as.numeric(`Area (mm^2)`),
               individual = str_extract(`Image`, "^[^-]+"),
               Treatment = as.factor(Treatment))
head(pig)        
View(pig)

by_fish_pig <- pig %>%
        group_by(Treatment, individual) %>%
        summarise(`Mean Gray Value` = mean(as.numeric(`Mean Gray Value`), na.rm = TRUE),
                  `cells.per.mm2` = mean(`cells.per.mm2`, na.rm = TRUE))
by_fish_pig$Treatment <- as.factor(by_fish_pig$Treatment)

head(by_fish_pig)
View(by_fish_pig)
ggplot(by_fish_pig, aes(x = Treatment, 
                      y = `cells.per.mm2`))+
        geom_boxplot(aes(col = Treatment)) + 
        geom_jitter(width = 0.1)+
        scale_color_manual(values = c("cornflowerblue", "tomato"))+
        theme_bw(base_size = 16)

ggplot(pig, aes(x = Treatment, 
                      y = `cells.per.mm2`))+
        geom_boxplot(aes(col = Treatment)) + 
        geom_jitter(width = 0.1)+
        scale_color_manual(values = c("cornflowerblue", "tomato"))+
        theme_bw(base_size = 16)

