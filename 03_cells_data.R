library(dplyr)
library(readr)
library(ggplot2)
library(stringr)

data_dir_23 <- '~/OneDrive - University of Helsinki/Thesis_Images/Phenotype/PigmInt-CR04/23'
csv_files_23 <- list.files(data_dir_23, pattern = "\\.csv$", full.names = TRUE)
data_list_23 <- lapply(csv_files_23, function(file) {
        df <- read_csv(file, show_col_types = FALSE)  
        df <- df %>% mutate(filename = rep(basename(file), n())) 
        return(df)
})


data_dir_29 <- '~/OneDrive - University of Helsinki/Thesis_Images/Phenotype/PigmInt-CR04/29'
csv_files_29 <- list.files(data_dir_29, pattern = "\\.csv$", full.names = TRUE)
data_list_29 <- lapply(csv_files_29, function(file) {
        df <- read_csv(file, show_col_types = FALSE)  
        df <- df %>% mutate(filename = rep(basename(file), n()))  
        
        return(df)
})




combined_data <- bind_rows(data_list_23, data_list_29)  %>% 
                mutate(Treatment = ifelse(grepl("23", filename), "23", "29"))  %>% 
                select(-FeretX, -FeretY, -FeretAngle, -MinFeret)
combined_data$filename <- gsub(".csv", "", combined_data$filename)
head(combined_data)
View(combined_data)

intbars <- combined_data %>%
        group_by(filename) %>%
        summarise(mean_area = mean(Area, na.rm = TRUE),
                  mean_gv = mean(Mean, na.rm = TRUE)) %>%
        mutate(temp = ifelse(grepl("23", filename), "23", "29"),
               individual = str_extract(filename, "^[^-]+"))
head(intbars)




by_fish_cells <- intbars %>%
        group_by(temp, individual) %>%
        summarise(
        mean_area_cells = mean(mean_area, na.rm = TRUE),
        mean_gv_cells = mean(mean_gv, na.rm = TRUE))
               
head(by_fish_cells)               
View(by_fish_cells)            

ggplot(combined_data, aes(x = Treatment, y = Area)) +
        geom_boxplot(aes(col = Treatment)) +
        geom_jitter(width = 0.1)+
        theme_bw()+
        scale_color_manual(values = c("cornflowerblue", "tomato"))

ggplot(combined_data, aes(x = Treatment, y = Mean)) +
        geom_boxplot(aes(col = Treatment)) +
        geom_jitter(width = 0.1)+
        theme_bw()+
        scale_color_manual(values = c("cornflowerblue", "tomato"))


ggplot(new_df, aes(x = temp , y = mean_area)) +
        geom_boxplot(aes(col = temp)) +
        geom_jitter(width = 0.1)+
        scale_color_manual(values = c("cornflowerblue", "tomato"))+
        theme_bw()

ggplot(by_fish_pig, aes(x = individual , y = mean_gv)) +
#geom_boxplot(data = new_df, aes(col = temp), size = 0.75) +
        geom_jitter(width = 0.1, , size = 2)+
        facet_wrap(~temp)+

        geom_boxplot(data = intbars, aes(col = temp), size = 0.75) +
        geom_jitter(data = intbars, aes(x = individual, y = mean_gv, col = temp), width = 0.1, size = 2, alpha = 0.7)+
        facet_wrap(~temp)+
        scale_color_manual(values = c("cornflowerblue", "tomato"))+
        
        theme_bw(base_size = 16)+
        theme(legend.position = "none")
        
        
write.csv(by_fish_pig, file = paste0("~/Documents/Thesis/Data/by_fish_pig", Sys.Date(), ".csv"))
write.csv(intbars, file = paste0("~/Documents/Thesis/Data/intbars", Sys.Date(), ".csv"))
        
write.csv(combined_data, file = paste0("~/Documents/Thesis/Data/combined_data", Sys.Date(), ".csv")) 





