# load relevant libraries and source temp_bar_data function (in utilities.R)
setwd("~/Documents/Thesis")
library(StereoMorph)
library(dplyr)
source("utilities.R")
library(ggplot2)


# load landmark files
shapes <- readShapes("~/Documents/Thesis/bar_landmarks_top")

# make a dataframe with 1 row per landmark
for (i in 1:dim(shapes$landmarks.pixel)[3]) {
  
  # get landmarks for that image
  landmarks <- shapes$landmarks.pixel[ , , i]
  
  # get dimension names (landmark names)
  nm <- dimnames(shapes$landmarks.pixel)[[3]][i]
  
  
  # match to the right curves object (I had an issue with this before so this is an extra step)
  idx <- match(nm, names(shapes$curves.pixel))
  if (is.na(idx)) { next }
  
  # get the dorsal curve out
  dors_ref <- shapes$curves.pixel[[idx]]$DOR_dorsal
  
  # if they're all NA values, skip this image
  if (sum(is.na(landmarks)) == prod(dim(landmarks))) { next }
  
  # otherwise, get lateral image info
  bar_data <- temp_bar_data(landmarks, dors_ref)
  
  # image name
  bar_data$image <- gsub("-top", "", dimnames(shapes$landmarks.pixel)[[3]][i])
  offsets <- calculate_offset(bar_data)
  bar_data$offset_sameside <- offsets$sameside_offset
  bar_data$offset_diffside <- offsets$diffside_offset
  
  # if it's the first iteration, start the big dataframe
  # otherwise append it
  if (i == 1) {
    by_bar <- bar_data
  } else {
    by_bar <- rbind(by_bar, bar_data)
  }
}
head(by_bar)

by_bar <- by_bar %>%
        mutate(
                temp = case_when(
                        grepl("23", image) ~ "23",
                        grepl("29", image) ~ "29",
                        TRUE ~ "28"
                ),
                treatment = case_when(
                        grepl("23", image) ~ "cool",
                        grepl("29", image) ~ "warm",
                        grepl("28", image) ~ "warm",
                ),
                
                bar_location = case_when(
                        grepl("DORS", lmk_name) ~ "dorsal",
                        grepl("LAT", lmk_name) ~ "lateral",
                        TRUE ~ NA_character_
                ),
                
                brood = case_when(
                        grepl("BR03", image) ~ "BR03",
                        grepl("BR04", image) ~ "BR04",
                        grepl("BR05", image) ~ "BR05",
                        grepl("BR06", image) ~ "BR06",
                        TRUE ~ "BR07"
                )
        )


# make a dataframe with per-fish values (e.g. overlap, average offset, barcount)
# 1 row = 1 individual
for (i in 1:length(unique(by_bar$image))) { 
  image_subset <- by_bar[by_bar$image == unique(by_bar$image)[i], ]
  image_subset$bar_number <- as.numeric(image_subset$bar_number)
  
  for (side in c("LEFT", "RIGHT")) {
    for (dv in c("DORS", "LAT")) {
      x <- table(image_subset$bar_number[image_subset$dors_vent == dv & image_subset$side == side])
      if( any(x != 2) ) {warning(paste(image_subset$image[1], "has an ant/post mismatch"))}
    }
  }
  
  # 1. number of bars, dorsal/ventral, left/right
  # in order: lhs_dor, lhs_mid, rhs_dor, rhs_mid
  barcounts <- as.data.frame(summarise(group_by(image_subset, side, dors_vent),
                                       max(bar_number), 
                                       .groups = "drop_last"))
  barcount_lhs_dor <- barcounts$`max(bar_number)`[1]
  barcount_lhs_mid <- barcounts$`max(bar_number)`[2]
  barcount_rhs_dor <- barcounts$`max(bar_number)`[3]
  barcount_rhs_mid <- barcounts$`max(bar_number)`[4]
  
  # 2. lateral overlap (dorsal/ventral)
  overlap_dor <- diffside_overlap(image_subset[which(image_subset$dors_vent == "DORS"), ])
  overlap_mid <- diffside_overlap(image_subset[which(image_subset$dors_vent == "LAT"), ])
  overlap_lhs <- sameside_overlap(image_subset[which(image_subset$side == "LEFT"), ])
  overlap_rhs <- sameside_overlap(image_subset[which(image_subset$side == "RIGHT"), ])
  
  
  # make a new row
  # dataframe format:
  # individual | lhs_dors_bars | rhs_dors_bars | lhs_vent_bars etc
  # 1 row = 1 fish
  new_row <- data.frame(id = image_subset$image[1],
                        barcount_lhs_dor = barcount_lhs_dor,
                        barcount_lhs_mid = barcount_lhs_mid,
                        barcount_rhs_dor = barcount_rhs_dor,
                        barcount_rhs_mid = barcount_rhs_mid,
                        overlap_dor = overlap_dor,
                        overlap_mid = overlap_mid,
                        overlap_lhs = overlap_lhs,
                        overlap_rhs = overlap_rhs,
                        offset_mean_sameside = mean(image_subset$offset_sameside),
                        offset_mean_diffside = mean(image_subset$offset_diffside))
  if (i == 1) {
    by_fish <- new_row
  } else {
    by_fish <- rbind(by_fish, new_row)
  }
 
}
head(by_fish)

# there are a LOT of overlap values, so let me break it down:

# abs_overlap = the absolute overlap (in units of standard length), e.g. how
# much of the fish is covered by overlapping bars

# rel_overlap = abs_overlap divided by the coverage of the smaller pattern
# (because the maximum overlap is the proportion of standard length covered by
# the smaller pattern)

# overlap_lhs / overlap_rhs: overlap between dorsal and midline bars on the left
# or righthand sides, respectively

# overlap_dor / overlap_mid: overlap between left and right bar patterns at the
# same height (dorsal or midline)

# plotting ####
library(ggplot2)

# we can look at offset as a function of standard length like before:
ggplot(by_bar, aes(x = location_sl, 
                   y = offset_diffside, 
                   group = treatment)) +
        geom_point(aes(color = treatment), alpha = 0.2) + 
        geom_smooth(aes(color = treatment), se = FALSE) +
        theme_classic(base_size = 16) +
        scale_color_manual(values = c("29" = "tomato", "28" = "goldenrod1", "23"="cornflowerblue"))


# for same side offset:
ggplot(by_bar, aes(x = location_sl, 
                   y = offset_sameside)) +
  geom_point(aes(color = temp), alpha = 0.2) + geom_smooth(aes(color = temp), se = FALSE) +
  theme_classic(base_size = 16) +
  scale_color_manual(values = c("29" = "tomato", "28" = "goldenrod1", "23"="cornflowerblue"))


# we can also use the ID string to get the treatment info:
img_info <- stringr::str_remove_all(by_fish$id, "_ch00|zoomed|left")
img_info <- stringr::str_split(img_info, "-", simplify = TRUE)
img_info <- as_tibble(img_info[, 1:4])

colnames(img_info) <- c("parents","temp", "brood", "individual")
by_fish <- cbind(by_fish, img_info)
by_fish$temp[by_fish$temp == "29" | by_fish$temp == "28"] <- "warm"
by_fish$temp[by_fish$temp == "23"] <- "cool"

# you'll probably plot by temperature and/or brood instead, but the example
# images are all from the same treatment, so I'm plotting by individual:
ggplot(by_fish, aes(x = temp, y = overlap_dor.rel_overlap)) +
  geom_jitter(height = 0, width = 0.05, aes(color=temp)) +
scale_color_manual(values = c("cool" = "cornflowerblue", "warm" = "tomato")) + 
        ylim(c(0.5, 1))+
        theme_bw(base_size = 16)

# ...ok, I stop now, this is boring with only three images...
ggplot(by_fish, aes(x = overlap_mid.rel_overlap, 
                    y = overlap_dor.rel_overlap)) +
  geom_point(aes(color=temp))+
        scale_color_manual(values = c("cool" = "cornflowerblue", "warm" = "tomato"))+
        theme_bw(base_size = 16)



ggplot(by_fish, aes(x = temp, 
                    y = overlap_mid.abs_overlap)) +
  geom_point(aes(color=temp))+
        geom_point(aes())+
        scale_color_manual(values = c("cool" = "cornflowerblue", "warm" = "tomato"))+
        geom_jitter(height = 0, width = 0.05, aes(color=temp), alpha = 0.2)+
        theme_bw(base_size = 16)


ggplot(by_fish, aes(x = barcount_rhs_dor,
                    y = barcount_lhs_dor)) +
        geom_point(aes(color=temp), alpha = 0.3)+
        geom_jitter(height = 0.05, width = 0.05, aes(color=temp), alpha = 0.3)+
        theme_bw(base_size = 16)+
        scale_color_manual(values = c("cool" = "cornflowerblue", "warm" = "tomato"))
        

write.csv(by_fish, file = paste0("~/Documents/Thesis/Data/by_fish", Sys.Date(), ".csv"))

write.csv(by_bar, file = paste0("~/Documents/Thesis/Data/by_bar", Sys.Date(), ".csv"))



head(by_fish)

