
# temp_bar_data ###
# makes a dataframe where 1 row = 1 landmark
temp_bar_data <- function(landmarks, dors_ref) {
        
        # get bar landmarks
        bar_lmk <- na.omit(landmarks[grep("_(RIGHT|LEFT)", rownames(landmarks)), ])
        
        # get midline length
        sl_px <- nrow(dors_ref)
        
        # get barname info
        barname_data <- stringr::str_split(rownames(bar_lmk), "_", simplify = TRUE)
        barname_data <- cbind(barname_data, stringr::str_extract(barname_data[, 2], "[0-9]"))
        barname_data[ , 2] <- gsub("[0-9]", "", barname_data[ , 2])
        
        # for every landmark...
        for (b in 1:nrow(bar_lmk)) {
                # print(b)
                
                # extract bar coords
                bar <- bar_lmk[b, ] 
                
                # find the closest reference point
                min_ref <- which.min(apply(dors_ref, 1, \(x) euc_dist(bar, x)))
                
                # make a new row:
                lmk_name <- rownames(bar_lmk)[b]
                new_row <- data.frame(lmk_name = lmk_name,
                                      view = "dorsal",
                                      side = barname_data[b, 2],
                                      dors_vent = barname_data[b, 1],
                                      ant_post = barname_data[b, 3],
                                      bar_number = barname_data[b, 4],
                                      location_px = min_ref,
                                      location_sl = min_ref / sl_px)
                
                if (b == 1) {
                        bar_info <- new_row
                } else {
                        bar_info <- rbind(bar_info, new_row)
                }
        }
        
        return(bar_info)
}

# euclidean distance ####
euc_dist <- function(p1, p2) {
        sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)
}
# overlap calculations ####
diffside_overlap <- function(fish_bars) {
        
        # make dataframe with just location, anterior/posterior
        events <- data.frame(location = fish_bars$location_sl,
                             ant_post = fish_bars$ant_post,
                             side = fish_bars$side)
        
        # convert to numeric - on (anterior = start of bar) is 1, off is -1
        events$state <- (1 - as.numeric(factor(events$ant_post))) * 2 + 1
        
        # order by location
        events <- events[order(events$location), ]
        
        
        if (any(is.na(events$ant_post))) stop("NA values detected in ant_post column")
        if (any(is.na(events$location))) stop("NA values detected in location column")
        
        
        # tracking variables for overlap
        current_overlap <- 0
        total_overlap_length <- 0
        last_coordinate <- events$location[1]
        
        # for every event...
        for (j in 1:nrow(events)) {
                #print(j)
                # if both vectors are 'on' (+1)
                if (current_overlap == 2) {
                        total_overlap_length <- total_overlap_length + (events$location[j] - last_coordinate)
                }
                
                # current overlap = last state plus new state
                # on -> on = 1 + 1 = 2
                # off -> on = -1 + 1 = 0
                # off -> off = -1 + -1 = -2
                # on -> off = 1 + -1 = 0
                current_overlap <- current_overlap + events$state[j]
                last_coordinate <- events$location[j]
                if (any(is.na(events$ant_post))) stop("NA values detected in ant_post column")
                if (any(is.na(events$location))) stop("NA values detected in location column")
                if (any(is.na(current_overlap))) stop("NA values detected in location column")
        }
        
        # theoretical max overlap = amount of time that the shorter vector spends in 'on' state
        lhs <- events[which(events$side == "LEFT"), ]
        if (table(lhs$ant_post)[1] != table(lhs$ant_post)[2]) {
                stop("Ant/post mismatch on LHS")
        }
        lsum <- sum(lhs$location[which(lhs$ant_post == "post")] - 
                            lhs$location[which(lhs$ant_post == "ant")])
        
        rhs <- events[which(events$side == "RIGHT"), ]
        rsum <- sum(rhs$location[which(rhs$ant_post == "post")] - 
                            rhs$location[which(rhs$ant_post == "ant")])
        max_overlap <- min(c(lsum, rsum))
        if (nrow(lhs) == 0 || nrow(rhs) == 0) stop("Empty subset detected for LEFT or RIGHT bars")
        return(list(abs_overlap = total_overlap_length,
                    rel_overlap = total_overlap_length / max_overlap))
}
sameside_overlap <- function(fish_bars) {
        
        # make dataframe with just location, anterior/posterior
        events <- data.frame(location = fish_bars$location_sl,
                             ant_post = fish_bars$ant_post,
                             side = fish_bars$dors_vent)
        
        # convert to numeric - on (anterior = start of bar) is 1, off is -1
        events$state <- (1 - as.numeric(factor(events$ant_post))) * 2 + 1
        
        # order by location
        events <- events[order(events$location), ]
        
        # tracking variables for overlap
        current_overlap1 <- 0
        total_overlap_length <- 0
        last_coordinate <- events$location[1]
        
        # for every event...
        for (j in 1:nrow(events)) {
                
                # if both vectors are 'on' (+1)
                if (current_overlap1 == 2) {
                        total_overlap_length <- total_overlap_length + (events$location[j] - last_coordinate)
                }
                
                # current overlap = last state plus new state
                # on -> on = 1 + 1 = 2
                # off -> on = -1 + 1 = 0
                # off -> off = -1 + -1 = -2
                # on -> off = 1 + -1 = 0
                current_overlap1 <- current_overlap1 + events$state[j]
                last_coordinate <- events$location[j]
        }
        
        # theoretical max overlap = amount of time that the shorter vector spends in 'on' state
        dor <- events[which(events$side == "DORS"), ]
        if (table(dor$ant_post)[1] != table(dor$ant_post)[2]) {
                stop("Ant/post mismatch on DORS")
        }
        dsum <- sum(dor$location[which(dor$ant_post == "post")] - 
                            dor$location[which(dor$ant_post == "ant")])
        
        mid <- events[which(events$side == "LAT"), ]
        msum <- sum(mid$location[which(mid$ant_post == "post")] - 
                            mid$location[which(mid$ant_post == "ant")])
        max_overlap <- min(c(dsum, msum))
        
        return(list(abs_overlap = total_overlap_length,
                    rel_overlap = total_overlap_length / max_overlap))
}
# offset calculations ####
calculate_offset <- function(fish_bars) {
        
        sameside_offset <- c()
        diffside_offset <- c()
        
        for (b in 1:nrow(fish_bars)) {
                bar <- fish_bars[b, ]
                
                # compare to bars at the same height, but on the other side:
                diffside_subset <- dplyr::filter(fish_bars, 
                                                 side != bar$side, 
                                                 dors_vent == bar$dors_vent,
                                                 ant_post == bar$ant_post)
                
                # compare to bars on a different plane, but on the same side:
                sameside_subset <- dplyr::filter(fish_bars, 
                                                 side == bar$side, 
                                                 dors_vent != bar$dors_vent,
                                                 ant_post == bar$ant_post)
                min_diffside <- min(abs(bar$location_sl - diffside_subset$location_sl))
                min_sameside <- min(abs(bar$location_sl - sameside_subset$location_sl))
                diffside_offset <- c(diffside_offset, min_diffside)
                sameside_offset <- c(sameside_offset, min_sameside)
        }
        
        return(data.frame(sameside_offset = sameside_offset, diffside_offset = diffside_offset))
}

