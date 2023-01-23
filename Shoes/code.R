## ----setup, include = F-----------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = F, message = F, warning = F, dpi = 300)


## ----initial----------------------------------------------------------------------------------------------------
library(tidyverse)
library(imager)

im <- load.dir("orig-images/", "2018.*6_1_1.*")
names(im) <- str_replace(names(im), "\\d{6}R_(\\d{8})_.*", "\\1")
im <- map_il(im, grayscale)
plot(im, "row")


## ----zoom-in----------------------------------------------------------------------------------------------------
# plot(as.cimg(im[[1]][850:1050,2100:2300]))
plot(as.cimg(im[[1]][380:580,2100:2300]))


## ----quantize---------------------------------------------------------------------------------------------------
#' Limit the number of colors in the image
#' 
#' @param shoe cimage
#' @param n number of colors in returned image
#' @param return cimg
#' @import imager
#' @import dplyr
#' @importFrom tidyr gather
#' @export
quantize_colors <- function(shoe, n = 16) {
  if (max(shoe) > 255) {
    shoe <- renorm(shoe)
  }
  # https://lumiamitie.github.io/r/imager-color-quantization/
  
  shoe_df <- shoe %>%
    as.data.frame(wide = 'c') %>%
    tbl_df()
  
  shoe_cluster <- suppressMessages(
    shoe_df %>% select(-x, -y) %>% kmeans(centers = n)
  )
  
  shoe_df %>%
    mutate(label = as.character(shoe_cluster$cluster)) %>%
    select(x, y, label) %>%
    left_join(
      shoe_cluster$centers %>% 
        tbl_df %>% 
        mutate(label = as.character(row_number())), 
      by = "label") %>%
    select(-label) %>%
    tidyr::gather(key = 'cc', value = 'value', starts_with('c.')) %>%
    mutate(cc = as.integer(gsub('c\\.', '', cc))) %>%
    as.cimg(dim = dim(shoe))
}

shoe <- map_il(im, quantize_colors, 8)
plot(shoe, "row")


## ---------------------------------------------------------------------------------------------------------------
#' Removes a label from the shoeprint scan
#' 
#' This algorithm assumes the label will be larger than 50 px and will be dark
#' 
#' @param shoe cimg image type
#' @param ... arguments to pass to imager::threshold
#' @return cimg with any label which may be present removed
#' @import imager
#' @export
remove_print_label <- function(shoe, ...) {
  
  if (max(shoe) > 255) {
    shoe <- renorm(shoe)
  }
  
  if (!exists("thr")) {
    thr <- "10%"
  }
  replace_color <- max(shoe) 
  
  # This identifies the label if it is bigger than 50 px and dark
  z <- threshold(shoe, thr = thr) %>%
    grow(100) %>%
    shrink(150)
  
  zlab <- label(z) # Label continuous areas
  chunks <- table(zlab) %>% sort() %>% rev() # Get frequencies/areas in correct order
  chunks <- chunks[-1] # Remove largest chunk - background
  
  largest_region <- zlab == as.numeric(names(chunks)[1])
  
  # Only remove label if it has a reasonable area
  islabeled <- mean(!largest_region) < .98
  
  if (islabeled) {
    shoe[largest_region] <- replace_color
  }
  
  return(shoe)
}

shoe_nolab <- map_il(shoe, remove_print_label)
plot(shoe_nolab, "row")


## ---------------------------------------------------------------------------------------------------------------
#' Remove local background 
#' 
#' @param shoe cimage
#' @param n number of sub-images to use. shoe is divided into an n x n grid for local background computation
#' @param threshold if a pixel is within threshold of the median pixel value, it should be set to white as well
#' @return a cimg
#' @export
#' @import imager
#' @importFrom purrr map_df
#' @importFrom dplyr data_frame
#' @importFrom purrr map
remove_local_background <- function(shoe, n = 10, threshold = 10, borderarea = 100, borderonly = F) {
  grid_shoe <- tibble(
    x = 1:n,
    xstr = imsplit(renorm(shoe), axis = "x", nb = n)
  ) %>%
    mutate(
      ydat = purrr::map(xstr, ~tibble(y = 1:n, xystr = as.list(imsplit(., "y", nb = n))))
    ) %>%
    unnest("ydat") %>%
    mutate(
      medvalue = map_dbl(xystr, median),
      normshoe = map2(xystr, medvalue, function(tmp, med) {
        tmp2 <- tmp
        tmp2[tmp >= med] <- 255
        as.cimg(tmp2)
      })
    )
    #     split(.$x) %>%
    # map_df(.x = ., .f = function(zz) {
    #   tmp <- imsplit(zz$xstr[[1]], axis = "y", nb = n)
    #   data_frame(
    #     x = zz$x,
    #     y = 1:n,
    #     xyshoe = tmp,
    #     # Get local background values
    #     medvalue = sapply(tmp, median),
    #     # Replace local background with white
    #     normshoe = lapply(tmp, function(x) {y <- x; y[(x >= median(x))] <- 255; y[abs(y - median(x)) < threshold] <- 255; y})
    #   )
    # })
  
  center_shoe <- !px.borders(shoe, n = borderarea)
  
  shoe_reassemble <- select(grid_shoe, x, y, normshoe) %>%
    split(.$x) %>%
    map(~imappend(imlist = .$normshoe, axis = "y")) %>%
    imappend(imlist = ., axis = "x") %>%
    renorm()
  
  if (borderonly) {
    shoe_fixpieces <- renorm(shoe)
    shoe_fixpieces[!center_shoe] <- shoe_reassemble[!center_shoe]
  } else {
    shoe_fixpieces <- shoe_reassemble
  }

  shoe_fixpieces %>% renorm()
}

shoe_lbg <- map_il(shoe_nolab, remove_local_background, n=20)
plot(shoe_lbg, "row")


## ---------------------------------------------------------------------------------------------------------------
#' Crop shoe border
#' 
#' Remove most of the white space around the shoe. This function splits the image
#' into two pieces at the middle, then uses the mean intensity (0-255) in the chosen dimension
#' to determine where the image should be cropped. The region to crop is the value closest to the 
#' center of the image which is within tol of the maximum mean intensity found in the
#' image. 
#' @param shoe cimg
#' @param axis either "x", "y", or "xy"
#' @param sigma radius to use for blurring the image
#' @param tol change the tolerance of the cropping function
#' @export
#' @import imager
crop_border <- function(shoe, axis = "xy", sigma = 5, tol = 0.04) {
  stopifnot(axis %in% c("x", "y", "xy"))
  
  if (axis == "xy") {
    axis = c("x", "y")
  }
  
  modimg <- shoe
  
  # Deal with x first
  if ("x" %in% axis) {
    tshoe_x <- imsplit(modimg, "x", 2)
    tshoe_x[[2]] <- mirror(tshoe_x[[2]], "x")
    
    tshoe_xfix <- lapply(tshoe_x, function(xx) {
      xxmod <- isoblur(xx, sigma)
      yy <- apply(xxmod, 1, mean) 
      zz <- which(abs(yy - max(yy)) < tol) %>% max()
      
      # Don't crop if removing more than 25% of the image
      if (zz > .5*length(yy)) {
        zz <-  -1
      }
      imsub(xx, x > (zz + 1)) 
    }) 
    
    tshoe_xfix[[2]] <- mirror(tshoe_xfix[[2]], "x")
    modimg <- imappend(tshoe_xfix, axis = "x")
    rm(tshoe_x, tshoe_xfix)
  }
  
  if ("y" %in% axis) {
    tshoe_y <- imsplit(modimg, "y", 2)
    tshoe_y[[2]] <- mirror(tshoe_y[[2]], "y")
    
    tshoe_yfix <- lapply(tshoe_y, function(xx) {
      xxmod <- isoblur(xx, sigma)
      yy <- apply(xxmod, 2, mean) 
      zz <- which(abs(yy - max(yy)) < tol) %>% max()   
      
      # Don't crop if removing more than 25% of the image
      if (zz > .5*length(yy)) {
        zz <-  -1
      }
      imsub(xx, y > (zz + 1)) 
    }) 
    
    tshoe_yfix[[2]] <- mirror(tshoe_yfix[[2]], "y")
    modimg <- imappend(tshoe_yfix, axis = "y")
    rm(tshoe_y, tshoe_yfix)
  }
  
  return(modimg)
}

shoe_crop <- map_il(shoe_lbg, crop_border)
plot(shoe_crop, "row")

