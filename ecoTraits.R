# List of packages
packages <- c("lidR", "lidRmetrics", "Rdimtools", "Lmoments", "fitdistrplus", 
              "alphashape3d", "tidyr", "purrr", "concaveman", "sf", "viewshed3d", 
              "pracma", "ggplot2", "fsbrain", "rayshader", "rgl", "terra")

# Function to check if a package is installed and install it if not
check_and_install <- function(pkg){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg, dependencies = TRUE)
  }
}

# Apply the function to the list of packages
sapply(packages, check_and_install)

# Load the required packages
library(lidR)
# library(devtools) # Development tools for R

# Install additional packages from GitHub
# devtools::install_github("ptompalski/lidRmetrics")
library(lidRmetrics) # Additional metrics for lidR

# devtools::install_github("kisungyou/Rdimtools")
library(Rdimtools) # Box dimensions

library(Lmoments) # Required for lidRmetrics package
library(fitdistrplus) # Weibull pdf
library(alphashape3d) # Alphashapes

library(tidyr) # Alphashapes
library(purrr) # Functional programming tools
library(concaveman) # LeaningDistance, canopy area & radius
library(sf) # Simple features - spatial data handling

# Install hyper.fit from archive https://cran.r-project.org/src/contrib/Archive/hyper.fit/
# remotes::install_github('Blecigne/viewshed3d', dependencies = TRUE)

# Load the viewshed script
source("viewshed.R")

library(pracma) # Practical numerical math functions
library(ggplot2) # Data visualization
library(fsbrain) # Managing and visualizing brain surface data

library(terra) # Spatial data analysis

# Define Functions

# Function to calculate canopy ratio
canopy_ratio <- function(Z) {
  # Calculate the 25th and 98th percentiles of Z
  rh <- quantile(Z, probs = c(0.25, 0.98))
  
  # Calculate the canopy ratio
  cr <- (rh[[2]] - rh[[1]]) / rh[[2]] # (98-25)/98
  
  # Return the canopy ratio
  return(list(CP = cr))
}

# Function to calculate the Effective Number of Layers (ENL)
ENL <- function(Z, minh = 1) { # default set to 2D and a minimum height 1m
  # Initialize a dataframe to store the pi values
  pi <- data.frame(p = NA)
  
  # Calculate the number of elements in Z that are greater than minh
  p0 <- length(which(Z > minh))

  # Loop to calculate pi
  for (i in (minh + 1):ceiling(max(Z))) {
    p <- (length(which(Z <= i & Z > i - 1)) / p0)
    frame <- data.frame(p = p)
    pi <- rbind(pi, c(frame))
  }

  # Remove the first row of pi (which is NA)
  pi <- pi[-1,]

  # Calculate the 2D ENL
  twoD <- 1 / sum(pi^2) #2D ENL formula
  
  # Return the 2D ENL
  return(list(ENL = twoD))
}

# Function to calculate the Box Dimension
bd_fun <- function(X, Y, Z) {
  # Combine X, Y, and Z into a matrix
  M <- cbind(X, Y, Z)
  
  # Calculate the box dimension
  Db <- list(box_dimension = Rdimtools::est.boxcount(M)$estdim)

  # Return the box dimension
  return(Db)
}

# Function to calculate the Alpha shape
ashape <- function(x, y, z, id) {
  # Combine x, y, and z into a dataframe
  p <- data.frame(x, y, z)

  # Center the points around 0,0,0
  p[, 1] <- p[, 1] - mean(p[, 1])
  p[, 2] <- p[, 2] - mean(p[, 2])
  p[, 3] <- p[, 3] - mean(p[, 3])

  # Create a unique matrix of the points
  xyz <- data.matrix(unique(p))

  # Calculate the alpha shape
  ashape3d_obj <- alphashape3d::ashape3d(xyz, alpha = ALPHA, pert = TRUE)

  # Get the volumes from your different parameters
  ashape <- volume_ashape3d(ashape3d_obj, indexAlpha = "all")

  # Convert ashape to a dataframe and spread it
  ashape <- as.data.frame(ashape) %>% spread(ashape, ashape)

  # Name the columns
  colnames(ashape) <- c('A1', 'A2', 'A3')

  # Calculate differences in the volumes
  ashape$a2a1 <- ashape$A2 - ashape$A1
  ashape$a3a1 <- ashape$A3 - ashape$A1
  ashape$a3a2 <- ashape$A3 - ashape$A2

  # Return the ashape dataframe
  return(ashape)
}

# Function to calculate the 2D shape based on (Owen et al., 2021)
canopyShape <- function(x, y, z) {
  # Check if GEOM is either "convex" or "concave"
  if (GEOM == "convex" | GEOM == "concave") {
    # Combine x, y, and z into a matrix
    p <- cbind(x, y, z)

    # Center the points around 0,0,0
    p[, 1] <- p[, 1] - mean(p[, 1])
    p[, 2] <- p[, 2] - mean(p[, 2])

    # Order the points by height
    m <- as.matrix(p)
    m <- m[order(m[, 3], decreasing = FALSE)]
    b <- m[1:10,]

    # Create 2D simple features of the tree and base
    tree_sf <- st_as_sf(data.frame(m), coords = 1:2)
    base_sf <- st_as_sf(data.frame(b), coords = 1:2)

    # Create concave polygons of the tree canopy and base
    canopy_sf <- concaveman::concaveman(tree_sf, concavity = 2)
    base_sf <- concaveman::concaveman(base_sf, concavity = 1)

    # Calculate the center of the canopy and base
    canopy_cop <- st_centroid(canopy_sf)
    base_cop <- st_centroid(base_sf)

    # Calculate the distance between the center of the canopy and the center of the base
    delta_d <- st_distance(canopy_cop, base_cop)

    # Calculate the area of the canopy
    area <- st_area(canopy_sf)

    # Calculate the verticals ("edges") of the canopy polygon
    geopoints <- st_cast(canopy_sf, "POINT")
    d <- st_distance(st_centroid(canopy_sf), geopoints)

    # Calculate the radius of the canopy as mean and max distance from the center of the canopy
    r_mean <- mean(d)
    r_max <- max(d)

    # Create a list with the calculated values
    canopy_shape <- list(CP = as.numeric(delta_d),
                         Ac = area,
                         rcmean = r_mean,
                         rcmax = r_max)

    # Return the canopy shape
    return(canopy_shape)
  } else {
    stop("GEOM must be either 'convex' or 'concave'")
  }
}

# Set Variables
DZ <- 1 # Set vertical layer height in meters
GEOM <- "concave" # 'crown_metrics' returns results as concave polygon
ALPHA <- c(0.5, 1, 2) # Set your alpha parameters here!

# Load Data
LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR") 
las <- readLAS(LASfile, filter = "-drop_z_below 1") # Read the file, dropping points below 1m

# Visualize Trees
plot(las, bg = "white", size = 4, color = "treeID") 

# Check LAS data
las_check(las)

# Clean Data
las <- filter_poi(las, treeID >= 1, Classification == 1) # Remove ground and understory points
las <- filter_duplicates(las) # Remove duplicate points

# Subset point cloud into big and small trees based on the number of points
POINT_NR_THRES <- 100 # Set threshold for number of points per tree
tree_list <- table(las$treeID) # Get a table of tree IDs

# Identify big and small trees based on the threshold
bigtree_list <- tree_list > POINT_NR_THRES
smalltree_list <- tree_list < POINT_NR_THRES

# Get the names of the big trees
big_trees_names <- as.numeric(names(tree_list[bigtree_list]))

# Subset all the trees with more than POINT_NR_THRES points
big_trees <- filter_poi(las, treeID %in% big_trees_names) 

# Only big trees
las <- big_trees

# Clean up the workspace
rm(big_trees)

# Calculate standard metrics for the crown and convert to a dataframe
trait_std <- crown_metrics(las, ~stdmetrics_z(Z, dz = DZ), geom = GEOM)
mainframe <- as.data.frame(trait_std)

# Reorder columns and select the first 5 rows
mainframe <- mainframe %>% dplyr::select(treeID, geometry, tidyselect::everything())
mainframe <- mainframe[1:5]

# Plot the maximum height of the crowns
terra::plot(trait_std["zmax"], pal = hcl.colors, at = c(10, 20, 30, 40, 50), main = "", axes = FALSE)

# Calculate L-moments for the crown and join with the mainframe
trait_Lmoments <- crown_metrics(las, ~metrics_Lmoments(Z), geom = GEOM)
mainframe <- dplyr::inner_join(mainframe, trait_Lmoments, by = c("treeID", "geometry"))

# Calculate metrics based on the leaf area density and join with the mainframe
trait_lad <- crown_metrics(las, ~metrics_lad(Z, dz = DZ, k = 0.25, z0 = 1), geom = GEOM)
mainframe <- dplyr::inner_join(mainframe, trait_lad, by = c("treeID", "geometry"))

# Plot the leaf area index of the crowns
terra::plot(trait_lad["lai"], pal = hcl.colors, main = "", axes = FALSE)

# Calculate voxel-based metrics for the crown and join with the mainframe
trait_voxel <- crown_metrics(las, ~metrics_voxels(X, Y, Z, vox_size = 1), geom = GEOM)
mainframe <- dplyr::inner_join(mainframe, trait_voxel, by = c("treeID", "geometry"))

# Remove unwanted columns from the mainframe
mainframe <- subset(mainframe, select = -c(vn, vFRall, vFRcanopy, vzrumple, vzsd, vzcv))

# Calculate the canopy ratio and join with the mainframe
trait_cp <- crown_metrics(las, ~canopy_ratio(Z), geom = GEOM)
mainframe <- dplyr::inner_join(mainframe, trait_cp, by = c("treeID", "geometry"))

# Calculate the effective number of layers and join with the mainframe
trait_enl <- crown_metrics(las, ~ENL(Z), geom = GEOM)
mainframe <- dplyr::inner_join(mainframe, trait_enl, by = c("treeID", "geometry"))

# Calculate box dimensions and join with the mainframe
trait_bd <- crown_metrics(las, ~bd_fun(X, Y, Z), geom = GEOM)
mainframe <- dplyr::inner_join(mainframe, trait_bd, by = c("treeID", "geometry"))

# Calculate alpha shapes and join with the mainframe
trait_avol <- crown_metrics(las, ~ashape(X, Y, Z, treeID), geom = GEOM)
mainframe <- dplyr::inner_join(mainframe, trait_avol, by = c("treeID", "geometry"))

# Calculate canopy size and canopy displacment and join with the mainframe
trait_2dShape <- crown_metrics(las, ~canopyShape(X, Y, Z), geom = GEOM)
mainframe <- dplyr::inner_join(mainframe, trait_2dShape, by = c("treeID", "geometry"))

# Plot the canopy displacment
terra::plot(trait_2dShape["CP"], pal = hcl.colors, nbreaks = 128, breaks = "pretty", at = c(round(min(trait_2dShape$CP), 0), 2, 4, 6, 8, round(max(trait_2dShape$leanDistance), 0)), main = NULL)

# Calculate tree visibility and join with the mainframe
trait_vis <- tree_visability(las)
mainframe <- dplyr::inner_join(mainframe, trait_vis, by = c("treeID", "geometry"))

# Create a map of tree visibility
trait_vis_map <- trait_std["zmax"]
trait_vis_map$zmax <- mainframe$visability

# Plot the tree visibility map
terra::plot(
  trait_vis_map["zmax"], 
  pal = hcl.colors, 
  main = "", 
  axes = FALSE, 
  mar = c(0, 0, 0, 0)
)
# Get the number of columns in the mainframe
nr_traits <- ncol(mainframe)

# Replace infinite values with 0 in the mainframe
mainframe[3:nr_traits] <- do.call(data.frame, lapply(mainframe[3:nr_traits], function(value) replace(value, is.infinite(value), 0)))

# Replace NA values with 0 in the mainframe
mainframe <- do.call(data.frame, lapply(mainframe, function(value) replace(value, is.na(value), 0)))

# Loop through each tree
for (i in seq_along(trait_2dShape$treeID)) {
  # Define the neighborhood radius as four times the canopy max radius
  r_cp <- trait_2dShape$canopymaxRadius[i] * 4
  if (r_cp < 5) {
    r_cp <- 5
  }

  # Define the neighborhood area and find the trees within this area
  neighbourhood <- st_buffer(st_centroid(trait_2dShape[i, 1]), r_cp)
  competition <- st_intersection(neighbourhood, trait_2dShape[, 1])

  # Select trees with height >90% of the focal tree
  competition <- competition[which(mainframe[which(mainframe$treeID %in% competition$treeID.1), "zmax"] > mainframe[which(mainframe$treeID == trait_2dShape[i, 1]$treeID), 3] * 0.9),]

  # Skip if there is only one tree in the competition
  if (length(competition$treeID) == 1) { next }

  # Remove the target tree from the competition
  competition <- competition[-c(which(competition$treeID == competition$treeID.1)),]

  # Loop through traits/columns and calculate the max and mean for each trait in the neighborhood
  for (t in 3:nr_traits) {
    trait_name_max <- paste0("max_", colnames(mainframe[t]))
    trait_value_max <- max(mainframe[which(mainframe$treeID %in% competition$treeID.1), t])
    trait_name_mean <- paste0("mean_", colnames(mainframe[t]))
    trait_value_mean <- mean(mainframe[which(mainframe$treeID %in% competition$treeID.1), t])

    # Add the max and mean values to the mainframe
    if (trait_name_max %in% colnames(mainframe)) {
      mainframe[i, which(colnames(mainframe) == trait_name_max)] <- trait_value_max
    } else {
      mainframe$max[i] <- trait_value_max
      colnames(mainframe)[ncol(mainframe)] <- trait_name_max
    }
    if (trait_name_mean %in% colnames(mainframe)) {
      mainframe[i, which(colnames(mainframe) == trait_name_mean)] <- trait_value_mean
    } else {
      mainframe$mean[i] <- trait_value_mean
      colnames(mainframe)[ncol(mainframe)] <- trait_name_mean
    }
  }

  # Calculate 2D and 3D competition, Tree Density Index, minimum distance to neighbors, and pressure
  total_area <- st_area(neighbourhood)
  comp_area <- st_area(st_union(competition))
  C2d <- comp_area / total_area
  total_volume <- max(las$Z) * pi * r_cp^2
  neighbors_vol <- sum(mainframe$A3[mainframe$treeID %in% competition$treeID.1])
  comp3D <- neighbors_vol / total_volume
  mainframe$C2d[i] <- C2d
  mainframe$C3d[i] <- C3d
  mainframe$TDI[i] <- (length(competition$treeID.1)) / (pi * r_cp^2)
  distance <- as.numeric(sf::st_distance(st_centroid(trait_2dShape[i, 1]), st_centroid(competition)))
  mainframe$min_dis_neighbors[i] <- min(distance)
  canopylevelindex <- (0.5 / (1 + exp(-(mainframe$zmax[i] - mainframe$zmax[mainframe$treeID %in% competition$treeID.1] + 2)))) * (exp(-((mainframe$zmax[i] - mainframe$zmax[mainframe$treeID %in% competition$treeID.1])^2 / (2 * 10^2))) + 1)
  mainframe$pressureN[i] <- sum(((mainframe$CP[mainframe$treeID %in% competition$treeID.1] + distance) * canopylevelindex) / distance^2)
}

# Clip the Tree Density Index data
trait_sd <- trait_std["zmax"]
trait_sd$zmax <- clip.data(mainframe$TDI, lower = 0.01, upper = 0.9)

# Plot the clipped data
terra::plot(trait_sd["zmax"], pal = hcl.colors, main = "", axes = FALSE, mar = c(0, 0, 0, 0))


# Remove 'treeID' and 'geometry' columns from the mainframe
mainframe <- subset(mainframe, select = -c(treeID, geometry))

# Replace infinite values with 0 in the mainframe
mainframe <- do.call(data.frame, lapply(mainframe, function(value) replace(value, is.infinite(value), 0)))

# Save the mainframe to a CSV file
write.csv(mainframe, "output.csv")
