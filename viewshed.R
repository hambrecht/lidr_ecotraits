## viewshed3d Lecigne et al. 2020


my_d_visibility <- function(data, position, angular_res, elevation_range, azimuth_range, scene_radius, store_points, accept_null) {

  #- declare variables to pass CRAN check as suggested by data.table maintainers
  N = theta = X = Y = Z = phi = r = . = Visibility = NULL

  #- default parameters
  if (missing(store_points)) store_points <- F
  if (missing(accept_null)) accept_null <- F
  if (missing(position)) position <- c(0, 0, 0)
  if (missing(angular_res)) angular_res <- 1
  if (missing(elevation_range)) elevation_range = c(0, 180)
  if (missing(azimuth_range)) azimuth_range = c(0, 360)

  #- check for possible problems in arguments
  if (missing(data)) stop("no input data provided.")
  if (class(data)[1] != "LAS") {
    stop("data must be a LAS. You can use lidr::LAS(data) to convert a
         data.frame or data.table to LAS.")
  }
  if (is.numeric(angular_res) == F) stop("angular_res must be numeric.")
  if (is.numeric(position) == F) stop("position must be numeric.")
  if (length(position) != 3) stop("position must be a vector of length 3.")
  if (is.logical(store_points) == F) stop("store_coord must be logical")
  if (length(elevation_range) != 2) stop("elevation_range must be a vector of length 2.")
  if (length(azimuth_range) != 2) stop("azimuth_range must be a vector of length 2.")
  if (min(azimuth_range) < 0 | max(azimuth_range) > 360) {
    stop("acceptable azimuth_range are between 0 and 360")
  }
  if (min(elevation_range) < 0 | max(elevation_range) > 180) {
    stop("acceptable elevation_range are between 0 and 180")
  }

  # sort elevation range, the smaller first
  elevation_range = sort(elevation_range)

  #- BUILD PARAMETER FILE
  #- with theta = the elevation angle, and N = the number of cells
  if (180 %% angular_res == 0) { # build the entire theta range
    param <- data.table::data.table(theta = seq(-90, 90, angular_res), N = NA)
  }else {
    param <- data.table::data.table(theta = seq(-90, 90 + angular_res, angular_res), N = NA)
  }
  # compute the number of cells in each latitudinal band following the
  # method described by Malkin (2016)
  L <- pi / (360 / angular_res) # cell side length at the equator
  N_cell <- round(pi / L) # number of cells at the equator
  # the number of decreases with increasing the elevation angle
  param[, N := round(N_cell * (cos(pracma::deg2rad(theta))))]

  #####################################################################
  ########## subset for the ranges in elevation and azimuth ###########
  #####################################################################

  # N = 1 for the two extremes and transform theta range from -90, 90 to 0, 180
  param[theta == 90 | theta == -90, N := 1]
  param[, theta := theta + 90]

  ########################################################################
  # computing the exact number of sightlines when azimuth_range is defined
  sphere <- generate_sphere(angular_res / 5, r = scene_radius)

  # follow the same process than with the data to compute the number of sightlines
  sphere[, ':='(theta = round((acos(Z * 100 / (sqrt(X^2 + Y^2 + Z^2) * 100)) *
    (180 / pi)) / angular_res) * angular_res, #-
                phi = acos(Y * 100 / (sqrt(Y^2 + X^2) * 100)) * (180 / pi),
                r = sqrt(X^2 + Y^2 + Z^2))]
  sphere[X < 0, phi := (180 - phi) + 180] # transform phi so its range is now 0-360#

  sphere = sphere[theta >= elevation_range[1] & theta <= elevation_range[2]]

  if (azimuth_range[1] < azimuth_range[2]) {
    sphere <- sphere[phi >= azimuth_range[1] & phi <= azimuth_range[2]]
  }else {
    sphere <- sphere[phi >= azimuth_range[1] | phi <= azimuth_range[2]]
  }

  data.table::setkey(param, theta)
  data.table::setkey(sphere, theta)
  sphere[, N := param[sphere, N]]

  sphere[, phi := round(phi / (360 / N)) * (360 / N)]
  sphere[phi >= 360, phi := 0]

  n_dir <- nrow(unique(sphere[, .(theta, phi)])) # the number of sightlines

  #####################################################################

  #- FIND THE NEAREST POINT IN EACH DIRECTION
  data <- data@data[, .(X, Y, Z)] # transform data into a data.table

  #- translate the data so the center is 0,0,0
  data[, ':='(X = X - position[1],
              Y = Y - position[2],
              Z = Z - position[3])]

  #- compute spherical coordinates: elevation angle (theta), azimuth (phi) and
  #- radius (r)
  data[, ':='(theta = round((acos(Z * 100 / (sqrt(X^2 + Y^2 + Z^2) * 100)) *
    (180 / pi)) / angular_res) * angular_res, #-
              phi = acos(Y * 100 / (sqrt(Y^2 + X^2) * 100)) * (180 / pi),
              r = sqrt(X^2 + Y^2 + Z^2))]
  data[X < 0, phi := (180 - phi) + 180] # transform phi so its range is now 0-360

  ################################################################
  # select the portion of the data

  # for elevation
  data <- data[theta >= elevation_range[1] & theta <= elevation_range[2]]

  # for azimuth
  if (azimuth_range[1] < azimuth_range[2]) {
    data <- data[phi >= azimuth_range[1] & phi <= azimuth_range[2]]
  }else {
    data <- data[phi >= azimuth_range[1] | phi <= azimuth_range[2]]
  }

  ################################################################

  if (missing(scene_radius) == F) { #- apply the scene radius if defined
    #- ensure scene.radius is valid
    if (is.numeric(scene_radius) == F) stop("scene.radius must be numeric.")
    data <- data[r <= scene_radius]
  }
  if (accept_null == T) {
    if (nrow(data) == 0) {
      near <- 0
      return(near)
    }
  }
  if (nrow(data) == 0) stop(paste("scene_radius is too small, the scene contains 0 point. TreeID"),id)

  #- match param and data
  data.table::setkey(param, theta)
  data.table::setkey(data, theta)
  data[, N := param[data, N]]

  #- compute the phi for each cell center
  data[, phi := round(phi / (360 / N)) * (360 / N)]
  data[phi >= 360, phi := 0]

  #- find the nearest point in each direction and store xyz coordinates
  near <- data[, .(r = min(r)), keyby = .(theta, phi)]
  near <- unique(near)

  #- store coordinates for export (needed for the cumview function)
  if (store_points == T) {
    data[, Visibility := 1]
    data.table::setkeyv(data, c("theta", "phi", "r"))
    data.table::setkeyv(near, c("theta", "phi", "r"))
    data[near, Visibility := 2]
    # data[,c("theta","phi","N"):=NULL]
    data[, N := NULL]
    data[, ':='(X = X + position[1], Y = Y + position[2], Z = Z + position[3])]
  }

  #- COMPUTE VISIBILITY
  #- number of nearest points at each distance
  near[, r := round(r, 2)]
  near <- near[, .(N = .N), keyby = r]

  #- add distances from 0 to the minimum recorded
  start <- data.table::data.table(r = seq(0, min(near$r), 0.1), N = 0)
  end <- data.table::data.table(r = seq(max(near$r), scene_radius, 0.1), N = 0)
  near <- rbind(start, near, end)

  #- compute visibility
  near[, visibility := (1 - cumsum(N) / n_dir) * 100]
  near[, N := NULL]

  if (store_points == T) {

    las = pkgcond::suppress_messages(lidR::LAS(data[, 1:3]))
    las = lidR::add_lasattribute_manual(las = las, x = data$r, name = "r", desc = "distance to the scene center", type = "float")
    las = lidR::add_lasattribute_manual(las = las, x = data$Visibility, name = "Visibility", desc = "is the point visible", type = "int")

    return(
      list(visibility = near, points = las)
    )
  }else {
    return(near)
  }
}

namevis <- function(z) { # random function
    zmean <- mean(z)
    return(list(visability = zmean))
  }

tree_visability <- function(las, scene_radius) {
  require(viewshed3d)



  trait_vis <- crown_metrics(las, ~namevis(Z), geom = "concave") # create template layer
  trait_loc <- crown_metrics(las, ~namevis(Z), geom = "point") # create layer to retrieve coordinates from
  tree_loc <- st_coordinates(trait_loc) # retrieve coordinates of trees, x,y are of center of tree and z is max height
  if (missing(scene_radius) == F) { #- apply the scene radius if defined
    #- ensure scene.radius is valid
    if (is.numeric(scene_radius) == F) stop("scene.radius must be numeric.")
  }
  if (missing(scene_radius) == T) {
    scene_radius <- max((las$"Max X" - las$"Min X"), (las$"Max Y" - las$"Min Y")) # if scence radius is missing, set it to length or width of the whole point cloud, what ever is longer
  }


  for (i in seq_along(trait_vis$treeID)) {
    vis <- viewshed3d::d_visibility(data = las, position = tree_loc[i,], angular_res = 1, scene_radius = scene_radius, store_points = F, elevation_range = c(0, 100)) # direction viewshed limited to 0--90 degrees
    visN <- viewshed3d::d_visibility(data = las, position = tree_loc[i,], angular_res = 1, scene_radius = scene_radius, store_points = F, azimuth_range = c(270, 90), elevation_range = c(0, 100)) # direction viewshed limited to 0--90 degrees elevation and 270--90 degrees azimuth
    # visN <- my_d_visibility(data = las, position = tree_loc[i,], angular_res = 1, scene_radius = scene_radius, store_points = F, azimuth_range = c(270, 90), elevation_range = c(0, 100), accept_null = TRUE) # direction viewshed limited to 0--90 degrees elevation and 270--90 degrees azimuth
    trait_vis$visability[i] <- pracma::trapz(x = vis$r, y = vis$visibility) # calculates the area under the curve of visibility
    if (length(visN) == 2) {
      trait_vis$visabilityN[i] <- pracma::trapz(x = visN$r, y = visN$visibility) # calculates the area under the curve of visibility
    }else { trait_vis$visabilityN[i] <- 0 }
  }
  return(trait_vis)
}