## ------------------------------------------------------------------------------------
## Description
## ------------------------------------------------------------------------------------

## Author: Renee Bichler, 2025
## Find me on GitHub: reneebichler

## This code uses the hourly 3D EQUATES data and extracts the surface concentrations
## of NO2 and converts it into a georeferenced GeoTIFF file.
## YOU CAN DO IT!! ╰( ^o^)╮╰( ^o^)╮

cat("##################################################################################\n")
cat("##\n")
cat("##                        Hourly CMAQ 3D EQUATES Data for CONUS\n")
cat("##                                    to GeoTIFF!\n")
cat("##                                 ╰( ^o^)╮╰( ^o^)╮\n")
cat("##\n")
cat("##################################################################################\n")

## ------------------------------------------------------------------------------------
## Libraries
## ------------------------------------------------------------------------------------

library(ncdf4)
library(terra)
library(sf)
library(dplyr)
library(scales)
library(rnaturalearth)
library(ggplot2)

## ------------------------------------------------------------------------------------
## Variables
## ------------------------------------------------------------------------------------

output_dir <- "/DATA/EQUATES"
polygon_path <- "/DATA/GEODATA/s_18mr25/s_18mr25.shp"

exclude <- c(
    "Alaska", "American Samoa", "Hawaii", "Puerto Rico", "Marshall Islands",
    "Fed States of Micronesia", "Rhode Island", "Virgin Islands", "Guam", "Palau",
    "Northern Mariana Islands"
)

## Set the EPSG code for the map
epsg_code <- 5071

## t_val and l_val stand for the time step and the layer
time_val <- 9
layer_val <- 1

## ------------------------------------------------------------------------------------
## Natural Earth data for map
## ------------------------------------------------------------------------------------

## States and Provinces
ne_10m_admin_1_states_provinces <- ne_states(returnclass = "sf")
ne_10m_admin_1_states_provinces <- st_transform(ne_10m_admin_1_states_provinces, crs = 4326)

## Populated places
ne_10m_populated_places <- ne_download(scale = 10, type = "populated_places", category = "cultural", returnclass = "sf")
ne_10m_populated_places <- st_transform(ne_10m_populated_places, crs = 4326)
ne_10m_populated_places <- subset(ne_10m_populated_places, subset = SOV_A3 == "USA")

## Airports
ne_10m_airports <- ne_download(scale = 10, type = "airports", category = "cultural", returnclass = "sf")
ne_10m_airports <- st_transform(ne_10m_airports, crs = 4326)

## Ports
ne_10m_ports <- ne_download(scale = 10, type = "ports", category = "cultural", returnclass = "sf")
ne_10m_ports <- st_transform(ne_10m_ports, crs = 4326)

## Lakes
ne_10m_lakes <- ne_download(scale = 10, type = "lakes", category = "physical", returnclass = "sf")
ne_10m_lakes <- st_transform(ne_10m_lakes, crs = 4326)

## Rivers
ne_10m_rivers_lake_centerlines <- ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical", returnclass = "sf")
ne_10m_rivers_lake_centerlines <- st_transform(ne_10m_rivers_lake_centerlines, crs = 4326)

## ------------------------------------------------------------------------------------
## Functions
## ------------------------------------------------------------------------------------

create_topmap <- function(polygon) {

    ne_10m_lakes$feature <- "Waterbodies"
    ne_10m_admin_1_states_provinces$feature <- "States"
    polygon$feature <- "U.S. States"

    ## Administrative boundaries for states, counties, and AOI
    topmap <- list(
        geom_sf(data = ne_10m_lakes, aes(color = feature, fill = feature, linetype = feature), linewidth = 0.2),
        geom_sf(data = ne_10m_admin_1_states_provinces, aes(color = feature, fill = feature, linetype = feature), linewidth = .2),
        geom_sf(data = polygon, aes(color = feature, fill = feature, linetype = feature), linewidth = .2),
        scale_color_manual(values = c("U.S. States" = "black", "States" = "grey75", "Waterbodies" = "steelblue3")),
        scale_fill_manual(values = c("U.S. States" = "transparent", "States" = "transparent", "Waterbodies" = "steelblue1")),
        scale_linetype_manual(values = c("U.S. States" = "solid", "States" = "solid", "Waterbodies" = "solid")),
        guides(
            color = guide_legend(title = "Area Features", ncol = 1, order = 1),
            fill = guide_legend(title = "Area Features", ncol = 1, order = 1),
            linetype = guide_legend(title = "Area Features", ncol = 1, order = 1)
        ),
        ggnewscale::new_scale_color(),
        ggnewscale::new_scale_fill()
    )
    return(topmap)
}

## ------------------------------------------------------------------------------------
## Main
## ------------------------------------------------------------------------------------

## List nc_files in the directory
cat("Listing NetCDF files in the directory...\n")
nc_files <- list.files(paste0(output_dir), pattern = "\\.nc$", full.names = TRUE)

geo_dir <- paste0(output_dir, "/geo")
map_dir <- paste0(output_dir, "/maps")

## Create "geo" folder if it doesn't exist
if (!dir.exists(geo_dir)) {
    cat("Creating output directory for GeoTIFF files...\n")
    dir.create(geo_dir, recursive = TRUE)
}

if (!dir.exists(map_dir)) {
    cat("Creating output directory for maps...\n")
    dir.create(map_dir, recursive = TRUE)
}

for (nc_file in nc_files) {
    
    cat("Processing file:", nc_file, "\n")

    ## Retrieve the basename and remove the .nc extension
    basename <- sub("\\.nc$", "", basename(nc_file))

    ## Open the NetCDF file
    cat("Opening NetCDF file using 'ncdf4' package...\n")
    nc <- nc_open(nc_file)

    ## Extract variables
    ## Start indices (COL=1, ROW=1, LAY=1, TSTEP=9)
    ## -1 means read all in that dimension; only 1 layer and timestep
    cat("Extracting variables from NetCDF file...\n")
    no2 <- ncvar_get(nc, varid = "NO2", start = c(1, 1, layer_val, time_val), count = c(-1, -1, 1, 1))  
    xcell <- ncatt_get(nc, 0, "XCELL")$value    # x resolution
    ycell <- ncatt_get(nc, 0, "YCELL")$value    # y resolution
    xorig <- ncatt_get(nc, 0, "XORIG")$value    # x origin
    yorig <- ncatt_get(nc, 0, "YORIG")$value    # y origin
    ncols <- ncatt_get(nc, 0, "NCOLS")$value    # xmax
    nrows <- ncatt_get(nc, 0, "NROWS")$value    # ymax
    palp <- ncatt_get(nc, 0, "P_ALP")$value     # xxx
    pbet <- ncatt_get(nc, 0, "P_BET")$value     # xxx
    ycent <- ncatt_get(nc, 0, "YCENT")$value    # y center
    xcent <- ncatt_get(nc, 0, "XCENT")$value    # x center
    tstep <- ncvar_get(nc, "TSTEP")

    ## Layer and time step information
    n_lay <- length(ncvar_get(nc, "LAY"))
    n_tstep <- length(tstep)

    ## Retrieve time information
    tstep_units <- ncatt_get(nc, "TSTEP")$units

    ## Close the NetCDF file
    cat("Closing NetCDF file...\n")
    nc_close(nc)

    ## Get time
    cat("Retrieving time information...\n")
    time_hour <- tstep[time_val]

    ## Retrieve the date and time stamp
    cat("Retrieving date and time stamp...\n")
    date_time_stamp <- format(as.POSIXct(tstep_units, format = "hours since %Y-%m-%d %H:%M:%S", tz = "UTC") + time_hour * 3600, "%Y-%m-%d %H:%M:%S")

    ## Convert to raster
    cat("Converting NO2 data to raster...\n")
    r <- rast(t(no2))

    ## Rename the layer
    r <- setNames(r, "layer")

    ## Change the extent of the raster
    cat("Setting raster extent...\n")
    ext(r) <- ext(0, ncols * xcell, 0, nrows * ycell)

    ## Retrieve the projection information
    cat("Setting raster CRS...\n")
    crs_lcc <- paste0("+proj=lcc +lat_1=", palp," +lat_2=", pbet, " +lat_0=", ycent, " +lon_0=", xcent, " +x_0=", xorig*(-1), " +y_0=", yorig*(-1), " +a=6370000 +b=6370000 +units=m +no_defs")
    crs(r) <- crs_lcc

    ## Reverse the order of the raster
    cat("Reversing the order of the raster...\n")
    r <- flip(r, direction = "vertical")

    ## Write to GeoTIFF
    cat("Writing raster to GeoTIFF...\n")
    writeRaster(r, paste0(geo_dir, "/", basename, "_geo_NO2_layer_", layer_val, ".tif"), overwrite = TRUE)

    ## Transform the raster layer to a different CRS
    cat("Transforming raster layer to EPSG:", epsg_code, "...\n")
    r_layer1 <- project(r, paste0("EPSG:", epsg_code))

    ## Load polygon
    polygon <- read_sf(paste0(polygon_path))

    ## Reproject polygon
    polygon <- st_transform(polygon, epsg_code)

    ## Remove all areas in exclude from polygon
    polygon <- polygon %>% filter(!polygon$NAME %in% exclude)

    ## Get the extent of the polygon
    extent <- as.numeric(st_bbox(polygon))

    ## Clip raster to polygon
    cat("Clipping raster to polygon...\n")
    r_layer1 <- crop(r_layer1, polygon)
    r_layer1 <- mask(r_layer1, polygon)

    ## Convert raster to data frame
    cat("Converting raster to data frame...\n")
    df <- as.data.frame(r_layer1, xy = TRUE, na.rm = TRUE)
    
    ## Convert data frame to sf object
    cat("Converting data frame to sf object...\n")
    sf <- st_as_sf(df, coords = c("x", "y"), crs = st_crs(epsg_code))

    ## Colorbar settings
    colors <- c(
        "magenta",          # ~0-
        "#08306B",          # 0
        "#08306B",          # 2
        "#5ea9ff",          # 6
        "#FFFF00",          # 8
        "#FFA500",          # 14
        "#FF0000",          # 28
        "#800080"           # 30
    )

    color_positions <- c(
        -0.000000000000001, # magenta
        0,                  # sharp transition to dark blue
        2,                  # dark blue
        6,                  # light blue
        8,                  # yellow starts
        14,                 # orange
        28,                 # red
        30                  # purple
    )

    limits <- c(0, 30)
    breaks <- c(0, 2, 4, 6, 8, 10, 15, 20, 25, 30)

    ## Normalize positions to 0–1 for gradientn
    rescaled_positions <- rescale(color_positions, to = c(0, 1), from = limits)

    ## Assign features to spatial objects
    ne_10m_lakes$feature <- "Waterbodies"
    ne_10m_admin_1_states_provinces$feature <- "States"
    polygon$feature <- "U.S. States"

    ## Extract the coordinates for geom_tile
    sf <- sf %>%
    st_centroid() %>%
    cbind(st_coordinates(.)) %>%
    as.data.frame()

    ## Create ggplot
    cat("Creating ggplot for visualization...\n")
    map <- ggplot() +
        geom_tile(data = sf, aes(x = X, y = Y, fill = layer)) +
        scale_fill_gradientn(
            name = expression("Surface NO"[2]*" (ppbV)"),
            colors = colors,
            values = rescaled_positions,
            limits = limits,
            breaks = breaks,
            oob = scales::squish,
            na.value = "transparent",
            guide = guide_colorbar(
                draw.ulim = TRUE,
                draw.llim = TRUE,
                barwidth = 12,
                barheight = .5,
                title.theme = element_text(size = 10),
                label.theme = element_text(size = 8),
                text = element_text(family = "serif")
            )
        ) +
        ggnewscale::new_scale_fill() +
        geom_sf(data = ne_10m_lakes, aes(color = feature, fill = feature, linetype = feature), linewidth = 0.2) +
        geom_sf(data = ne_10m_admin_1_states_provinces, aes(color = feature, fill = feature, linetype = feature), linewidth = .2) +
        geom_sf(data = polygon, aes(color = feature, fill = feature, linetype = feature), linewidth = .2) +
        scale_color_manual(values = c("U.S. States" = "black", "States" = "grey75", "Waterbodies" = "steelblue3")) +
        scale_fill_manual(values = c("U.S. States" = "transparent", "States" = "transparent", "Waterbodies" = "steelblue1")) +
        scale_linetype_manual(values = c("U.S. States" = "solid", "States" = "solid", "Waterbodies" = "solid")) +
        guides(
            color = guide_legend(title = "Area Features", ncol = 1, order = 1),
            fill = guide_legend(title = "Area Features", ncol = 1, order = 1),
            linetype = guide_legend(title = "Area Features", ncol = 1, order = 1)
        ) +
        ggnewscale::new_scale_color() +
        ggnewscale::new_scale_fill() +
        geom_sf(data = polygon, fill = "transparent", color = "black", na.rm = TRUE) +
        coord_sf(
            crs = st_crs(epsg_code),
            xlim = c(round(extent[1]-1, digits = 0), round(extent[3]+1, digits = 0)),
            ylim = c(round(extent[2]-1, digits = 0), round(extent[4]+1, digits = 0))
        ) +
        labs(
            title = paste("CMAQ EQUATES -", date_time_stamp, "- L:", layer_val),
            subtitle = paste(basename),
            x = NULL, y = NULL
        ) +
        theme_bw() +
        theme(
            legend.position = "bottom",
            legend.title = element_text(size = 10),
            legend.title.position = "top",
            legend.text = element_text(size = 8),
            legend.key = element_rect(color = NA, fill = NA),
            legend.key.size = unit(0.2, "cm"),
            legend.key.height = unit(0.2, "cm"),
            legend.key.width = unit(0.4, "cm"),
            legend.background = element_blank(),
            text = element_text(family = "serif"),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 8),
            panel.grid.major = element_line(linetype = "dashed", color = "grey60", linewidth = 0.1),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank()
        )

    ## Save ggplot
    cat("Saving ggplot to PNG...\n")
    ggsave(
        filename = paste0(map_dir, "/", basename, "_NO2_layer1.png"),
        plot = map, width = 7, height = 4.5, bg = "transparent"
    )
}

cat("##################################################################################\n")
cat("##\n")
cat("##                                      Success!\n")
cat("##                                  Job completed!\n")
cat("##                                     ( ┘^o^)┘\n")
cat("##\n")
cat("##################################################################################\n")
