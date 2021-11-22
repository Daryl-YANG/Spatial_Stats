###########################################################################################
#
#        bootstrap spatial variogram for raster images
#
#    --- Last updated:  2021.11.20 By Daryl Yang <dediyang@bnl.gov>
###########################################################################################

#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#*****************************************************************************************#

#****************************** load required libraries **********************************#
### install and load required R packages
list.of.packages <- c("ggplot2", 'tidypaleo', 'raster', 'rgdal', 'gstat', 'sp', 'matrixStats')  
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=c("Depends", "Imports",
                                                                       "LinkingTo"))
version_requirements <- c("3.3.2")
if (!packageVersion("ggplot2") >= version_requirements[1]) {
  remotes::install_version(package="ggplot2", version=paste0(">= ", version_requirements), 
                           dependencies=c("Depends", "Imports", "LinkingTo"), upgrade="ask",
                           quiet=TRUE)
}
# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#*****************************************************************************************#

#**************************** define semivariance function *******************************#
bootstrap_variogram <- function(raster, iter, width)
{
  # convert raster file to dataframe for semivariance
  semivar.DF <- as(raster, 'SpatialPointsDataFrame')
  # remove na values in the dataframe
  semivar.DF <- na.omit(semivar.DF)
  
  # perform an iteration of semivariance based on the number defined in iter
  range.DF <- c()
  dist.orig.DF <- c()
  dist.pred.DF <- c()
  semivar.orig <- c()
  semivar.pred <- c()
  for (i in 1:iter)
  {
    print(paste0('current iteration ', i))
    # sample data from the raster dataframe using defined perc
    sample.index <- sample(1:nrow(semivar.DF), round(nrow(semivar.DF)*perc))
    sample.data <- semivar.DF[sample.index, ]
    # calculate variogram
    var <- gstat::variogram(spatialX ~ 1, data = sample.data, cutoff = 8000, 
                            width = width)
    # store semivar values
    dist.orig.DF <- cbind(dist.orig.DF, var$dist)
    semivar.orig <- cbind(semivar.orig, var$gamma)
    print(semivar.orig)
    # fit varigram
    # if want the variogram model fit psill, nugget, and range, please delete
    # them 
    # potential models: c('Sph', 'Exp', 'Mat', 'Ste', 'Gau')
    var.fit <- gstat::fit.variogram(var, vgm(psill = 0.7, nugget = 0.02, range = 1000,
                                             model = 'Sph',
                                             kappa = 0.25))
    # store range parameter
    range.DF <- rbind(range.DF, var.fit$range)
    
    # predict semivariance
    semivar.p <- variogramLine(var.fit, maxdist = max(var$dist))
    # store predicted semivar values
    dist.pred.DF <- cbind(dist.pred.DF, semivar.p$dist)
    semivar.pred <- cbind(semivar.pred, semivar.p$gamma)
  }
  # calculate mean and variance of original semivariance and predicted semivariance
  semivar.orig.DF <- data.frame(rowMeans(dist.orig.DF), rowMeans(semivar.orig), rowSds(semivar.orig))
  names(semivar.orig.DF) <- c('dist', 'mean', 'sd')
  semivar.pred.DF <- data.frame(rowMeans(dist.pred.DF), rowMeans(semivar.pred), rowSds(semivar.pred))
  names(semivar.pred.DF) <- c('dist', 'mean', 'sd')
  # put all data into a list and return
  output <- list('range' = range.DF[, 2], 'semivar_orig' = semivar.orig.DF, 
                 'semivar_pred' = semivar.pred.DF)
  return(output)
}
#*****************************************************************************************#

#************************************ user parameters ************************************#
# define output directory
outDIR <- "/Volumes/GoogleDrive/My Drive/manuscripts/Yang_et_al_2021_Tall_Shru_ELF/figures/figure_3"

# define the percentage of data to sample from the raster file for semivariance calculation
# this parameter is set based on the size of your raster file. the resulting data points 
# used for semivariance better to be less than 100000.
perc <- 0.01

# define the number of iterations for calculate mean semivariance model
iter <- 10

# define the step width used to calucate semivariance
width = 50 
#*****************************************************************************************#

#*************************************** load data ***************************************#
### load data
# load raster data
raster.DIR <- "/Volumes/data2/dyang/projects/shrub_elf/council/analysis/pft_fcover_map/fcover_mosaic/fcover_mosaic_5m"
raster.RST <- brick(raster.DIR)

### preprocess the imagery before calculating semivariance (!!!not required)
### clip the raster based on a shp file
shp.DIR <- "/Volumes/data2/dyang/projects/shrub_elf/council/analysis/boundary/aviris_ng_small.shp"
shp.VCT <- readOGR(shp.DIR)
# project the vector file to target projection system to match with raster files
shp.VCT <- spTransform(shp.VCT, crs(raster.RST))
# clip the raster file
raster.CLIPPED <- crop(raster.RST, shp.VCT)
raster.MASKED <- mask(raster.CLIPPED, shp.VCT)

# assign the processed data back to raster.RST
raster.RST <- raster.MASKED
#*****************************************************************************************#

#************************************ semivariance ***************************************#
# calculate semivariance for alnus
alnus.RST <- raster.RST$Resize..Band.2.fcover_mosaic.tif.
names(alnus.RST) <- 'spatialX'
semivar.alnus <- bootstrap_variogram(alnus.RST, iter, width)
# calculate semivariance for salix
salix.RST <- raster.RST$Resize..Band.3.fcover_mosaic.tif.
names(salix.RST) <- 'spatialX'
semivar.salix <- bootstrap_variogram(salix.RST, iter, width)
#*****************************************************************************************#

#************************************** make plot ****************************************#
ggplot(data = NULL) +
  geom_point(data = semivar.alnus$semivar_orig, aes(x = dist, y = mean), size = 1,
            color = 'navyblue') +
  geom_line(data = semivar.alnus$semivar_pred, aes(x = dist, y = mean), size = 1,
            color = 'navyblue') +
  geom_ribbon(data = semivar.alnus$semivar_pred, aes(x = dist, ymin = mean-sd, ymax = mean+sd),
              fill = 'black', alpha = 0.3) +
  geom_vline(xintercept = mean(semivar.alnus$range), color = 'navyblue', linetype = 'dashed',
             size = 1) +
  #geom_ribbonh(aes(xmin = min(semivar.alnus$range), xmax = max(semivar.alnus$range), 
  #                 y = semivar.alnus$semivar_pred$mean), alpha = 0.3) +
  geom_point(data = semivar.salix$semivar_orig, aes(x = dist, y = mean), size = 1,
             color = 'green4') +
  geom_line(data = semivar.salix$semivar_pred, aes(x = dist, y = mean), size = 1,
            color = 'green4') +
  geom_ribbon(data = semivar.salix$semivar_pred, aes(x = dist, ymin = mean-sd, ymax = mean+sd),
              fill = 'black', alpha = 0.3) +
  geom_vline(xintercept = mean(semivar.salix$range), color = 'green4', linetype = 'dashed',
             size = 1) +
  #geom_ribbonh(aes(xmin = min(semivar.salix$range), xmax = max(semivar.salix$range), 
  #                 y = semivar.alnus$semivar_pred$mean), alpha = 0.3) +
  theme(legend.position = 'none') + ylim(0.02, 0.13) + xlim(0, 6000) +
  ylab('FCover Semivariance') + xlab('Distance (m)') +
  theme(axis.text = element_text(size=11, color = 'black'),
        axis.title=element_text(size=12)) +
  scale_y_continuous(breaks = seq(0.02, 0.13, 0.02)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pngName = paste0(outDIR, "/",'semivariance.pdf')
ggsave(pngName, plot = last_plot(), width = 10, height = 8, units = 'cm')












