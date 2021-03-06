Instructions for usage:

The code and data files here accompany the manuscript “An integrated occupancy and home-range model to predict abundance of a wide-ranging, territorial vertebrate” by Morgan W. Tingley, Robert L. Wilkerson, Christine A. Howell, and Rodney B. Siegel.

Any questions should be directed to Morgan Tingley at <morgan.tingley@uconn.edu>

The code predicts total abundance and relative density of Black-backed Woodpeckers within post-fire landscapes of the Sierra Nevada and southern Cascade mountains of California. The code predicts abundance using three component models (an occupancy model, a snag model, and a telemetry model), and integrates the models to predict abundance onto a new landscape. In this case, an example fire is illustrated using GIS data from the Reading fire in California. 

The two codes (“Area_abundance_predict” and “Fire_density_predict”) do two similar but slightly different operations. Both apply the same model but have different outputs. 

The first ("Area_abundance_predict”) outputs the total abundance of woodpecker pairs expected within a defined region (defined by a single shapefile polygon). In the example, it calculates the total number of woodpecker pairs expected within the entire Reading fire. The output is a posterior distribution of this number. Thus, the most likely estimate is the mode, but the mean can also be informative, as can be the 95% interval.

The second (“Fire_density_predict”) calculates the expected density per pixel (in pairs), as well as the standard deviation (i.e. uncertainty) of that density. As most pixels are quite small relative to a home range size, the density is often far less than 1. The code outputs both modal density (i.e., the expected density) and the standard deviation of density as raster layers (GeoTiff) that can be easily imported into other GIS programs. As we have found that areas within fires with higher predicted densities are more likely to actually contain detectable woodpeckers, we recommend that the modal density raster can be used for intra-fire management planning.
