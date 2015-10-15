# A primer to close-range image analysis with R
Charles A. Martin  
Updated on 2015-10-15

For any questions, please email me at charles.martin1@uqtr.ca

##1. Install the required packages

This primer relies on a couple of CRAN packages, as well as some packages only available on GitHub or BitBucket.
It is thus necessary to use the devtools package to install the last ones.


```r
  install.packages("raster")
  install.packages("rgdal")
  
  install.packages("devtools")
  library("devtools")  
  devtools::install_github("cmartin/LAI")
  devtools::install_github("cmartin/EXIFr")
  devtools::install_bitbucket("persican/imagemetrics", ref="36d34cc")
```
N.B. If you have the necessary tools to install packages that use Rcpp code, you can remove the "ref" 
  argument and use the much faster calculateHistoCPP function.
  
##2. Pictures taken parallel to the ground
###2.1. Mean Information Gain (MIG)

MIG is a measure of texture, based information theory, proposed by Proulx and Parrott (2008). MIG represents the amount of spatial randomness in the patterns of an image channel. This heterogeneity measure is usually calculated with a single diagonal neighbor for each cell.

MIG values are normalized between 0 and 1, where MIG value close to 0 are associated with structural patterns and values close to 1 are more random.

MIG values tend to average around 0.4 for natural scenes.

In this example, we will measure MIG on Hue, Saturation and Value channels, as well on a grey-scale one and the number of bins used for the histograms will be k=15

###2.2. Anisotropy (orientation)

Anisotropy is a measure used to determine if heterogeneity patterns, as measured by MIG, are more heterogeneous along a particular direction. This is usually done by calculating the ratio of MIG~horizontal~/MIG~vertical~, that is, MIG calculated with the neighbor to the right of the cell divided by the one calculated from the neighbor below the cell.

This measure should give you some hindsights about the presence of major horizontal or vertical structures.

###2.3. Color indices

Color indices are metrics commonly used in agriculture to assess crop biomass. The main idea is to extract pixels of a specific color (in our case green) from the image, and calculate what percentage of the whole image these pixels represent. The function used here (getBinaryVegetationMask) implements the ExG-ExR method (Meyer and Neto 2008).

Additional methods also exist to extract additional colors, for example yellow to measure flowering phenology.

###2.4. Code example

Note : this code calculates MIG on multiple channels only as an example of how to do such calculations.
The resulting measures will be highly correlated. In order to avoid spurious results, one
should, based on strong hypotheses, decide _a priori_ which color channels to use.


```r
  library(raster)
  library(imagemetrics)
  
  path <- "TestImages/0/IMG_5401.jpg" 
  
  #############################################################################
  #### STEP 1 : Load or calculate all image layers we need to work with.
  #############################################################################
  
  # Load 3 rasters from the target image, 
  # one for red, one for green and one for the blue channel
  R <- raster(path, band = 1)
  G <- raster(path, band = 2)
  B <- raster(path, band = 3)
  
  # Combine the RGB channels to create a grayscale image
  RGB <- brick(R,G,B) 
  r_grey <- mean(RGB)
  
  # Convert RGB bands to HSV channels
  HSV <- rgb2hsv(getValues(R),getValues(G),getValues(B))  	
  r_H <- r_S <- r_V <- raster(ncols = ncol(R), nrows = nrow(R))
  extent(r_H) <- extent(r_S) <- extent(r_V) <- extent(R)
  values(r_H) <- HSV[1,]
  values(r_S) <- HSV[2,]
  values(r_V) <- HSV[3,]  
  
  # Create a binary image seperating vegetation from the background
  # This functions needs the whole 3-color stack
  r_green <- getBinaryVegetationMask(RGB) 
  
  #############################################################################
  #### STEP 2 : Calculate neighbor matrices and histograms on each layer. 
  ####          We will then feed these into the actual analysis functions
  #############################################################################
  
  # On the four channels, find either right (1), diagonal (2) or below (3) 
  # neighbors for histogram calculations
  v_grey_1 <- getImagePixels(r_grey, side = 1)
  v_H_1 <- getImagePixels(r_H, side = 1)
  v_S_1 <- getImagePixels(r_S, side = 1)
  v_V_1 <- getImagePixels(r_V, side = 1)
  
  v_grey_2 <- getImagePixels(r_grey, side = 2)
  v_H_2 <- getImagePixels(r_H, side = 2)
  v_S_2 <- getImagePixels(r_S, side = 2)
  v_V_2 <- getImagePixels(r_V, side = 2)
  
  v_grey_3 <- getImagePixels(r_grey, side = 3)
  v_H_3 <- getImagePixels(r_H, side = 3)
  v_S_3 <- getImagePixels(r_S, side = 3)
  v_V_3 <- getImagePixels(r_V, side = 3)
  
  # Calculate histograms from neighbor vectors
  nbins <- 15
  
   prob_grey_1 <- calculateHisto(
     reference_vector = v_grey_1$reference_vector,
     neighbour_vector = v_grey_1$neighbour_vector,
     nbins = nbins
  )
  prob_H_1 <- calculateHisto(
    reference_vector = v_H_1$reference_vector,
    neighbour_vector = v_H_1$neighbour_vector,
    nbins = nbins
  )
  prob_S_1 <- calculateHisto(
    reference_vector = v_S_1$reference_vector,
    neighbour_vector = v_S_1$neighbour_vector,
    nbins = nbins
  )
  prob_V_1 <- calculateHisto(
    reference_vector = v_V_1$reference_vector,
    neighbour_vector = v_V_1$neighbour_vector,
    nbins = nbins
  )
  
  prob_grey_2 <- calculateHisto(
    reference_vector = v_grey_2$reference_vector,
    neighbour_vector = v_grey_2$neighbour_vector,
    nbins = nbins
  )
  prob_H_2 <- calculateHisto(
    reference_vector = v_H_2$reference_vector,
    neighbour_vector = v_H_2$neighbour_vector,
    nbins = nbins
  )
  prob_S_2 <- calculateHisto(
    reference_vector = v_S_2$reference_vector,
    neighbour_vector = v_S_2$neighbour_vector,
    nbins = nbins
  )
  prob_V_2 <- calculateHisto(
    reference_vector = v_V_2$reference_vector,
    neighbour_vector = v_V_2$neighbour_vector,
    nbins = nbins
  )
  
  prob_grey_3 <- calculateHisto(
    reference_vector = v_grey_3$reference_vector,
    neighbour_vector = v_grey_3$neighbour_vector,
    nbins = nbins
  )
  prob_H_3 <- calculateHisto(
    reference_vector = v_H_3$reference_vector,
    neighbour_vector = v_H_3$neighbour_vector,
    nbins = nbins
  )
  prob_S_3 <- calculateHisto(
    reference_vector = v_S_3$reference_vector,
    neighbour_vector = v_S_3$neighbour_vector, 
    nbins = nbins
  )
  prob_V_3 <- calculateHisto(
    reference_vector = v_V_3$reference_vector,
    neighbour_vector = v_V_3$neighbour_vector,
    nbins = nbins
  )
  
  #############################################################################
  #### STEP 3 : Calculate the metrics.
  #############################################################################
  
  # Calculate the green index
  sum(getValues(r_green)) / length(getValues(r_green))
  
  # Calculate mean information gain
  meanInformationGain(prob_grey_2)
  meanInformationGain(prob_H_2)
  meanInformationGain(prob_S_2)
  meanInformationGain(prob_V_2)
  
  # Calculate anisotropy
  meanInformationGain(prob_grey_1) / meanInformationGain(prob_grey_3)
  meanInformationGain(prob_H_1) / meanInformationGain(prob_H_3)
  meanInformationGain(prob_S_1) / meanInformationGain(prob_S_3)
  meanInformationGain(prob_V_1) / meanInformationGain(prob_V_3)
```

```
[1] 0.1817417
[1] 0.6257871
[1] 0.241573
[1] 0.7049987
[1] 0.6305297
[1] 1.097495
[1] 1.208036
[1] 1.065547
[1] 1.105934
```

##3. Pictures taken at an angle

###3.1 Leaf Area Index (LAI)

LAI is defined as the ratio of total leaf surface divided by the area the canopy is covering on the ground. The higher the LAI number, the more stratified the vegetation is.

Here, we use an indirect LAI measure which, to be technically correct, actually measures plant area index (i.e. it cannot differentiate between trunks and leafs).

It it based on the fact that (theoretically), when looking up in a forest at an angle of 57.5°, the leaf orientation (relative to the observer) is random. Then, through a couple of mathematical hoops, one can translate the gap fraction (% of open sky) at that angle to an LAI value (MacFarlane 2011).

In order to do these calculations, the image is first cropped around the 57.5° angle (+/- 5°) and converted to a binary mask, where each pixel is assigned either to vegetation or sky. This is done by implementing Rosin's unimodal thresholding (2001), adding modifications suggested by MacFarlane (2011).

###3.2 Code example

N.B. This code assumes a couple of things about the lens and picture properties which need to be properly checked before running your own analysis (i.e. Camera FOV = 73.7° and focal angle = 45°)

```r
  library(LAI)
  path <- "TestImages/45/IMG_5428.jpg" 
  LAI_from_gf_at_57(path)
```

```
[1] 1.795992
```

##4. Canopy images

###4.1 Gap fraction (canopy opening)

Using Rosin's (2001) binary thresholding algorithm, one can also simply extract the ratio of sky pixels divided by the total number of pixels, thus obtaining the canopy opening.

###4.2 Code example

```r
  library(LAI)
  path <- "TestImages/90/IMG_7822.JPG" 
  
  # Unimodal thresholding is designed to work on the blue channel only
  B <- raster(path, band = 3)
  binary_image <- unimodal_threshold(B)
  gap_fraction(binary_image)
```

```
[1] 0.2976806
```

##5. Using camera settings as additional information

Using information from the EXIF header, one can add aditional insights into the image analysis process

###5.1 How to analyze variation in light conditions

If all but one of the camera settings were fixed during the experiments (e.g. fixing ISO, aperture and focal length but allowing exposure time to vary) while using automatic exposure modes, one can directly interpret variations in that "free" setting as changes in light conditions.

###5.2 Code example

```r
  library(EXIFr)
  image1 <- "TestImages/90/IMG_0009.JPG" 
  image2 <- "TestImages/90/IMG_7822.JPG" 
  
  rational_to_numeric(read_exif_tags(image1)[["ApertureValue"]])
  rational_to_numeric(read_exif_tags(image2)[["ApertureValue"]])
  
  read_exif_tags(image1)[["ISOSpeedRatings"]]
  read_exif_tags(image2)[["ISOSpeedRatings"]]
  
  rational_to_numeric(read_exif_tags(image1)[["FocalLength"]])
  rational_to_numeric(read_exif_tags(image2)[["FocalLength"]])
  
  read_exif_tags(image1)[["ExposureTime"]]
  read_exif_tags(image2)[["ExposureTime"]]
```

```
[1] 5.375
[1] 5.375
[1] 800
[1] 800
[1] 15
[1] 15
[1] "1/2500"
[1] "1/1600"
```


##6. How to put all these calculations in a loop you can stop and resume at will?

Image analysis can be a lengthy process that runs for many hours (and even days) for large image sets.
Here, we introduce a way to code analysis loops that can be stopped and resumed at any time.
After your first run, you will need to erase "results.csv" if you wish to reanalyse the first images.


```r
  library(EXIFr)
  library(LAI)
  
  # This is the path to the folder containing images to analyze
  images_folder <- "TestImages/0"
  
  # This is where the results will be written
  results_file <- "results.csv"
  
  # This code checks to see if the loop needs to resume or start from scratch
  if (file.exists(results_file)) {
    existing <- read.csv(results_file)
  	start <- max(existing$i) + 1
  } else {
  	start <- 1
  }
  
  files <- dir(images_folder)
  nb_photos <- length(files)
  
  if (start <= nb_photos) {
  
    for (i in start:nb_photos) {
    	file_to_analyze <- files[[i]]
    	path <- paste(images_folder,file_to_analyze,sep = "/")
      
    	B <- raster(path, band = 3)
      
    	binary_image <- unimodal_threshold(B)

      # Write to the CSV file after every image is analyzed
    	write.table(
    		data.frame(
    			i = i,
    			ID = file_to_analyze,
    			GF = gap_fraction(binary_image),
    			ExposureTime = read_exif_tags(path)[["ExposureTime"]]
    		),
    		file = results_file,
    		append = i != 1,
    		col.names = i == 1,
    		row.names = FALSE,
    		sep = ","
    	)	
    }
  }
  
  # Peek at the results
  cat(readLines(results_file), sep = "\n")
```

```
"i","ID","GF","ExposureTime"
1,"IMG_4955.jpg",0.00660960960960961,"1/15"
2,"IMG_5401.jpg",0.0642522522522523,"1/25"
3,"IMG_5691.jpg",0.777821321321321,"1/160"
4,"IMG_6198.jpg",0.0455780780780781,"1/4"
```

##7. A word of caution

Most of the above techniques have only been tested in a small sample of forests, thus they are still, at best, good case studies. 

From our experience, although they are not direct 1:1 equivalent to field sampling, they should still provide interesting ecological insights.

##8. Key references

Macfarlane, Craig. "Classification method of mixed pixels does not affect canopy metrics from digital images of forest overstorey." Agricultural and forest meteorology 151.7 (2011): 833-840.

Meyer, George E., and João Camargo Neto. "Verification of color vegetation indices for automated crop imaging applications." Computers and Electronics in Agriculture 63.2 (2008): 282-293.

Proulx, Raphaël, and Lael Parrott. "Measures of structural complexity in digital images for monitoring the ecological signature of an old-growth forest ecosystem." Ecological Indicators 8.3 (2008): 270-284.

Rosin, Paul L. "Unimodal thresholding." Pattern recognition 34.11 (2001): 2083-2096.
