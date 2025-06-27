# save_channel.ipynb 

We want to construct composite images for each GeoMx spot 

1. A read folder for metadata - "../data/"
2. A read folder for channels - "../channel/"
3. A write folder for crops - "../crops/"
4. A write folder for composites - "../composites/"

## read in protein channels 
PTB 21.1: DAPI, CD68, CD3, CD8

PTB 22.1/2/3: DAPI, CD68, CD3, CD20 

## read in coordinates, write window around each spot, build composites and save 
1. scale raw to denom and clip
2. take exponential and scaling factor 
3. combine shared channels and construct composite 

| Sample     | window | red_denom | blue_denom | green_denom | yellow_denom | blue_adjust | yellow_adjust | red_adjust | green_adjust | gamma |
|------------|-----------|-----------|------------|-------------|--------------|-------------|---------------|------------|--------------|-------|
| **ptb 22.1** | 512 | 2**16     | 2**12      | 2**16       | 2**12        | 0.7         | 0.7           | 1.9        | 1.6          | 1.2   |
| **ptb 22.3** | 512 |  2**16     | 2**12      | 2**16       | 2**12        | 0.7         | 0.7           | 1.9        | 1.6          | 1.2   |
| **ptb 22.2** | 512 | 2**16     | 2**14      | 2**15       | 2**13        | 0.9         | 0.8           | 1.3        | 1.6          | 1.1   |
| **ptb 21.1** | img.shape[0]/2 | 2**13     | 2**14      | 2**12       | 2**13        | 1.2         | 1.1           | 1.8        | 1.8          | 1.4   |
## additive RGB
https://en.wikipedia.org/wiki/RGB_color_model

### PTB 22.1/2/3
red = red + yellow

green = green + yellow

blue = blue

* DAPI is mapping to blue as the blue file
* CD68 is mapping to red as the red file
* CD3 is mapping to green as the green file
* CD20 is mapping to yellow as the yellow file 
### PTB 21.1 
red = red 

green = blue (green file) + yellow (blue file)

blue = green (yellow file) + yellow (blue file)

* DAPI is mapping to green as the yellow file
* CD68 is mapping to cyan as the blue file
* CD3 is mapping to blue as the green file
* CD8 is mapping to red as the red file

# save_bb.ipynb 
We want to find the bounding boxes for PTB 22.1/3, whereas PTB 21.1 and PTB 22.2 have moderately high quality ROIs already saved. PTB 22.1/3 have bounding box annotations with separate coordinates, hence we have to find a scale factor mapping the GeoMx spots to their new coordinates.

$$
annotated_{XY} = metadata_{XY} 	\cdot sf + constant
$$

## Find white pixels 
PTB 22.1/2/3: [255, 255, 255, 255]

PTB 21.1: [254/5, 254/5, 254/5, 255] 
* A little less sensitive for this batch

## Drawing the box 
We find the white pixels within the GeoMx spot window, and take the point nearest to the center as the start. Note that PTB 21.1 and PTB 22.2 have numeric annotations already provided, but we want to write our own. 

Add white points to a set such that they are all nearest neighbors and stop when the set stops increasing. 

Construct the box as a scatter plot of the subset of white points, which are then smoothed over. 

## Saving masks and images
1. Save the mask for bounding box
2. Save composite with mask
3. Save the original annotation
4. Save the subset of the composite within the mask

# load_211/221/222/223_mat.ipynb 

Using the same intensities as save_channel.ipynb, we save full mosaic plots of the tissue sections with spots highlighted with and without labels. 

Because of the size of the channel files, individual matrices are deleted and memory is freed up. This was run on a computer with 64 GB of RAM. 

## Fill in bound box
Using contours to preserve outer shape

## Place text and validate 
Text should be adjacent to each ROI, and we run through several plots of the downsized matrix to make sure this happens. 

## Saving downsized figure
We take 1/16th of the full img ~100 MB

