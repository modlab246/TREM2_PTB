# save_channel.ipynb 

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
