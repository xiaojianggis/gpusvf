## Mapping heat exposure using GPU parallel computing
Urban microclimate modeling helps to quantitatively understand the fine level urban heat exposure distribution, which would provide an important reference for urban heat management. Fine level continuous sky view factor (SVF) maps are usually needed for modeling how the solar radiation fluxes reach the urban surface and impact the urban microclimate dynamics. However, the estimation of continuous SVF maps is very time-consuming, which limits the urban microclimate modeling to small geographical areas. In this study, we proposed to use graphics processing unit (GPU) parallel computing to accelerate the computing of SVF. The high-resolution digital surface models were used as the input to generate the continuous SVF maps of Philadelphia, Pennsylvania, USA. This repository includes the scripts to compute the SVFs using the GPU-accelerated algorithms.

by Xiaojiang Li, Temple University

### 1. Data preparation
The datasets needed in this repo is digital surface model and land cover dataset. In Philadelphia, both of these dataset can be accessed from PASDA website [Link](https://www.pasda.psu.edu). Just type in the keyword of `land cover` and `digital surface model`, you will be able to find the dataset and download them. 

### 2. Prepare the computing environment
In order to use GPU to accelearte the computing the SVF, you need to setup the computing environment. First have the Anaconda installed on your command. We are going to use Anaconda to install the required modules. Once the Anaconda install, then you can setup the environment using the following command, 
`conda create --name climategpu numpy shapely matplotlib rasterio fiona pandas ipython pyproj gdal jupyter`

### 3. Compute the SVF based on building height model
You can compute the SVF based on the building height model. There are two most commonly used algorithm for compute SVF values, shadow casting-based algorithm, ray-tracing algorithm. 

1. Shadow casting algorith [link](cpu/shadow-casting-svf.py)
2. Ray-tracing algorithm [link](cpu/ray-tracing-svf.py)

### 4. GPU parralel computing for SVF
Using GPU parallel computing, you can make the SVF computing much much faster!!!!!!!!!!!!! Here are two GPU-accelerated algorithms, GPU-based shadow casting algorithm, and GPU-based ray-tracing algorithm. 

1. GPU-based shadow casting [link](gpu/pycuda-svf-shadowcasting.ipynb)
2. GPU-based ray-tracing algorithm [link](gpu/pycuda-svf-raytracing.ipynb)


