## Mapping heat exposure using GPU parallel computing
Urban microclimate modeling helps to quantitatively understand the fine level urban heat exposure distribution, which would provide an important reference for urban heat management. Fine level continuous sky view factor (SVF) maps are usually needed for modeling how the solar radiation fluxes reach the urban surface and impact the urban microclimate dynamics. However, the estimation of continuous SVF maps is very time-consuming, which limits the urban microclimate modeling to small geographical areas. In this study, we proposed to use graphics processing unit (GPU) parallel computing to accelerate the computing of SVF. The high-resolution digital surface models were used as the input to generate the continuous SVF maps of Philadelphia, Pennsylvania, USA. Based on the computed SVF maps, this study further calculated and examined the spatial distribution of mean radiation temperature (Tmrt) as an indicator of outdoor heat exposure using the SOlar and LongWave Environmental Irradiance Geometry (SOLWEIG) model. 

by Xiaojiang Li, Temple University

## 1. Data preparation
The datasets needed in this repo is digital surface model and land cover dataset. In Philadelphia, both of these dataset can be accessed from PASDA website. Just type in the keyword of `land cover` and `digital surface model`, you will be able to find the dataset and download them. 

