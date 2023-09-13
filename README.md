# Density and Triangulation-based Clustering (DTC)

This repository is the official implementation of the Density and Triangulation-based Clustering algorithm from the paper "Triangulation based Spatial Clustering for Adjacent Data with Various Density". This algorithm aims at solving clustering difficulties faced by complex data with irregular shapes, various densities, and adjacent connections between multiple clusters.

## DTC data experiment notebook:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1-u9kJkxvodytEU4nta3stBDyj3BR-Cvi?usp=sharing)

## Usage
```python
alg = triangulation_dbscan(data[['x','y']], kde = True)
clusters_df = alg.tri_dbscan()
```

### Parameters:

**data**: dataframe
- The input dataframe with column names 'x' and 'y' representing the coordinates of the spatial data

**minPts**: int, default=5 
- The number of samples in a neighborhood for a point to be considered as a core point. This does not include the point itself.

**kde**: bool, defult=False
- Create initial clusters based on the kernel density estimation of the data before using triangulation-based DBSCAN. It aims at solving the various density problem in complex data.
- If the density of the data is uniformly distributed, this option can be set to False to improve efficiency.

**local_std**: float, default = 2.5
- The standard deviation used as threshold to locally remove outlying triangles.

**figsize**: (float, float), default = (8,8)
- Width, height in inches of the final output.

**progress**: bool, defult = False
- Allow the displayment of the detailed algorithm progress.
