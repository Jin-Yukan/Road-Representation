# Road-Representation

This repository provides tools to simplify road network geometries and represent their topology in different forms.

## Features

### Spatial datasets (road units)
Three types of spatial datasets are generated to represent different road units:

1. **Segments**  
   - Roads are split at intersections.  

2. **Strokes**  
   - Road units are constructed based on visual continuity.  

3. **Mixed map**  
   - Road units are first aggregated by street names.  
   - For roads without names, strokes are generated to ensure completeness.

### Network datasets
Four types of network datasets are generated, including one primal graph and three dual graphs corresponding to different road units:

1. **Primal graph**  
   - Intersections are represented as nodes, and road segments as edges.  

2. **Dual graphs**  
   - Road units (Segments, Strokes, Mixed map) are represented as nodes.  
   - Edges indicate intersection relationships between road units.  

## Usage

1. Install and configure **geopandas** and **networkx** within the **ArcGIS Pro 3.0** Python environment. Detailed versions are as follows:
   - GeoPandas 0.9.0
   - Pandas 1.3.5
   - Shapely 2.0.1
   - Fiona 1.8.21
   - PyProj 2.6.1
   - NetworkX 2.7.1
   - NumPy 1.20.1
2. Prepare your road network dataset (e.g., from [OpenStreetMap](https://www.openstreetmap.org/)).  
3. Run the scripts to process the road network. An example is provided in the `main` function.  

## Acknowledgements

- Stroke generation code is adapted from [momepy](https://github.com/pysal/momepy). The source code is licensed under the BSD 3-Clause License (see LICENSES file).  
- Built on top of **arcpy**, **geopandas** and **networkx**.  
