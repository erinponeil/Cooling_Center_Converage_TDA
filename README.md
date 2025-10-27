This repository provides the data and Python notebooks needed to recreate the analysis of the coverage of cooling centers in Boston, MA as presented in the paper "Evaluating Holes in Cooling Center Coverage Using Persistent Homology of a Filtered Witness Complex" which can be found here: https://arxiv.org/abs/2410.09067.

It containes two folders, each with their own subfiles/folders:
  1) **Data**
      - <ins>Shapefile and Centroid (Landmark) Data:</ins> Contains the Shapefiles for plotting maps (from which we created a dataset of landmarks using the .centroid attribute of the Shapely library)
      - <ins>Cooling Center Location (Witness) Data:</ins> Contains the notebook needed to scan OpenStreetMap to collect the latitude and longitude coordinates of cooling centers within a given domain
      - <ins>HVI Score Data and Shapefile: </ins> Contains the cleaned demographic data and Shapefiles to make the HVI maps seen in the paper
      - <ins>Processed Data</ins>
  2) **Figs**
      - <ins>HVI_Figs</ins>


