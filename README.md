This repository provides the data and Python notebooks needed to recreate the analysis of the coverage of cooling centers in Boston, MA as presented in the paper "Evaluating Holes in Cooling Center Coverage Using Persistent Homology."

It containes two folders, each with their own subfiles/folders:
  1) **Data**
      - <ins>Shapefile and Centroid (Landmark) Data:</ins> Contains the Shapefiles for plotting maps, the notebook to find the centroids of census tracts, and the final dataset of landmarks used
      - <ins>Cooling Center Location (Witness) Data:</ins> Contains the notebook needed to scan OpenStreetMap to collect the latitude and longitude coordinates of cooling centers within a given domain and the final dataset of witnesses used
  3) **Python Notebooks**
      - The notebook that creates the witness complex of the two datastes and performs persistent homology
