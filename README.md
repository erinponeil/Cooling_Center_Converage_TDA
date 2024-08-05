This repository provides the data and Python notebooks needed to recreate the analysis of the coverage of cooling centers in Boston, MA as presented in the paper "Evaluating Holes in Cooling Center Coverage Using Persistent Homology of a Filtered Witness Complex."

It containes two folders, each with their own subfiles/folders:
  1) **Data**
      - <ins>Shapefile and Centroid (Landmark) Data:</ins> Contains the Shapefiles for plotting maps (from which we created a dataset of landmarks using the .centroid attribute of the Shapely library)
      - <ins>Cooling Center Location (Witness) Data:</ins> Contains the notebook needed to scan OpenStreetMap to collect the latitude and longitude coordinates of cooling centers within a given domain
      - <ins> HVI Score Data and Shapefile: </ins> Contains the cleaned demographic data, notebook, and Shapefiles to make the HVI maps seen in the paper
  3) **Python Notebooks**
      - The notebook that creates the witness complex and performs persistent homology

