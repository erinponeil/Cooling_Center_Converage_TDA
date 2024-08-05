This repository provides the data and Python notebooks needed to recreate the analysis of the coverage of cooling centers in Boston, MA as presented in the paper "Evaluating Holes in Cooling Center Coverage Using Persistent Homology of a Filtered Witness Complex."

It containes two folders, each with their own subfiles/folders:
  1) **Data**
      - <ins>Shapefile and Centroid (Landmark) Data:</ins> Contains the Shapefiles for plotting maps, the notebook to find the centroids of census tracts, and the final dataset of landmarks used
      - <ins>Cooling Center Location (Witness) Data:</ins> Contains the notebook needed to scan OpenStreetMap to collect the latitude and longitude coordinates of cooling centers within a given domain and the final dataset of witnesses used
      - <ins> HVI demographic Data: </ins> Contains the data, notebook and Shapefiles to make the HVI maps seen in the paper
  3) **Python Notebooks**
      - The notebook that creates the witness complex of the two datastes and performs persistent homology

Note: Users must ensure they have verion **3.1** or later for the *networkx* package
