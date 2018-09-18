# Low resolution scans can provide a sufficiently accurate, cost- and time-effective alternative to high resolution scans for 3D shape analyses - data and code
Authors: Ariel E. Marcy, Carmelo Fruciano, and Vera Weisbecker

This repository contains all of the data and analyses needed to recreate the results from [the manuscript](https://peerj.com/articles/5032/)

To cite the paper and code:
> Marcy AE, Fruciano C, Phillips MJ, Mardon K, Weisbecker V. (2018) Low resolution scans can provide a sufficiently accurate, cost- and time-effective alternative to high resolution scans for 3D shape analyses. PeerJ 6:e5032 https://doi.org/10.7717/peerj.5032

## Data
Landmarking data:
[3D meshes for all skulls](https://www.morphosource.org/Detail/ProjectDetail/Show/project_id/458) in the study available via MorphoSource
[Raw_Coordinates.csv](3d-vs-ct-scanning/data/Raw_coordinates.csv) - the shape coordinates from landmarking both 3D and uCT scanned skulls in Viewbox 

Running generalized Procrustes analysis (GPA) with bilateral symmetry requires several files to establish the relationship between landmarks:
* [Smatrix.csv](3d-vs-ct-scanning/data/Smatrix.csv) - slider matrix required by geomorph to slide curve semi-landmarks properly
* [Bilateral_Landmarks.csv](3d-vs-ct-scanning/data/Bilateral_Landmarks.csv) - table identifies bilaterally symmetrical landmarks

Re-running GPA with different sets of landmarks required different versions of the bilateral landmark table:
* [Bilateral_LM_SM_only.csv](3d-vs-ct-scanning/data/Bilateral_LM_SM_only.csv) - table for dataset without patch points (i.e. only fixed landmarks and curve semi-landmarks)
* [Bilateral_LM_only.csv](3d-vs-ct-scanning/data/Bilateral_LM_only.csv) - table for dataset with only fixed landmarks (i.e. no curve or patch semilandmarks)

Sex data:
* [Pse_sex.csv](3d-vs-ct-scanning/data/Pse_sex.csv) - table identifying which skulls are male or female 
    
## Analyses
All analyses are contained in the Rscript, [Rscript_Marcyetal_2018.R](3d-vs-ct-scanning/Rscript_Marcyetal_2018.R). Comments within explain each stage and the corresponding Figure and Table numbers in the manuscript. 

Before starting remember to either set your working directory to the 3d-vs-ct-scanning folder on your computer, or open an RStudio project from that folder.
