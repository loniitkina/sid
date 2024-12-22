# sid - Sea Ice Deformation from imaging remote sensing data
Sea ice floe edges become indistinguishable inside the winter pack ice of the central Arctic. 
They freeze together into clusters or groups of floes - Coherent Dynamic Clisters (CDCs), that move relative to each other only during sea ice deformation events driven by weather patterns.
Consequently, the boundaries and shapes of these clusters change on a synoptic scale. 
These boundaries are typically several kilometers long and nearly straight – commonly named Linear Kinematic Features (LKFs).
Geometric characteristics of these CDCs between the LKFs and their spatio-temporal development can be therefore tracked by studying sea ice deformation concentrated along LKFs.

In this code we use sea ice drift retrievals from displacements between SAR image pairs. The method is applied to a collection SAR imagery of Sentinel-1.
We derive a deformation-rate thresholding method that enables relatively high spatial resolution of  ~800m and bellow.
We use the sequences of deformation fields to: 
1) track sea ice cover than underwent sea ice deformation (damage parcles), and 
2) derive LKF and CDC geometry.

