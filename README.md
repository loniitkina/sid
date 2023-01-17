# sid
Sea ice floe edges become indistinguishable inside the winter pack ice of the central Arctic. 
They freeze together into clusters or groups of floes that move relative to each other only during sea ice deformation events driven by weather patterns.
Consequently, the boundaries and shapes of these clusters change on a synoptic scale. 
These boundaries are typically several kilometers long and nearly straight – commonly named linear kinematic features (LKFs).
Geometric characteristics of these ‘dynamic elements’ between the LKFs and their spatio-temporal development can be therefore tracked by studying sea ice deformation concentrated along LKFs.

In this code we use sea ice drift retrievals from displacements between SAR image pairs. The method is applied to a collection SAR imagery of Sentinel-1.
We derive a deformation-rate thresholding method that enables relatively high spatial resolution of  ~200m and bellow.
We use the sequences of deformation fields to: 
1) track sea ice cover than underwent sea ice deformation (damaged ice), and 
2) derive LKF area, intersection angle, and active floe area and shape. 
