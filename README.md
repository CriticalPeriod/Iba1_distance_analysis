# VGLUT2 PV PNN 

* **Developed for:** Camille
* **Team:** Prochiantz
* **Date:** Feb 25
* **Software:** Fiji



### Images description

3D images taken with a x60 objective

4 channels:
  1. *Alexa Fluor 647:* PV cells
  2. *Alexa Fluor 546:* PNN cells
  3. *EGFP:* foci protein 
  4. *Hoechst 33342:* nuclei

### Plugin description

* Detect nuclei, PV and PNN cells with Cellpose
* Compute their colocalization
* Keep PV or PNN cells with a nucleus only
* Detect DAPI  PV/PNN nucleus with Stardist
* Detect GFP foci and intensity in PV/PNN membrane with Stardist
### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **StarDist** *pmls2.zip* (homemade) model
* **Cellpose** conda environment + *cyto2*, *cyto_PV1* (homemade) and *livecell_PNN1* (homemade) models

### Version history

Version 1 released on November 17, 2022.

