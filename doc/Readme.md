# DeuteMetCon
This toolbox is intended to map deuterium metabolites with different offresonance with either full signal model equation for bSSFP or just phase evolution in case of spoiled acquisition.

# features 
- zero padding
- Signal models
    - Phase evolution
    - bSSFP signal model
    - 
- flexible solvers 
    - pinv
    - IDEAL (only for spoiled)
- Coil Combination
    - ESPIRIT
    - adaptive combine
    - wsvd?
- retrospective undersampling

# analysis
- Residual
- condition number
- simulation

# dependencies
- adaptive Coil Combination
- SPM (image registration)


# Todo
- [ ]  Treat CSI data same as bSSFP
- [ ]  can we include voxel chemical shift into the fit
- [ ]  