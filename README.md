# Automated script for magnitude-and-phase-based T2 mapping method

Authors: Difei Wang, DZNE Bonn

### Processing pipeline

Data processing is performed by the script *T2T1AM.py*. To call the script use

```
 *T2T1AM.py <magnitude> <phase> <mask> <b1map> -tr TR -fa FA -phi delta_phi -outputdir path*
```
where 
- *magnitude* is the filename of the magnitude image
- *phase* is the filename of the phase image
- *mask* is the filename of the mask
- *b1map* is the filename of the B1 map
- *TR* is the repetition time in ms
- *FA* is the flip angle in degrees
- *delta_phi* is the list of RF phase increments
- *path* is the path to the output folder.

The pipeline expects the following data:

- B<sub>1</sub><sup>+</sup> map in percentage unit
- mask of the brain
- background-removed 4D magnitude images
- background-removed 4D phase images

The outputs are T<sub>2</sub>, T<sub>1</sub>, amplitude maps.



### Demo Dataset

An example dataset of the phantom can be found in [here](https://osf.io/fknyh/). 


### Required package

The EPG package can be found in [here](https://github.com/mrphysics-bonn/EPGpp).

The Python packages are: 
- numpy
- nibabel
- scipy.optimize
- os
- glob
- pymp
- argparse
- sys
- time / datetime (optional)
