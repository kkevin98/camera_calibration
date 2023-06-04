# camera_calibration

Project for the "Computer vision and pattern recognition 2020/21" course at University of Trieste.
Goals of the project, given a set of photos, are:
1. Calibrate the camera that took the photos using Zhang procedure
2. Refine the results obtained in point 1 taking into account radial distorsion
3. Superimpose an object in the photos

## Structure

- `/images`: folder that store the photos used in the project
- `Kevin_Marzio_project_function.m`: Matlab script that solves the 3 goals
- `CVPR_Kevin_Marzio.pdf`: report containing problem statement, description of the used approach, implementation
choices, results and references
- Other `.m` files are Matlab functions used inside `Kevin_Marzio_project_function.m`

## Usage

Execute `Kevin_Marzio_project_function.m` one chunk at a time to get the results  
