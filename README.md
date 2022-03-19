# Nanomechanical_by_machine_learning

## Description
This respository is the machine learning implementation to predict mechanical properties (yiled stress, yield strain) of carbon nano-coils (CNC). Find more information in our paper "TBD".

## Data
The data for this project are taken from molecular dynamics (MD) simulations that has been done in ref. [1]. The data can be found in the "Data/model" folder for both elastic and plastic deformation.

The elastic data contains these columns (from left to right):
1. Mean radius of the CNC.
2. Coil diameter of the CNC
3. Pitch angle of the CNC
4. Pitch length	of the CNC
5. Number of turns
6. Total length of the CNC
7. Yiled stress
8. Yiled strain

The plastic data contains these columns (from left to right):
1. Mean radius of the CNC.
2. Coil diameter of the CNC
3. Pitch angle of the CNC
4. Pitch length	of the CNC
5. Number of turns
6. Total length of the CNC
7. Yiled stress
8. Fracture strain
9. Toughness of the CNC

"Data/RAW" contains the unprocessed data which is not cleaned.

## Code
The 

If you are interested in coordinates of a CNC without the lenghty optimization process you can define six inputs in the helix.xyz MATLAB m file in the main folder.
For a complete definition of these inputs refer to ref. [2]










[1] Shahini, E., Rangriz, F., Karimi Taheri, A., & Abdi-Jalebi, M. (2021). Optimizing Structural and Mechanical Properties of Coiled Carbon Nanotubes with NSGA-II and Reactive Molecular Dynamics Simulation. The Journal of Physical Chemistry C, 125(11), 6237-6248.
[2] Chuang, Chern, Yuan-Chia Fan, and Bih-Yaw Jin. "On the structural rules of helically coiled carbon nanotubes." Journal of Molecular Structure 1008 (2012): 1-7.
