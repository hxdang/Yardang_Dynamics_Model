# Yardang_Dynamics_Model
Yardang_Dynamics_Model is a comprehensive, fully coupled numerical modele for the flow field and morphology evolution dynamics of Yardangs basaed on Boussinesq types of equations and stochastic process parameterization.

***

The source code includes the following files:
* Yardang_Dynamics_Model.f90
    - This is the main program for the Yardang Dynamics Model.
    - The runtime and time step of the calculation are defined here.
    - The result data is output in text format or other specified formats.
* grid.f90
    - Generate the staggered grid.
* topography.f90
    - Customize the initial topography.
* ini_condition.f90
    - Set the initial conditions for the calculation.
* flow_field.f90
    - Calculate the flow field across the topography.
* sand_transport.f90
    - Calculte the sand transport.
* parameter_flow.txt
    - Set the parameters related to the flow field calculation.
* parameter_sand.txt
    - Set the parameters related to the sand transport calculation.
 
*The gfortran compiler must be installed before using the model source code.*

***

@Copyright Key Laboratory of Earth System Numerical Modeling and Application, University of Chinese Academy of Sciences, Beijing, China.

All rights reserved.

Created by Haoxuan Dang, 2024.

When using the code, please cite the associated papers.


