# BC-RTI

This repository contains supporting material for the paper: "Block Condensing for Fast Nonlinear MPC with the Dual Newton Strategy". More precisely, the chain of masses benchmark problem is solved with block condensing and qpDUNES within ACADO. Users can also compare the results with the dense solver qpOASES or the structure exploiting solver FORCES. For the latter, a license for the software FORCES Pro is required.

## Instructions

First clone the *master* branch of ACADO and compile the MATLAB interface (following the instructions in http://www.acadotoolkit.org). To ensure that the code will run without problems despite any changes in the ACADO software, a zipped folder with the version used for the results of the paper is also provided.

To run the example, simply open the script NMPC_chain_mass.m in MATLAB and edit the simulation options in the first lines according to your interests.
