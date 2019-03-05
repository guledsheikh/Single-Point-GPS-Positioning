# Introduction
Accurate determination of position has become very important in recent years. To estimate the GNSS receiver’s position at any time, one determines the four unknown parameters, that is, the X, Y, Z coordinates of the receiver and its clock error. 

This project focuses on the implementation of static single point positioning (SPP). In static SPP, the receiver is stationary throughout the GPS data collection period. 

A software module for conducting SPP using Least Squares approximation is given here using GPS pseudorange measurements. The software utility includes the following main functions which needs to be added to ones project for succesful execution:
  1. SPP.cpp---performs epoch-by-epoch least square estimation of the position
  2. Satelliteposition.cpp---computes the satellite positions epoch-by-epoch
  3. Reader.cpp---reads the observation and navigation files

# Installation
To run this software, you will nead Visual Studio 2015 and above and also Eigen library which can be found on their website http://eigen.tuxfamily.org/index.php?title=Main_Page. 

# Usage
Add the single point positioning files in your current project or create a new one. Make sure you include the environment file provided and and the sample datasets in your working (directory).

Steps...
  1. Create a new project or use existing project
  2. Add the single point positioning files to your project
  3. Put the environment file and observation and navigation files in the working directory.
  3. Import Eigen library.
  4. Right click on project in visual studio...then go to Project Properties > C/C++ > General > Additional Include Directories. Browse to      the folder where you saved the Eigen library. Then click apply.
  5. After all this build your project then run.

# Sample Results

coming soon.........

# License
The software is free for all to use. Disclaimer: Use it at your own risk
