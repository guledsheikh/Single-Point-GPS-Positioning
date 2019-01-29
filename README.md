# Introduction
Accurate determination of position has become very important in recent years. To estimate the GNSS receiver’s position at any time, one determines the four unknown parameters, that is, the X, Y, Z coordinates of the receiver and its clock error. 

This project focuses on the implementation of static single point positioning (SPP). In static SPP, the receiver is stationary throughout the GPS data collection period. 

A software module for conducting SPP using Least Squares approximation is given here using GPS pseudorange measurements. The software utility includes the following main functions which needs to be added to ones project for succesful execution:
  1. SPP.cpp
  2. Satelliteposition.cpp
  3. Reader.cpp

# Installation
You nead Eigen library which can be found on their website http://eigen.tuxfamily.org/index.php?title=Main_Page. 
