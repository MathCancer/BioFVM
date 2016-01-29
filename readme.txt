BioFVM: Finite Volume Solver for Biological Problems. 

Version:      1.0.4
Release date: 28 January 2016

Homepage:     http://BioFVM.MathCancer.org
Downloads:    http://BioFVM.sf.net
Mirror:       https://github.com/MathCancer/BioFVM
			  
Summary: 
This update removes a mistaken file from the main directory 
and tests github as a project mirror  

New features:
none 

Bugfixes: 
+ removed a convergence source file from the main directory 


Version:      1.0.3
Release date: 21 January 2016

Homepage:     http://BioFVM.MathCancer.org
Downloads:    http://BioFVM.sf.net

Summary: 
This update contains small bug fixes for BioFVM, a template 
project (a C++ with all the syntax of starting a BioFVM 
application and a sample makefile), and matlab files to 
facilitate visualization and  

New features:
+ Added template_project.cpp and Makefile in the new 
  template_project subdirectory. Copy this template when 
  making a new BioFVM project.
  
+ Added read_microenvironment.m to the matlab subdirectory 
  to read BioFVM data and arrange it into intutitive arrays.

+ Added plot_microenvironment.m to the matlab subdirectory 
  to display filled contour maps of BioFVM data.

Bugfixes: 
+ Updated Microenvironment::display_information() to fix 
  minor display bugs and typos. 
