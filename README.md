README 

1. Background 

This project prototypes a genetic algorithm (GA) coupled to a numerical estuarine circulation model of Old Tampa Bay, FL (OTB). The GA is designed to generate modifications (cut-throughs) of causeways in OTB and rank them based on their impact on modeled hydrodynamic flushing. The prototype was limited to 2, and only 2, modifications in the northwestmost causeway span. 

The project was conceptualized and implemented by 
Steven D. Meyers, Ph.D. in collaboration with Mark E. Luther, Ph.D., both at the College of Marine Science
University of South Florida. 

Support was provided by the Tampa Bay Estuary Program. 

The hydrodynamic model of OTB was based on the Estuarine and Coastal Ocean Model (ECOM), itself a derivative of the Princeton Ocean Model. ECOM is FORTRAN-based, so the GA was largely written in FORTRAN to faciliate compatibility. A control file written in BASH was used to enable optimal usage of the multi-CPU workstation by the GA. Text files were also used to pass information to and between subroutines.

Files in this repository are research-level, not professional.  
They have not been cleaned. NOTE: paths names in scripts are original and will need to be updated after download.

2. Outline 

This is a sketch of essential steps executed by code in this repository

A. create random generation of unique causeway configurations
B. assign copies of ECOM code to processors
C. each copy reads own causeway configuration and 
   computes flushing as prescribed
D. fitness rank each model output by impact on flushing
E. identify top-ranking configuration (elite) 
F. keep highest ranking and mutate them to create new generation
G. if elite not represented in new generation, randomly replace 
   one member with the elite configuration
H. check if all copies are finished executing
I. repeat cycle from B for fixed number of generations

3. Major Software Components

This is a list of most important software files in the algorithm

modelgroup: (BASH) the primary control file (2B)
genmain.f: (FORTAN) executes 2A and calls modelgroup
Makefile: (BASH) creates the executables for Linux environment 
otbpt: (.exe) ECOM-based executable created by Makefile from 
        ~44 .f files as listed in Makefile
main.f: (FORTRAN) primary ECOM routine
ftestN.f: (FOTRAN) N is from 1 to 4 (so far) are options for the 
      fitness test
mutate2.f: (FORTRAN) most recent version of the mutation 
         algorithm
roulette.f: (FORTRAN) randomizes decision on which solution to 
         mutate
cutpasses.f: (FORTRAN) alters ECOM bathymetric grid for 
         each causeway configuration
 
   

4. Contact Information
   Steven Meyers
   smeyers@usf.edu 
   phone: 727-553-1188



