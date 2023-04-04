# Genetic-Algorithm-for-Causeway-Modification

1. Background 

     This project prototypes a genetic algorithm (GA) coupled to a numerical estuarine circulation model of Old Tampa Bay, FL (OTB). The GA is designed to generate modifications (cut-throughs) of causeways in OTB and rank them based on their impact on modeled hydrodynamic flushing. The prototype was limited to 2, and only 2, modifications in the northwestmost causeway span. 

     The project was conceptualized and implemented by Steven D. Meyers, Ph.D. in collaboration with Mark E. Luther, Ph.D., both at the [College of Marine Science University of South Florida](https://www.usf.edu/marine-science/). 

     Support was provided by the [Tampa Bay Estuary Program](https://tbep.org). 

     The hydrodynamic model of OTB was based on the Estuarine and Coastal Ocean Model (ECOM), itself a derivative of the Princeton Ocean Model. ECOM is FORTRAN-based, so the GA was largely written in FORTRAN to faciliate compatibility. A control file written in BASH was used to enable optimal usage of the multi-CPU workstation by the GA. Text files were also used to pass information to and between subroutines.

     Files in this repository are research-level, not professional.  They have not been cleaned. NOTE: paths names in scripts are relative and will need to be updated after download.

2. Outline 

     This is a sketch of essential steps executed by code in this repository

     A. randomly generate set of unique causeway configurations
     B. assign copies of ECOM code to processors
     C. each copy reads unique causeway configuration and computes flushing
     D. compute fitness rank of each model output based on changes in hydrodynamic flushing
     E. Check if all copies are finished
     F. identify top-ranking configuration (elite) 
     G. keep highest ranking and mutate them to create new generation
     H. if elite not represented in new generation, randomly replace one member with the elite configuration
     I. repeat cycle B-H for fixed number of generations

3. Major Software Components

     This is a list of most important software files in the application

     * genmain.f: (FORTAN) executes 2A and calls modelgroup          
     * modelgroup: (BASH) the primary control file (2B)
     * Makefile: (BASH) creates the executables for Linux environment 
     * otbpt: (.exe) ECOM-based executable created by Makefile from ~44 .f files as listed in Makefile
     * main.f: (FORTRAN) primary ECOM routine
     * ftestN.f: (FOTRAN) N is 1 to 4 (so far) are versions of the fitness test (2D)
     * mutate2.f: (FORTRAN) most recent version of the mutation algorithm (2G)
     * roulette.f: (FORTRAN) randomizes decision on which solution to mutate
     * cutpasses.f: (FORTRAN) alters ECOM bathymetric grid to conform to generated causeway configuration
 
4. Contact Information
 
   Steven Meyers, [smeyers@usf.edu](mailto:smeyers@usf.edu), phone: 727-553-1188

<a rel='license' href='http://creativecommons.org/licenses/by/4.0/'><img src="https://i.creativecommons.org/l/by/4.0/88x31.png" alt="Creative Commons License" style="border-width:0"/></a>  This repository is licensed under a <a rel='license' href='http://creativecommons.org/licenses/by/4.0/'>Creative Commons Attribution 4.0 International License</a>.
