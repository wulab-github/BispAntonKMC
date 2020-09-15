# BispAntonKMC
This program uses rigid-body-based Monte-Carlo simulations to study the binding kinetics of MHC/PD-L1 bispecific biologics
with membrane receptors TCR and PD-1 on the 2D cell surface. The rotational and translational variations derived from Anton2 all-atom
molecular dynamic simulations are used as input to regulate the conformational fluctuations of biologics in the Monte-Carlo simulations.

The input files are listed as follows:

    1: GS15_RtTotFl_ProbDistribution.dat  >>>> the probability parameters of rotational variations for linker GS15
    
    2: GS15_TlTotFl_ProbDistribution.dat  >>>> the probability parameters of translational variations for linker GS15

    3: PLPII15_RtTotFl_ProbDistribution.dat  >>>> the probability parameters of rotational variations for linker PLP15
    
    4: PLPII15_TlTotFl_ProbDistribution.dat  >>>> the probability parameters of translational variations for linker PLP15
    
    5: PLPII30_RtTotFl_ProbDistribution.dat  >>>> the probability parameters of rotational variations for linker PLP30
    
    6: PLPII30_TlTotFl_ProbDistribution.dat  >>>> the probability parameters of translational variations for linker PLP30
    
   
The source codes that use the kinetic Monte-Carlo simulation algorithm with rigid-body representation were written in the Fortran77 format.
They are listed as follows:


    1: BispecificRBKMC_main.f  >>>>  the main program for Monte-Carlo simulation using parameters listed above as inputs
    
    2: sub_rot_along_axis.f   >>>>  the subroutine for 2D rotation of cell surface proteins along the membrane normal
    

This program also provides the compiled excutable file: BispRBKMC

The outputs from the program is also provided, using the linker GS15 as a test system:    

    1: BispRBKMCresult_simupara_SamplTry_TgCel_GS15_PDWT0.txt  >>>>  the list of all the simulation parameters including the diffusion constants, binding rates and concentration
    
    2: BispRBKMCresult_ene_SamplTry_TgCel_GS15_PDWT0.dat   >>>> the records of interactions formed along the first 20000 simulation steps
    
    3: BispRBKMCresult_trj_SamplTry_TgCel_GS15_PDWT0.pdb   >>>> the snapshots from the simulation with the PDB format
    
In order to test other linker systems, one need to change the input files in the source codes.
    
The format of the output file with the records of interactions is shown as following:

         100     0     0     0
         200     0     0     0
         300     0     0     0
         400     0     0     0
         500     0     0     0
         600     0     0     0
         700     0     0     0
         800     0     0     0
         900     0     0     0
        1000     0     0     0
        
    The 1st colume is the index of simulaiton steps
    The 2nd colume is the number of interactions formed between TCR and MHC
    The 3rd colume is the number of interactions formed between PD-1 and PD-L1
    The 4th colume is the nubmer of biologics interacting with both TCR and PD-1
