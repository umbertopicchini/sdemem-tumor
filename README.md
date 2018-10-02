# SDE and ODE mixed-effects models for tumor growth in mice

This is accompanying code for Picchini, U. and Forman J.L. "Bayesian inference for stochastic differential equation mixed effects models of a tumor xenography study", arXiv:1607.02633

Contact Umberto Picchini for info. Currently at https://umbertopicchini.github.io/ or just search me on the internet using the usual tools.

-----------

The code performs inference for a specific type of SDE mixed effects models (Matlab code) and ODE mixed-effects models (Stan code, via the Rstan interface). In the following
we use the acronyms SDEMEMs and ODEMEMs to denote the two types of models above.

The inference methods considered are:
- Bayesian synthetic likelihoods (BSL) for SDEMEMs.
- particle marginal MCMC (PMCMC) using the auxiliary particle filter (APF), for SDEMEMs.
- exact Bayesian inference for ODEMEMs.

Code is organised as follows:

- folder "utilities": this folder should be loaded into Matlab's search path before running any script/run file. It also contain the data file "tumorlong_nofollowup.txt".

- folder Synlike: contains code using BSL for the three groups of subjects considered in the paper; it contains corresponding subfolders "group_x", where x=1,3,5.
                  There is considerable redundancy of files in the three subfolders, meaning that copies of some specific files are unnecessarily copied in all the group_x folders.
                  However some files having the same names are actually quite different,  as for example the underlying model to fit group 5 is different from the one used for group 1 and 3,
                  so files with the same name but in different subfolders might actually behave differently.
                  Therefore, as a safety measure, it is better to keep things this way.

- folder PMCMC-APF: contains code for particle MCMC (PMCMC) using the auxiliary particle filter (APF). Same as above, there is some redundancy.

- folder ODEMEM-stan: produces exact Bayesian inference for ODEMEMs using Rstan, the R interface to Stan http://mc-stan.org/.
   

** INSTALLATION:

- load the folder "utilities" into Matlab's search path.
- the file mytruncgaussrandraw.m contains a call to a Matlab MEX-file normcdf_fex. Unless Windows is your OS, the file normcdf_fex.c should be compiled for your architecture (run "mex -setup" for more info).
  If you are on a Windows machine you should be fine as I have included compiled versions for x32 and x64 architectures.
  The MEX file is used for computational efficiency. However if the compilation of the C-file normcdf_fex.c causes you headache, you can replace
  the lines 75-76 in mytruncgaussrandraw.m, that is

     PHIl = normcdf_fex((a-mu)/sigma);  
     PHIr = normcdf_fex((b-mu)/sigma); 

  with the equivalent
 
     PHIl = normcdf((a-mu)/sigma);
     PHIr = normcdf((b-mu)/sigma);

** USAGE: 

For SDEMEM inference: access one of the subfolders into the Synlike or PMCMC-APF folders and run the corresponding tumor_run file.
For ODEMEM inference: acces folder ODEMEM-stan and run fitdata1.R (fits group 1) or fitdata3.R (fits group 3).

** Info on the data-files

The data file used for SDEMEM inference in Matlab is "tumorlong_nofollowup.txt", contained in the "utilities" folder. The first column report the day of the measurement. The second column is the ID for each mouse. The third column is the tumor volume [mm^3]. The fourth column is the group the mouse is assigned to. There are 5 groups and we only analysed mice in groups 1, 3 and 5. Also, as explained in the paper, some mice have been discarded from the analyses, these being mice from group 1 having ID 11,13,15,17,18. As mentioned in the paper, with the exception of mice in group 5 we have considered data from day 6 onward.    

Folder ODEMEM-stan contains also two datafiles (datagroup1 and datagroup3) which are subsets of tumorlong_nofollowup.txt, but also differ from tumorlong_nofollowup.txt in the following ways: (i) here the first column gives "normalised" days starting from day 6. However, we normalize those as described in section 5 of the paper, where times are divided by the largest observed time across all groups (example: the largest observed time is day 39; hence since we consider measurements from day 6, we have that for each subject the first reported time is 6/39 = 0.1539).  (ii) the second column gives the natural logarithm of the volumes (instead of the volumes). (iii) the third column is the mouse ID. (iv) Measurements are from day 6 onward (i.e. from time 6/39=0.1539).

Data was kindly provided by the research team at the Center for Nanomedicine and Theranostics (DTU Nanotech, Denmark).

