# Codes #
In the codes, the path to the "Archive" folder that contains this "finalized_codes" are written as "~/Archive". Change adjust this path to match your folder organization.
A run of each code contains 100 stochastic simulations. Each result out of the 100 stochastic simulations will be added as a line in the output file.
1-1. WLmodel7-4_5nets.jl: Main simulation with default setting; competition coefficient c = constant (between 0 and 1), the environment in each patch fluctuates differently, though they are spatially autocorrelated (rho_fluc determins the level of spatial autocorrelation in this environmental noise).
1-2. WLmodel7-4_5nets_randomc.jl: Use a 2-d random matrix (C) to define competition for resource between species; c_ij is sampled from a uniform distribution between 0 and 1, when c=asymmetric (specified in the bash script below), whereas when c=symmetric, this 2-d matrix will be modified into a symmetric matrix.
1-3. WLmodel7-4_5nets_syncnoise.jl: Instead of function 'noises' in the main model 1-1, use 'synchnoises' or 'synchnoises2' to to allow environments in some patches to fluctuates in the same way over time. In 'synchnoises2', the patches that have the same mean environmental conditions fluctuates (defined as environmental file (env_rhoeta....txt) in environments folder) in the same way over time.
1-4. WLmodel7-4_5nets_syncnoise_randomc.jl: Combination of 1-2 and 1-3.
1-5. environment.py: a code to generate environmnetal conditions for n_com patches, out of which n_diffenv number of patches have different number of mean environmental conditions. rho_eta: a paremeter for spatial autocorrelation of the mean environmental conditions. I ran this code (with random numbers) for n_com = 25, 50, 100, 500, 750, and n_diffenv = n_com or 25, until it generates a file with environmental conditions for each of these parameters, which is stored in Archive/environments.


# Bash scripts to run the corresponding codes above #
In the following scripts, aid=$SLURM_ARRAY_TASK_ID is for array job in slurm system. Replace this with aid=[a line number] (counting from 0) to specify which parameter set in the file parameters_5nets_all.txt to use. Or simply define parameters (rho_eta_hetero, fname, rho_eta_fluc,r_std,d,n_com,n_diffenv) in the script and run the code.
2-1. run_julia_WLmodel7-4_5nets_params.slurm
2-2. run_julia_WLmodel7-4_5nets_randomc_params.slurm
2-3. run_julia_WLmodel7-4_5nets_syncnoise_params.slurm
2-4. run_julia_WLmodel7-4_5nets_syncnoise_randomc_params.slurm
2-5. (environment.py runs on its own, with input arguments n_com and n_diffenv.)
