#!/bin/bash
#SBATCH --job-name=5nets
#SBATCH -p compute
#SBATCH -t 1-0
#SBATCH --mem-per-cpu=5G

#aid=${list[$SLURM_ARRAY_TASK_ID]}
aid=$SLURM_ARRAY_TASK_ID

rho_eta_hetero=0.8
fname=Complete_undir
rho_eta_fluc=0.95
r_std=1.0
d=0.5
n_com=25
n_diffenv=$n_com

netdir=5nets_networks_N/N${n_com}

module load julia/1.4.1

cseed=0 #0 or 10
h=0.1
c=asymmetric

maxtime=4000
r=0.7
# create a temporary directory for this job and save the name
tempdir=$(mktemp -d ~/job/XXXXXX)
# enter the temporary directory
cd $tempdir
outfile=sp5_${fname}_${rho_eta_hetero}_${rho_eta_fluc}_${r_std}_${d}
echo StartOfFile > $outfile 
outdir=WLmodel7-4_5nets_N${n_com}_e${n_diffenv}/h${h}/r${r}/${c}/randc_seed${cseed}/maxgen$maxtime
julia --project=/home/y/yuka-suzuki/stability/stabilityProject ~/Archive/finallized_codes/WLmodel7-4_5nets_randomc_plots.jl $fname $netdir $rho_eta_hetero $d $c $rho_eta_fluc $maxtime $r_std $h $r $n_com $n_diffenv $outfile $cseed

scp $tempdir/* ~/Archive/$outdir
rm -r $tempdir

