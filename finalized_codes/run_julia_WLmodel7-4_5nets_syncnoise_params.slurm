#!/bin/bash
#SBATCH --job-name=5nets
#SBATCH -p compute
#SBATCH -t 2-0
#SBATCH --mem-per-cpu=5G
#SBATCH --array=0-719

module load julia/1.4.1
aid=$SLURM_ARRAY_TASK_ID

input="~/Archive/finalized_codes/parameters_5nets_all.txt"
count=0
while IFS= read -r line
do
        if [[ $count -eq $aid ]]; then
                l="$line"
                params=(`echo $l`)
                rho_eta_hetero=${params[0]}
                fname=${params[1]}
                rho_eta_fluc=${params[2]}
                r_std=${params[3]}
                d=${params[4]}
		n_com=${params[5]}
		n_diffenv=${params[6]}
                break
        fi
        count=$( expr $count + 1 )
done < "$input"

netdir=5nets_networks_N/N${n_com}

c=0
maxtime=4000
r=0.7

h=1
# create a temporary directory for this job and save the name
tempdir=$(mktemp -d ~/job/XXXXXX)
cd $tempdir

# main simulation run
outfile=Result_${fname}_${rho_eta_hetero}_${rho_eta_fluc}_${r_std}_${d}
echo StartOfFile > $outfile 
outdir=WLmodel7-4_5nets_synchnoise/WLmodel7-4_5nets_N${n_com}_e${n_diffenv}/h${h}/r${r}/c$c/maxgen$maxtime
julia --project=~/Archive/stabilityProject ~/Archive/finalized_codes/WLmodel7-4_5nets_syncnoise.jl $fname $netdir $rho_eta_hetero $d $c $rho_eta_fluc $maxtime $r_std $h $r $n_com $n_diffenv $outfile

# move the output files from the temporary folder to the specified output folder
scp $tempdir/* ~/Archive/$outdir
rm -r $tempdir

