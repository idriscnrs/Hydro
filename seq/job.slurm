#!/bin/bash
#SBATCH --job-name=Hydro                # Nom du job
#SBATCH --ntasks=1                      # Nombre de processus MPI
#SBATCH --gres=gpu:1                    # Nombre de GPUs
# /!\ Attention, la ligne suivante est trompeuse mais dans le vocabulaire
# de Slurm "multithread" fait bien référence à l'hyperthreading.
#SBATCH --hint=nomultithread            # 1 thread par coeur physique (pas d'hyperthreading)
#SBATCH --time=00:30:00                 # Temps d’exécution maximum demandé (HH:MM:SS)
#SBATCH --output=Hydro.%j.out     # Nom du fichier de sortie
#SBATCH --error=Hydro.%j.out      # Nom du fichier d'erreur (ici commun avec la sortie)
#SBATCH -A for@gpu
#SBATCH --qos=qos_gpu-dev 
 
# on se place dans le répertoire de soumission
cd ${SLURM_SUBMIT_DIR}
 
# nettoyage des modules charges en interactif et herites par defaut
module purge
 
# chargement des modules
module load linaro-forge/23.0.4
module load module load cuda/11.2

# 
export ALLINEA_MPI_INIT=main
export ALLINEA_HOLD_MPI_INIT=1

# echo des commandes lancées
set -x
 
# exécution du code
ddt --connect srun ./hydro input_sedov_noio_250x250.nml

