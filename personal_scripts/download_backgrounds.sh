#!/bin/bash
#SBATCH --job-name=background_download
#SBATCH --output=log_files/background_download_%j.out
#SBATCH --error=log_files/background_download_%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=15G

source ~/.bashrc

cd /storage/group/izg5139/default/lefteris/background_pop

aria2c -i ftp_files.txt \
    -j 10 \
    -x 4 \
    -c \
    -d . \
    --max-tries=0 \
    --retry-wait=10 \
    --ftp-pasv=true