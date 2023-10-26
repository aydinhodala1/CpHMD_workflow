#!/bin/bash --login

#$ -N ola_NAME_VAR_fin_analysis
#$ -hold_jid ola_NAME_VAR_analysis_*
#$ -cwd

#Run final analysis script
python final_titration.py
