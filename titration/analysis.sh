#!/bin/bash --login

#$ -N ola_NAME_VAR_analysis_REP_VAR
#$ -hold_jid ola_NAME_VAR_REP_VAR
#$ -cwd

#Run first analysis script
python titration.py
