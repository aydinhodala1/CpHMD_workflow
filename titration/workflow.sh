#!/bin/bash --login

for repeat in {1..3}
do
	#File generation
	mkdir iter$repeat
	cp -r setup ./iter$repeat 
	sed -i s/REP_VAR/$repeat/g ./iter$repeat/setup/prep.sh
	sed s/REP_VAR/$repeat/g production.sh > ./iter$repeat/production.sh
        sed s/REP_VAR/$repeat/g analysis.sh > ./iter$repeat/ analysis.sh
        cp titration.py ./iter$repeat
	
	#System generation
	cd iter$repeat/setup
	qsub prep.sh
	cd ..
	
	#Production runs
	mkdir data
	qsub production.sh
        
        #Analysis jobs
        qsub analysis.sh
	cd ..
done

#Final analysis job
qsub final_analysis.sh
