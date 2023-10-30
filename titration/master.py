########### Import modules ###########

import numpy as np
import os
import re
import shutil

########### Read input file ###########

main_sec = False

with open("input.in", "r") as input_file:
    for line in input_file:
        line = line[:-1] 
        #Check if in MAIN
        if re.match("^MAIN", line):
            main_sec = True

        elif main_sec == True:
                    if re.match('^Repeats', line, re.IGNORECASE):
                        num_iter = int(re.split('=|#', line)[1])
                        break

########### Generate initial directories for each repeat ###########

for iteration in range(1, num_iter+1):
    dirname = f"iter{iteration}"
    if os.path.exists(dirname):
            shutil.rmtree(dirname)
    os.mkdir(dirname)
    shutil.copytree("templates", f"{dirname}/setup")