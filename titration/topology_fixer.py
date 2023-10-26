import glob
import re
path = glob.glob("*.top")[0]
dummy = '1 \n'

with open(path, "r") as topol:
   data = topol.read()

tracker = 0 #0 for start of file, 1 for angles, 2 for dihedrals, 3 for end

for i, section in enumerate(re.split('angles|dihedrals|system', data)):
    new_section = section
    
    if i == 1:
        new_section = section.replace(dummy, '5 \n')
    elif i == 2:
        new_section = section.replace(dummy, '9 \n')
    print(new_section)
    data = data.replace(section, new_section)

#with open(path, 'w') as file:
#    file.write(data)
