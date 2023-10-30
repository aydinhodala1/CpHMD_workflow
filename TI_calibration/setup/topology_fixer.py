import glob
import re

path = glob.glob("*.top")[0]
dummy = '1 \n'

with open(path, "r") as topol:
   data = topol.read()

tracker = 0 #0 for start, 1 for atoms, 2 for bonds, 3 for angles, 4 for dihedrals, 5 for end

for i, section in enumerate(re.split('atoms|bonds|angles|dihedrals|system', data)):
    new_section = section

    if i == 1:
        split_sec = section.split('\n')
        for line in split_sec:
            new_line = line
            try:

                if (int(line.split()[0])-1) % 54 == 35 or (int(line.split()[0]) -1) % 54 == 36:
                    new_line = line.replace('HA2', 'HE1')
                    print(new_line)
                elif (int(line.split()[0]) -1) % 54 == 51 or (int(line.split()[0]) -1) % 54 == 52 or (int(line.split()[0])-1) %54 == 53:
                    new_line = line.replace('HA2', 'HA3')
                    print(new_line)
                new_section = new_section.replace(line, new_line)

            except:
                pass

    elif i == 3:
        new_section = section.replace(dummy, '5 \n')

    elif i == 4:
        new_section = section.replace(dummy, '9 \n')

    data = data.replace(section, new_section)

with open(path, 'w') as file:
    file.write(data)
