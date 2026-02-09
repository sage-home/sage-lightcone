import os
import re

script_dir = os.path.dirname(os.path.abspath(__file__))
# script_dir is root/tests/sage-model-tests/utils
# We want to go to root/tests/sage-model-tests/input/millennium.par
input_path = os.path.join(script_dir, '../input/millennium.par')
# We want to output to root/tests/sage-model-tests/mypar_files/millennium_settings.txt
output_path = os.path.join(script_dir, '../mypar_files/millennium_settings.txt')

# Ensure the output directory exists
os.makedirs(os.path.dirname(output_path), exist_ok=True)

found_marker = False
with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
    for line in infile:
        if found_marker:
            if line.strip().startswith('SimulationDir') or line.strip().startswith('FileWithSnapList'):
                 # Match any path ending in /input/ and replace it with ./input/
                 line = re.sub(r' +/\S+/input/', ' ./input/', line)
            elif line.strip().startswith('OutputDir'):
                 # Match any path ending in /output/ and replace it with ./output/
                 line = re.sub(r' +/\S+/output/', ' ./output/', line)
            
            outfile.write(line)
        elif 'OutputFormat' in line:
            found_marker = True
