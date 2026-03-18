import os
import re

# Use CWD-relative paths: callers always cd to tests/sage-model-tests/ before invoking
input_path = 'input/millennium.par'
output_path = 'mypar_files/millennium_settings.txt'

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
