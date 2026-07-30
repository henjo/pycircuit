import re

with open('pycircuit/circuit/jaxtransient.py', 'r') as f:
    lines = f.readlines()

new_lines = []
unindent = False
for line in lines:
    if line.startswith('        def cond_fun(nr_state: NewtonState):'):
        unindent = True
        
    if unindent:
        if line.startswith('    '):
            new_lines.append(line[4:])
        else:
            new_lines.append(line)
    else:
        new_lines.append(line)
        
with open('pycircuit/circuit/jaxtransient.py', 'w') as f:
    f.writelines(new_lines)
