import re
with open('pycircuit/circuit/_jaxtoolkit.py', 'r') as f:
    lines = f.readlines()
    
# Find where jax = True is
idx = lines.index("jax = True\n")
for i in range(idx+1, len(lines)):
    if lines[i].startswith("    "):
        lines[i] = lines[i][4:]
        
with open('pycircuit/circuit/_jaxtoolkit.py', 'w') as f:
    f.writelines(lines)
print("Fixed indent")
