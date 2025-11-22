import numpy as np

outKey = "322"

n_levels = [2, 3]
# ignore_lvls = ["311"]
ignore_lvls = []

S1 = []
S2 = []
S3 = []
S4 = []

file1 = open("load_rate_input.txt", 'r')
Lines = file1.readlines()
  
final_out = r"\dot{\epsilon}_{"+outKey+"} = \gamma_{"+outKey+ r"}^{\rm SR} \epsilon_{" + outKey + "} "
count = 0

for line in Lines:
    line_parts = line.split()
    
    if (line_parts[0] in ignore_lvls)or(line_parts[1] in ignore_lvls)or(line_parts[2] in ignore_lvls):
        continue
    
    if (float(line_parts[0][0]) in n_levels)and(float(line_parts[1][0]) in n_levels)and(float(line_parts[2][0]) in n_levels):
        
        
        
        signL = "-"
        preF = " "
        if (outKey != line_parts[0])and(outKey != line_parts[1])and(outKey != line_parts[2]):
            continue
            
        if outKey == line_parts[2]:
            signL = "+"
        
        if (line_parts[0] == line_parts[1])and(line_parts[0]==outKey):
            preF = "2"
            epsPre = r"\epsilon_{" + line_parts[0] + r"}^2 "
        else:
            epsPre = r"\epsilon_{" + line_parts[0] + r"} \epsilon_{" + line_parts[1] + "} "
            
        fnlstate = r"{\rm BH}"
        if line_parts[3] == "Inf":
            fnlstate = "\infty"
        
        final_out += signL + preF + "\gamma_{" + line_parts[0] + r"\times " + line_parts[1] + "}^{" + line_parts[2] + r" \times " + fnlstate + "} " + epsPre + "\epsilon_{" + line_parts[2] + "} "
        
        
print(final_out)
