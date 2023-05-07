f = open("twist.mol","r").readlines()
g = open("twist2.mol","w")

for i in range(len(f)):
    line = f[i]
    v = line.split()
    if i <= 15:
       g.write(line)
    else:
        
        if v[2] == '1':
          g.write(v[0]+' '+v[1]+' '+ v[2] +' 0.82275 '+v[3]+' '+v[4]+' '+v[5]+' '+ v[6]+' ' + v[7]+' '+v[8] +'\n')
        elif v[2] == '2':
          g.write(v[0]+' '+v[1]+' '+ v[2] +' -0.82275 '+v[3]+' '+v[4]+' '+v[5]+' '+ v[6]+' '+ v[7]+' '+v[8] +'\n')
