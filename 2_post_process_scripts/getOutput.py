f = open("twist.mol", "r")

for line in f:
    a = line.split()
    if ("xlo xhi" in line):
       a1x = a[1]
    if ("ylo yhi" in line):
       a2x = a[1]
    if ("xy xz yz" in line):
       a2y = a[0]

g = open("unitVectors.dat", "w")

g.write(a1x + " 0.0 0.0\n")
g.write(a2y + " " + a2x + " 0.0\n")
g.write("0.0 0.0 35.0")
 
f.close()
g.close()

def tail(filename, n=10):
    'Return the last n lines of a file'
    return deque(open(filename), n)


f2 = open("dump.minimization", "r")
g3 = open("visualizeInit.xyz", "w")
g4 = open("visualize.xyz", "w")


i = 0
numberOfLines = 0

for line in f2:
    i = i+1
    if (i==4):
        numberOfLines = int(line.split()[0])
    if (i>9 and i<=9+numberOfLines):
        g3.write(line)
    if (i>9+numberOfLines):
        break

f2.close()
g3.close()

import os

from collections import deque

#endLines = tail(f2, numberOfLines, offset=0)
endLines = tail("dump.minimization", numberOfLines+4)
i = 0
for el in endLines:
    i = i+1
    if (i==1):
        xlo_bound = float(el.split()[0])
        xhi_bound = float(el.split()[1])
        xy = float(el.split()[2])
    elif (i==2):
        ylo_bound = float(el.split()[0])
        yhi_bound = float(el.split()[1])
        xz = float(el.split()[2])
    elif (i==3):
        zlo_bound = float(el.split()[0])
        zhi_bound = float(el.split()[1])
        yz = float(el.split()[2])
    elif (i==4):
        pass
    else:
        g4.write(el)

g4.close()

g5 = open("generateWithBN.xyz", "w")
f1 = open("unitVectors.dat", "r")
f2 = open("visualize.xyz", "r")

for line in f1:
    g5.write(line)

print(numberOfLines)

g5.write("\n")
g5.write(str(numberOfLines) + "\n")

for line in f2:
    a = line.split()
    if a[1] == "1":
       g5.write("B " + str(a[2] + "   " + str(a[3]) + "  " + str(a[4]) + "\n"))
    elif a[1] == "2":
       g5.write("N " + str(a[2] + "   " + str(a[3]) + "  " + str(a[4]) + "\n"))
    elif a[1] == "3":
       g5.write("C " + str(a[2] + "   " + str(a[3]) + "  " + str(a[4]) + "\n"))


f1.close()
f2.close()
g5.close()

g5 = open("generateInitWithBN.xyz", "w")
f1 = open("unitVectors.dat", "r")
f2 = open("visualizeInit.xyz", "r")

for line in f1:
    g5.write(line)

print(numberOfLines)

g5.write("\n")
g5.write(str(numberOfLines) + "\n")

for line in f2:
    a = line.split()
    if a[1] == "1":
       g5.write("B " + str(a[2] + "   " + str(a[3]) + "  " + str(a[4]) + "\n"))
    elif a[1] == "2":
       g5.write("N " + str(a[2] + "   " + str(a[3]) + "  " + str(a[4]) + "\n"))
    elif a[1] == "3":
       g5.write("C " + str(a[2] + "   " + str(a[3]) + "  " + str(a[4]) + "\n"))

f1.close()
f2.close()
g5.close()
