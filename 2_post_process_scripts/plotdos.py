import numpy as np
import matplotlib.pyplot as plt
import __tBG_dkl as tBG

LDOS = tBG.Load_object("./LDOS")
Egrid = tBG.Load_object("./E_LDOS")
top = 0
for i in range(2773):
    top += LDOS[i]

plt.plot(Egrid,top)

bot= 0
for j in range(2774,5548):
    bot += LDOS[j]
plt.plot(Egrid,bot)

plt.show()
