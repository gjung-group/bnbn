#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
top = open('top.xyz','r')
tN_pos = open('tN','w')
tB_pos = open('tB','w')
lines = top.readlines()
type(lines)
for i in range(0,len(lines)):
    if (i%2!=0):
        tN_pos.write(lines[i])
    else:
        tB_pos.write(lines[i])
    
top.close()
tN_pos.close()
tB_pos.close()

bot = open('bot.xyz','r')
bN_pos = open('bN','w')
bB_pos = open('bB','w')
lines = bot.readlines()
type(lines)
for i in range(0,len(lines)):
    if (i%2!=0):
        bN_pos.write(lines[i])
    else:
        bB_pos.write(lines[i])
    
bot.close()
bN_pos.close()
bB_pos.close()


# In[ ]:




