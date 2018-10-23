# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 10:29:19 2017
name: class1gvh_hvg.py
@author: changliu
This script reads in output from HLAmatchmaker to quantify the number of mismatched epitopes in the graft-versus-host (GVH) and host-versus-graft (HVG) directions respectively.
input files: class1 analysis result from Hmm, tab 4 and tab 5 exported as class1rec.csv and class1don.csv respective.
output: class1gvh_hvg --> need to pivot using the type of ep-mm.
Note: 
1. no need to delete all the blank rows at the end of each input file. 
2. NEW: dose of mm is considered with this algorithm. (if a host-specific ep appears twice, counted twice)
3. NEW: ep crossreactive between two loci is not counted (for example 26L is shared by DQ and DR)
*4. for the don file: Aab=[8:68]+[118:178], Aot=[68:118]+[178:228], Bab=[228:288]+[338:398], Bot=[288:338]+[398:448], Cab=[448:508]+[558:618], Cot=[508:558]+[618:668]
*5. for the rec file: Aab=[10:70]+[120:180], Aot=[70:120]+[180+230], Bab=[230:290]+[340:400], Bot=[290:340]+[400:450], Cab=[450:510]+[560:620], Cot=[510:560]+[620:670]
*6. ABC are gene names, "ab" means antibody confirmed epitope, "ot" means other. 
Validation: 
Categorical counts add up to the total ep mm count in each direction. 
HvG directed mmEp count correlates very well between original Hmm and modified Hmm. correlation coefficient = 0.9861282253
(CORREL function), as a measure of the linear correlation of the two data sets.
"""

from __future__ import division
import sys
import os
import re

don=open('class1don.csv','r')
rec=open('class1rec.csv','r')
out=open('class1gvh_hvg','w')
out.write('type of ep mm'+'\t'+'case row#'+'\t'+'Ep mm count'+'\n') #column titles

DAab=[]
DAot=[]
DBab=[]
DBot=[]
DCab=[]
DCot=[]
DonEp=[]
for line in don:
    if line.split(',')[2].startswith('A'):
        line=line.upper()
        DonEp.append(line.split(',')[8:668])
        DAab.append(line.split(',')[8:68]+line.split(',')[118:178])
        DAot.append(line.split(',')[68:118]+line.split(',')[178:228])
        DBab.append(line.split(',')[228:288]+line.split(',')[338:398])
        DBot.append(line.split(',')[288:338]+line.split(',')[398:448])
        DCab.append(line.split(',')[448:508]+line.split(',')[558:618])
        DCot.append(line.split(',')[508:558]+line.split(',')[618:668])
        
RAab=[]
RAot=[]
RBab=[]
RBot=[]
RCab=[]
RCot=[]
RecEp=[]
for line in rec: 
    if line.split(',')[4].startswith('A'): 
        line=line.upper()
        RecEp.append(line.split(',')[10:670])
        RAab.append(line.split(',')[10:70]+line.split(',')[120:180])
        RAot.append(line.split(',')[70:120]+line.split(',')[180:230])
        RBab.append(line.split(',')[230:290]+line.split(',')[340:400])
        RBot.append(line.split(',')[290:340]+line.split(',')[400:450])
        RCab.append(line.split(',')[450:510]+line.split(',')[560:620])
        RCot.append(line.split(',')[510:560]+line.split(',')[620:670])

catego=['Aab','Aot','Bab','Bot','Cab','Cot']

#TOTAL GVH MM COUNT (counting ep in recipient, but not in donor)
i=0
for case in RecEp:
    gvhcount=0
    for ep in case:
        if ep.strip() not in map(str.strip, DonEp[i]):
            gvhcount=gvhcount+1
    #column 1: type of mm, 2: case row#, 3: graft-versus-host Ep mm             
    out.write('GvH total (in host, not in donor)'+'\t'+str(i+1)+'\t'+str(gvhcount)+'\n')
    i=i+1


#DIFFERENT CATEGORIES OF GVH MM COUNT
j=0        
for REpType in [RAab,RAot,RBab,RBot,RCab,RCot]: #four blocks of results DRab --> DRot --> DQab --> DQot
    i=0
    for case in REpType:
        reccount=0
        for ep in case:
            if ep.strip() not in map(str.strip, DonEp[i]):
                reccount=reccount+1
        #columns 1: case row#, 2: HVG mm, 3: GVH mm, 4: neither
        out.write('GvH '+catego[j]+'\t'+str(i+1)+'\t'+str(reccount)+'\n') 
        i=i+1
    j=j+1

#TOTAL HVG MM COUNT (counting ep in donor, but not in recipient)
i=0
for case in DonEp:
    hvgcount=0
    for ep in case:
        if ep.strip() not in map(str.strip, RecEp[i]):
            hvgcount=hvgcount+1
    #column 1: type of mm, 2: case row#, 3: host-versus-graft Ep mm             
    out.write('HvG total (in donor, not in host)'+'\t'+str(i+1)+'\t'+str(hvgcount)+'\n')
    i=i+1


#DIFFERENT CATEGORIES OF HVG MM COUNT
j=0
for DEpType in [DAab,DAot,DBab,DBot,DCab,DCot]: #four blocks of results DRab --> DRot --> DQab --> DQot
    i=0
    for case in DEpType:
        doncount=0
        for ep in case:
            if ep.strip() not in map(str.strip, RecEp[i]):
                doncount=doncount+1
        #columns 1: case row#, 2: HVG mm, 3: GVH mm, 4: neither
        out.write('HvG '+catego[j]+'\t'+str(i+1)+'\t'+str(doncount)+'\n') 
        i=i+1
    j=j+1
        
don.close()
rec.close()
out.close()               
