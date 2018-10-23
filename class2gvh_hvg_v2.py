# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 10:29:19 2017
name: class2gvh_hvg.py
@author: changliu
This script reads in output from HLAmatchmaker to quantify the number of mismatched epitopes in the graft-versus-host (GVH) and host-versus-graft (HVG) directions respectively.
input files: class2 analysis result from Hmm, tab 4 and tab 5 exported as class2rec.csv and class2don.csv respective.
output: class2gvh_hvg --> need to pivot using the type of ep-mm.
Note: 
1. no need to delete all the blank rows at the end of each input file. 
2. NEW: dose of mm is considered with this algorithm. (if a host-specific ep appears twice, counted twice)
3. NEW: ep crossreactive between two loci is not counted (for example 26L is shared by DQ and DR)
*4. for the don file: DRab=[16:41]+[76:101], DRot=[41:76]+[101:136], DQBab=[256:286]+[315:345], DQBot=[286:315]+[345:374]
*5. for the rec file: DRab=[17:42]+[77:102] , DRot=[42:77]+[102:137], DQBab=[257:287]+[316:346], DQBot=[287:316]+[346:375]
*6. ABC are gene names, "ab" means antibody confirmed epitope, "ot" means other.
Validation: 
Categorical counts add up to the total ep mm count in each direction. 
HvG directed mmEp count correlates very well between original Hmm and modified Hmm. correlation coefficient = 0.9960510417
(CORREL function), as a measure of the linear correlation of the two data sets.
"""

from __future__ import division
import sys
import os
import re

don=open('class2don.csv','r')
rec=open('class2rec.csv','r')
out=open('class2gvh_hvg','w')
out.write('type of ep mm'+'\t'+'case row#'+'\t'+'Ep mm count'+'\n') #column titles

DDRab=[]
DDRot=[]
DDQab=[]
DDQot=[]
DonEp=[]
for line in don:
    if line.split(',')[3].startswith('DRB1'):
        line=line.upper()
        DonEp.append(line.split(',')[16:])
        DDRab.append(line.split(',')[16:41]+line.split(',')[76:101])
        DDRot.append(line.split(',')[41:76]+line.split(',')[101:136])
        DDQab.append(line.split(',')[256:286]+line.split(',')[315:345])
        DDQot.append(line.split(',')[286:315]+line.split(',')[345:374])
        
RDRab=[]
RDRot=[]
RDQab=[]
RDQot=[]
RecEp=[]
for line in rec: 
    if line.split(',')[4].startswith('DRB1'): 
        line=line.upper()
        RecEp.append(line.split(',')[17:])
        RDRab.append(line.split(',')[17:42]+line.split(',')[77:102])
        RDRot.append(line.split(',')[42:77]+line.split(',')[102:137])
        RDQab.append(line.split(',')[257:287]+line.split(',')[316:346])
        RDQot.append(line.split(',')[287:316]+line.split(',')[346:375])

catego=['DRab','DRot','DQab','DQot']

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
for REpType in [RDRab,RDRot,RDQab,RDQot]: #four blocks of results DRab --> DRot --> DQab --> DQot
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
for DEpType in [DDRab,DDRot,DDQab,DDQot]: #four blocks of results DRab --> DRot --> DQab --> DQot
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
