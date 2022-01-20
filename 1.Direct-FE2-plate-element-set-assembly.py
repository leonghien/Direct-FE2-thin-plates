# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 19:47:24 2018

@author: a0132600
"""

# -*- coding: utf-8 -*-


ele_connect=raw_input("name of element_ connectivity .dat file?")
FE_2_model_gen=raw_input("name of Fe_square abaqus model generation .py file?")
integration=raw_input("input 0 for full integration or 1 for reduced integration FE_2 element:")
order=raw_input("Enter 1 for linear super elements or Enter 2 for quadratic super elements:")
import string
import math
import sys
import os
pi=3.1415926

Nodes=[];
nodal_connectvity=[];

GPs=([[-3**-0.5,-3**-0.5,0],[3**-0.5,-3**-0.5,0],[3**-0.5,3**-0.5,0],[-3**-0.5,3**-0.5,0]],#Full
[[0,0,0]])# Reduced 

os.chdir('C:\Users\zyliu\Desktop\FE2_Paper\Submit file final')
#C:\Temp\Karthic_temp\FE_2_FYP\cantilever_FE_2_generation\1x10_lin.txt
keyword=['*Node','*Element, type=S4'];
f1=open('RPs-definition.txt','r+')#input file name
inp=[]
#////////////////////////    read in all the lines in the input file
while 1:
	line=f1.readline()
	if not line:
		break
	line=line.strip()
	inp.append(line)
#///////////////////////     serach the keyword in inp, and assign the variables
def Search_simple(inp,line,array):
	array[:]=[]
	a=inp.index(line)
	
	for line_temp in inp[(a+1):]:
		if line_temp=='':
			break
		line_temp=(line_temp.replace(',','')).split()
		array.append(line_temp)
        
#////////
def Search(inp,line,array):
    array[:]=[]
    a=inp.index(line);
    
    if (inp[a].count("generate")==0):
        #print('if')
        for line_temp in inp[(a+1):]:
            if (line_temp=='') or (line_temp.count("*")!=0):
                    break
            line_temp=(line_temp.replace(',','')).split()
            array.append(line_temp)
    else:
        #print('else')
        for line_temp in inp[(a+1):]:
            if (line_temp=='') or (line_temp.count("*")!=0):
                    break
            line_temp=(line_temp.replace(',','')).split()
            for i in range(int(line_temp[0]),(int(line_temp[1])+1),int(line_temp[2])):   
                array.append(i)
#////
def tsi_eta_shape_fn(tsi,eta):
    N1=float(0.25*(1-tsi)*(1-eta))
    N2=float(0.25*(1+tsi)*(1-eta))
    N3=float(0.25*(1+tsi)*(1+eta))
    N4=float(0.25*(1-tsi)*(1+eta))
    return N1,N2,N3,N4

def tsi_eta_quad_shape_fn(tsi,eta):
    N1=float(-0.25*(1-tsi)*(1-eta)*(1+tsi+eta))
    N2=float(-0.25*(1+tsi)*(1-eta)*(1-tsi+eta))
    N3=float(-0.25*(1+tsi)*(1+eta)*(1-tsi-eta))
    N4=float(-0.25*(1-tsi)*(1+eta)*(1+tsi-eta))
    N5=float(0.25*2*(1-tsi**2)*(1-eta))
    N6=float(0.25*2*(1+tsi)*(1-eta**2))
    N7=float(0.25*2*(1-tsi**2)*(1+eta))
    N8=float(0.25*2*(1-tsi)*(1-eta**2)) 
    return N1,N2,N3,N4,N5,N6,N7,N8  

Search(inp,keyword[0],Nodes);#node no, x,y,x
Search(inp,keyword[1],nodal_connectvity);#element no,N1,N2,N3,N4

model=raw_input("Which model?")
Partname=raw_input("Which part?")

RVE=open('%s.py'%(FE_2_model_gen),'w')

print>>RVE, "from abaqus import *"
print>>RVE, "from abaqusConstants import *"
a='''
#==========import modulus of abaqus
import part
import regionToolset
import displayGroupMdbToolset as dgm
import material
import section
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import displayGroupMdbToolset as dgm
'''
RVE.write(a)

#=====================   construction of the RVE cell
a='''
import numpy as np
import math
a1 = mdb.models['%s'].rootAssembly
'''%(model)
RVE.write(a)
for i in range(0,len(Nodes)):
    print>>RVE,'a1.ReferencePoint(point=(%E, %E, %E))'%(float(Nodes[i][1]),float(Nodes[i][2]),float(Nodes[i][3]))


for i in range(0,len(nodal_connectvity)):
    ele_no=int(nodal_connectvity[i][0]);
    
    n1=int(nodal_connectvity[i][3])
    n2=int(nodal_connectvity[i][4])
    n3=int(nodal_connectvity[i][1])
    n4=int(nodal_connectvity[i][2])
    
    X1=float(Nodes[n1-1][1])
    Y1=float(Nodes[n1-1][2])
    Z1=float(Nodes[n1-1][3])
    
    X2=float(Nodes[n2-1][1])
    Y2=float(Nodes[n2-1][2])
    Z2=float(Nodes[n2-1][3])

    X3=float(Nodes[n3-1][1])
    Y3=float(Nodes[n3-1][2])
    Z3=float(Nodes[n3-1][3])  
    
    X4=float(Nodes[n4-1][1])
    Y4=float(Nodes[n4-1][2])
    Z4=float(Nodes[n4-1][3])
    
    if order==2:
        n5=int(nodal_connectvity[i][5])
        n6=int(nodal_connectvity[i][6])
        n7=int(nodal_connectvity[i][7])
        n8=int(nodal_connectvity[i][8])
        
        X5=float(Nodes[n5-1][1])
        Y5=float(Nodes[n5-1][2])
        Z5=float(Nodes[n5-1][3])    
    
        X6=float(Nodes[n6-1][1])
        Y6=float(Nodes[n6-1][2])
        Z6=float(Nodes[n6-1][3])
    
        X7=float(Nodes[n7-1][1])
        Y7=float(Nodes[n7-1][2])
        Z7=float(Nodes[n7-1][3])    
    
        X8=float(Nodes[n8-1][1])
        Y8=float(Nodes[n8-1][2])
        Z8=float(Nodes[n8-1][3])    
    

    GP_RVEs=GPs[int(integration)]

    print 'N=',i+1
    for ii in range(1,(len(GP_RVEs)+1)):
        N1,N2,N3,N4=tsi_eta_shape_fn(GP_RVEs[ii-1][0],GP_RVEs[ii-1][1]);
        print 'shape function',N1,N2,N3,N4
        print 'X_coordinate',X1,X2,X3,X4
        print 'Y_coordinate',Y1,Y2,Y3,Y4
        position_x=N1*X1+N2*X2+N3*X3+N4*X4
        position_y=N1*Y1+N2*Y2+N3*Y3+N4*Y4
        position_z=0
        print 'RVE_coordinate',position_x,position_y,position_z

      
       
            
        print>>RVE, '#=========When assembling  '+'Super ele-'+str(ele_no)+'-RVE @ GP-'+str(ii)
        a='''a = mdb.models['%s'].rootAssembly
#a.DatumCsysByDefault(CARTESIAN)
'''%(model)
        RVE.write(a)
        #Partname='';Partname=str(Part-1-ele-%d-GP-%d
        print>>RVE, "p = mdb.models['%s'].parts['"%(model)+str(Partname)+"']"
        
        print>>RVE, "a.Instance(name='Part-1-ele-%s-GP-%s', part=p, dependent=ON)" %(str(ele_no),str(ii))	
        print>>RVE, "p2 = a.instances['Part-1-ele-%s-GP-%s']" %(str(ele_no),str(ii))	
        a='''a = mdb.models['%s'].rootAssembly
a = mdb.models['%s'].rootAssembly
'''%(model,model)
        RVE.write(a)	
        fibre=i
        
        angle_x=0.0
        angle_y=0.0
        
        
        angle_z=0.0
        print>>RVE,'p2.rotateAboutAxis(axisPoint=(0.0,0.0,0.0),axisDirection=(1.0,0.0,0.0),angle=%s)'   %(str(angle_x))
        print>>RVE,'p2.rotateAboutAxis(axisPoint=(0.0,0.0,0.0),axisDirection=(0.0,1.0,0.0),angle=%s)'   %(str(angle_y))
        print>>RVE,'p2.rotateAboutAxis(axisPoint=(0.0,0.0,0.0),axisDirection=(0.0,0.0,1.0),angle=%s)'   %(str(angle_z))

        print>>RVE, "p2.translate(vector=(%s,%s,%s))"   %(str(position_x),str(position_y),str(position_z))
        print>>RVE, "a = mdb.models['%s'].rootAssembly"%(model)
        print>>RVE, "p2 = a.instances['Part-1-ele-%s-GP-%s']" %(str(ele_no),str(ii))
        print>>RVE,''

	
        
 
RVE.close()

NC=open('%s.dat'%(str(ele_connect)),'w') 
for i in range(0,len(nodal_connectvity)):
    if int(order)==1:#linear
        print>>NC,'[%d, %d, %d, %d],'%(int(nodal_connectvity[i][3]),int(nodal_connectvity[i][4]),
                   int(nodal_connectvity[i][1]),int(nodal_connectvity[i][2]))
    if int(order)==2:#quadratic
        print>>NC,'[%d, %d, %d, %d, %d, %d, %d, %d],'%(int(nodal_connectvity[i][1]),int(nodal_connectvity[i][2]),
                   int(nodal_connectvity[i][3]),int(nodal_connectvity[i][4]),int(nodal_connectvity[i][5]),
                   int(nodal_connectvity[i][6]),int(nodal_connectvity[i][7]),int(nodal_connectvity[i][8]))

NC.close()

