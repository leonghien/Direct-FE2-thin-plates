# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# script to apply MPCs for FE2 method by using the isoparametric formulation

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

import numpy as np
import math

session.journalOptions.setValues(recoverGeometry = COORDINATE)

#Parameter


modelName,Num_NRP,Num_element,Num_Gauss,Macro_thickness,RVE_thickness= getInputs(
	fields=(('Model Name:', 'Model-1'),('Num_NRP:', '49'),('Num_element:', '36'),('Num_Gauss:', '4'),('Macro thickness :', '1'),('RVE thickness:', '1.5874010666389')),
	label='Enter Parameter',
	dialogTitle='By XUJUNHAO in NUS')


Num_NRP=int(Num_NRP)  #number of reference point
Num_element=int(Num_element)  #number of macro elements
Num_Gauss=int(Num_Gauss)+1 # full or reduce
Macro_thickness=float(Macro_thickness)   #thickness of macro element by thickness of RVE model
RVE_thickness=float(RVE_thickness)
Scale=float(Macro_thickness/RVE_thickness)

# modelName='Model-1'
macroscale_instance_name='Part-1-1'



def LabelName(index, root):
    return root + str(index)


def Get_Cube_dimension(modelName, instanceName):
    node = mdb.models[modelName].rootAssembly.instances[instanceName].nodes
    Xmin = 1e+08
    Xmax = -1e+08
    Ymin = 1e+09
    Ymax = -1e+09
    Zmin = 1e+11
    Zmax = -1e+10
    for i in range(len(node)):
        x = node[i].coordinates[0]
        y = node[i].coordinates[1]
        z = node[i].coordinates[2]
        if Xmin > x:
            Xmin = x
        elif Xmax < x:
            Xmax = x
        
        if Ymin > y:
            Ymin = y
        elif Ymax < y:
            Ymax = y
        
        if Zmin > z:
            Zmin = z
            continue
        if Zmax < z:
            Zmax = z
            continue
    
    return (Xmin, Xmax, Ymin, Ymax, Zmin, Zmax)





a=mdb.models[modelName].rootAssembly
for i in range(0,Num_NRP):
    a.Set(name='N'+str(i+1)+'-RP', nodes=(a.instances[str(macroscale_instance_name)].nodes[i:int(i+1)],))


def Horizontal_shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4):
    X=P.nodes[0].coordinates[0]
    Y=P.nodes[0].coordinates[1]
    
    Dx=float(4*X-(X1+X2+X3+X4))
    Dy=float(4*Y-(Y1+Y2+Y3+Y4))
    
    Ax=float(-X1+X2+X3-X4)
    Ay=float(-Y1-Y2+Y3+Y4)
    tsi=float(Dx/Ax);
    eta=float(Dy/Ay);
    
    assert abs(tsi)<=1 and abs(eta)<=1
    
    N1_H=float(0.25*(1-tsi)*(1-eta))
    N2_H=float(0.25*(1+tsi)*(1-eta))
    N3_H=float(0.25*(1+tsi)*(1+eta))
    N4_H=float(0.25*(1-tsi)*(1+eta))
    
    N1_H_dx=float((1/Ax)*(-1)*(1-eta))
    N2_H_dx=float((1/Ax)*(1)*(1-eta))
    N3_H_dx=float((1/Ax)*(1)*(1+eta))
    N4_H_dx=float((1/Ax)*(-1)*(1+eta))
    
    N1_H_dy=float((1/Ay)*(-1)*(1-tsi))
    N2_H_dy=float((1/Ay)*(-1)*(1+tsi))
    N3_H_dy=float((1/Ay)*(1)*(1+tsi))
    N4_H_dy=float((1/Ay)*(1)*(1-tsi))
    
    
    return N1_H,N2_H,N3_H,N4_H,N1_H_dx,N2_H_dx,N3_H_dx,N4_H_dx,N1_H_dy,N2_H_dy,N3_H_dy,N4_H_dy





def shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4):               
    X=P.nodes[0].coordinates[0]
    Y=P.nodes[0].coordinates[1]
    
    Dx=float(4*X-(X1+X2+X3+X4))
    Dy=float(4*Y-(Y1+Y2+Y3+Y4))
    
    Ax=float(-X1+X2+X3-X4)
    Ay=float(-Y1-Y2+Y3+Y4)
    tsi=float(Dx/Ax);
    eta=float(Dy/Ay);
    
    assert abs(tsi)<=1 and abs(eta)<=1
    N1=float(0.125*(1-tsi)*(1-eta)*(2-tsi-eta-tsi**2-eta**2))
    N2=float(0.125*(1+tsi)*(1-eta)*(2+tsi-eta-tsi**2-eta**2))
    N3=float(0.125*(1+tsi)*(1+eta)*(2+tsi+eta-tsi**2-eta**2))
    N4=float(0.125*(1-tsi)*(1+eta)*(2-tsi+eta-tsi**2-eta**2))
    
    NX1=float(-0.25*Ay*-1*(1-tsi)*(1-eta)*(1-eta**2)*0.125)
    NX2=float(-0.25*Ay*-1*(1+tsi)*(1-eta)*(1-eta**2)*0.125)
    NX3=float(-0.25*Ay*1*(1+tsi)*(1+eta)*(1-eta**2)*0.125)
    NX4=float(-0.25*Ay*1*(1-tsi)*(1+eta)*(1-eta**2)*0.125)
    
    NY1=float(0.25*Ax*-1*(1-tsi)*(1-eta)*(1-tsi**2)*0.125)
    NY2=float(0.25*Ax*1*(1+tsi)*(1-eta)*(1-tsi**2)*0.125)
    NY3=float(0.25*Ax*1*(1+tsi)*(1+eta)*(1-tsi**2)*0.125)
    NY4=float(0.25*Ax*-1*(1-tsi)*(1+eta)*(1-tsi**2)*0.125)
    
    return N1,N2,N3,N4,NX1,NX2,NX3,NX4,NY1,NY2,NY3,NY4


def dw_dy_shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4):             
    X=P.nodes[0].coordinates[0]
    Y=P.nodes[0].coordinates[1]
    
    Dx=float(4*X-(X1+X2+X3+X4))
    Dy=float(4*Y-(Y1+Y2+Y3+Y4))
    
    Ax=float(-X1+X2+X3-X4)
    Ay=float(-Y1-Y2+Y3+Y4)
    
    tsi=float(Dx/Ax);
    eta=float(Dy/Ay);
    assert abs(tsi)<=1 and abs(eta)<=1
    
    N1_dy=float((0.5/Ay)*(1-tsi)*((-1)*(2-tsi-eta-tsi**2-eta**2)+(1-eta)*(-1-2*eta)))
    N2_dy=float((0.5/Ay)*(1+tsi)*((-1)*(2+tsi-eta-tsi**2-eta**2)+(1-eta)*(-1-2*eta)))
    N3_dy=float((0.5/Ay)*(1+tsi)*((1)*(2+tsi+eta-tsi**2-eta**2)+(1+eta)*(1-2*eta)))
    N4_dy=float((0.5/Ay)*(1-tsi)*((1)*(2-tsi+eta-tsi**2-eta**2)+(1+eta)*(1-2*eta)))
    
    
    NX1_dy=float(-0.125*(-1)*(1-tsi)*((-1)*(1-eta**2)+(1-eta)*(-2*eta)))
    NX2_dy=float(-0.125*(-1)*(1+tsi)*((-1)*(1-eta**2)+(1-eta)*(-2*eta)))
    NX3_dy=float(-0.125*(1)*(1+tsi)*((1)*(1-eta**2)+(1+eta)*(-2*eta)))
    NX4_dy=float(-0.125*(1)*(1-tsi)*((1)*(1-eta**2)+(1+eta)*(-2*eta)))
    
    
    NY1_dy=float(0.125*(Ax/Ay)*(-1*-1)*(1-tsi)*(1-tsi**2))
    NY2_dy=float(0.125*(Ax/Ay)*(1*-1)*(1+tsi)*(1-tsi**2))
    NY3_dy=float(0.125*(Ax/Ay)*(1*1)*(1+tsi)*(1-tsi**2))
    NY4_dy=float(0.125*(Ax/Ay)*(-1*1)*(1-tsi)*(1-tsi**2))
    
    
    return N1_dy,N2_dy,N3_dy,N4_dy,NX1_dy,NX2_dy,NX3_dy,NX4_dy,NY1_dy,NY2_dy,NY3_dy,NY4_dy



def minus_dw_dx_shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4):      ##       without minus
    X=P.nodes[0].coordinates[0]
    Y=P.nodes[0].coordinates[1]
    
    Dx=float(4*X-(X1+X2+X3+X4))
    Dy=float(4*Y-(Y1+Y2+Y3+Y4))
            
    Ax=float(-X1+X2+X3-X4)
    Ay=float(-Y1-Y2+Y3+Y4)
    
    tsi=float(Dx/Ax);
    eta=float(Dy/Ay);
    assert abs(tsi)<=1 and abs(eta)<=1
    
    N1_dx=float(1*(0.5/Ax)*(1-eta)*((-1)*(2-tsi-eta-tsi**2-eta**2)+(1-tsi)*(-1-2*tsi)))
    N2_dx=float(1*(0.5/Ax)*(1-eta)*((1)*(2+tsi-eta-tsi**2-eta**2)+(1+tsi)*(1-2*tsi)))
    N3_dx=float(1*(0.5/Ax)*(1+eta)*((1)*(2+tsi+eta-tsi**2-eta**2)+(1+tsi)*(1-2*tsi)))
    N4_dx=float(1*(0.5/Ax)*(1+eta)*((-1)*(2-tsi+eta-tsi**2-eta**2)+(1-tsi)*(-1-2*tsi)))
    
    NX1_dx=float(-0.125*(Ay/Ax)*(-1*-1)*(1-eta)*(1-eta**2))
    NX2_dx=float(-0.125*(Ay/Ax)*(1*-1)*(1-eta)*(1-eta**2))
    NX3_dx=float(-0.125*(Ay/Ax)*(1*1)*(1+eta)*(1-eta**2))
    NX4_dx=float(-0.125*(Ay/Ax)*(-1*1)*(1+eta)*(1-eta**2))
    
    NY1_dx=float(0.125*(-1)*(1-eta)*((-1)*(1-tsi**2)+(1-tsi)*(-2*tsi)))
    NY2_dx=float(0.125*(1)*(1-eta)*((1)*(1-tsi**2)+(1+tsi)*(-2*tsi)))
    NY3_dx=float(0.125*(1)*(1+eta)*((1)*(1-tsi**2)+(1+tsi)*(-2*tsi)))
    NY4_dx=float(0.125*(-1)*(1+eta)*((-1)*(1-tsi**2)+(1-tsi)*(-2*tsi)))
    
    return N1_dx,N2_dx,N3_dx,N4_dx,NX1_dx,NX2_dx,NX3_dx,NX4_dx,NY1_dx,NY2_dx,NY3_dx,NY4_dx



def d2w_dx2_shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4):    
    X=P.nodes[0].coordinates[0]
    Y=P.nodes[0].coordinates[1]
    
    Dx=float(4*X-(X1+X2+X3+X4))
    Dy=float(4*Y-(Y1+Y2+Y3+Y4))
            
    Ax=float(-X1+X2+X3-X4)
    Ay=float(-Y1-Y2+Y3+Y4)
    
    tsi=float(Dx/Ax);
    eta=float(Dy/Ay);
    assert abs(tsi)<=1 and abs(eta)<=1
    
    N1_d2x=float(-12.0/(Ax**2)*(1-eta)*(-1)*tsi)
    N2_d2x=float(-12.0/(Ax**2)*(1-eta)*(1)*tsi)
    N3_d2x=float(-12.0/(Ax**2)*(1+eta)*(1)*tsi)
    N4_d2x=float(-12.0/(Ax**2)*(1+eta)*(-1)*tsi)
    
    NX1_d2x=float(0)
    NX2_d2x=float(0)
    NX3_d2x=float(0)
    NX4_d2x=float(0)
    
    NY1_d2x=float(-1/Ax*(1-eta)*(-1+3*tsi))
    NY2_d2x=float(-1/Ax*(1-eta)*(1+3*tsi))
    NY3_d2x=float(-1/Ax*(1+eta)*(1+3*tsi))
    NY4_d2x=float(-1/Ax*(1+eta)*(-1+3*tsi))
    
    return N1_d2x,N2_d2x,N3_d2x,N4_d2x,NX1_d2x,NX2_d2x,NX3_d2x,NX4_d2x,NY1_d2x,NY2_d2x,NY3_d2x,NY4_d2x



def d2w_dy2_shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4):    
    X=P.nodes[0].coordinates[0]
    Y=P.nodes[0].coordinates[1]
    
    Dx=float(4*X-(X1+X2+X3+X4))
    Dy=float(4*Y-(Y1+Y2+Y3+Y4))
            
    Ax=float(-X1+X2+X3-X4)
    Ay=float(-Y1-Y2+Y3+Y4)
    
    tsi=float(Dx/Ax);
    eta=float(Dy/Ay);
    assert abs(tsi)<=1 and abs(eta)<=1
    
    N1_d2y=float(-12.0/(Ay**2)*(1-tsi)*(-1)*eta)
    N2_d2y=float(-12.0/(Ay**2)*(1+tsi)*(-1)*eta)
    N3_d2y=float(-12.0/(Ay**2)*(1+tsi)*(1)*eta)
    N4_d2y=float(-12.0/(Ay**2)*(1-tsi)*(1)*eta)
    
    NX1_d2y=float(1/Ay*(1-tsi)*(-1+3*eta))
    NX2_d2y=float(1/Ay*(1+tsi)*(-1+3*eta))
    NX3_d2y=float(1/Ay*(1+tsi)*(1+3*eta))
    NX4_d2y=float(1/Ay*(1-tsi)*(1+3*eta))
    
    NY1_d2y=float(0)
    NY2_d2y=float(0)
    NY3_d2y=float(0)
    NY4_d2y=float(0)
    
    return N1_d2y,N2_d2y,N3_d2y,N4_d2y,NX1_d2y,NX2_d2y,NX3_d2y,NX4_d2y,NY1_d2y,NY2_d2y,NY3_d2y,NY4_d2y


def d2w_dxdy_shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4):    
    X=P.nodes[0].coordinates[0]
    Y=P.nodes[0].coordinates[1]
    
    Dx=float(4*X-(X1+X2+X3+X4))
    Dy=float(4*Y-(Y1+Y2+Y3+Y4))
            
    Ax=float(-X1+X2+X3-X4)
    Ay=float(-Y1-Y2+Y3+Y4)
    
    tsi=float(Dx/Ax);
    eta=float(Dy/Ay);
    assert abs(tsi)<=1 and abs(eta)<=1
    
    N1_dxdy=float(2/(Ay*Ax)*(-1)*(-1)*(4-3*tsi**2-3*eta**2))
    N2_dxdy=float(2/(Ay*Ax)*(1)*(-1)*(4-3*tsi**2-3*eta**2))
    N3_dxdy=float(2/(Ay*Ax)*(1)*(1)*(4-3*tsi**2-3*eta**2))
    N4_dxdy=float(2/(Ay*Ax)*(-1)*(1)*(4-3*tsi**2-3*eta**2))


    NX1_dxdy=float(-0.5/Ax*(-1)*(-1)*(-1-2*eta-3*(-1)*eta**2))
    NX2_dxdy=float(-0.5/Ax*(1)*(-1)*(-1-2*eta-3*(-1)*eta**2))
    NX3_dxdy=float(-0.5/Ax*(1)*(1)*(1-2*eta-3*(1)*eta**2))
    NX4_dxdy=float(-0.5/Ax*(-1)*(1)*(1-2*eta-3*(1)*eta**2))

    NY1_dxdy=float(0.5/Ay*(-1)*(-1)*(-1-2*tsi-3*(-1)*tsi**2))
    NY2_dxdy=float(0.5/Ay*(1)*(-1)*(1-2*tsi-3*(1)*tsi**2))
    NY3_dxdy=float(0.5/Ay*(1)*(1)*(1-2*tsi-3*(1)*tsi**2))
    NY4_dxdy=float(0.5/Ay*(-1)*(1)*(-1-2*tsi-3*(-1)*tsi**2))


    return N1_dxdy,N2_dxdy,N3_dxdy,N4_dxdy,NX1_dxdy,NX2_dxdy,NX3_dxdy,NX4_dxdy,NY1_dxdy,NY2_dxdy,NY3_dxdy,NY4_dxdy





nodal_connectvity=np.array([[9, 8, 1, 2],
[10, 9, 2, 3],
[11, 10, 3, 4],
[12, 11, 4, 5],
[13, 12, 5, 6],
[14, 13, 6, 7],
[16, 15, 8, 9],
[17, 16, 9, 10],
[18, 17, 10, 11],
[19, 18, 11, 12],
[20, 19, 12, 13],
[21, 20, 13, 14],
[23, 22, 15, 16],
[24, 23, 16, 17],
[25, 24, 17, 18],
[26, 25, 18, 19],
[27, 26, 19, 20],
[28, 27, 20, 21],
[30, 29, 22, 23],
[31, 30, 23, 24],
[32, 31, 24, 25],
[33, 32, 25, 26],
[34, 33, 26, 27],
[35, 34, 27, 28],
[37, 36, 29, 30],
[38, 37, 30, 31],
[39, 38, 31, 32],
[40, 39, 32, 33],
[41, 40, 33, 34],
[42, 41, 34, 35],
[44, 43, 36, 37],
[45, 44, 37, 38],
[46, 45, 38, 39],
[47, 46, 39, 40],
[48, 47, 40, 41],
[49, 48, 41, 42]])





####################################################
###Define mesh set
####################################################



for ii in range(0,Num_element): #loop to go around element
    ele=ii+1;#element number
    nodes=nodal_connectvity[ii]
    X1=mdb.models[modelName].rootAssembly.features['RP-'+str(nodes[0])].xValue
    X2=mdb.models[modelName].rootAssembly.features['RP-'+str(nodes[1])].xValue
    X3=mdb.models[modelName].rootAssembly.features['RP-'+str(nodes[2])].xValue
    X4=mdb.models[modelName].rootAssembly.features['RP-'+str(nodes[3])].xValue
    Y1=mdb.models[modelName].rootAssembly.features['RP-'+str(nodes[0])].yValue
    Y2=mdb.models[modelName].rootAssembly.features['RP-'+str(nodes[1])].yValue
    Y3=mdb.models[modelName].rootAssembly.features['RP-'+str(nodes[2])].yValue
    Y4=mdb.models[modelName].rootAssembly.features['RP-'+str(nodes[3])].yValue
    
    
    for jj in range(1,Num_Gauss):# loop to go around the gauss point RVEs
        instanceName = 'Part-1-ele-'+str(ele)+'-GP-'+str(jj)
        Dimension=Get_Cube_dimension(modelName,instanceName)
        Xmin = Dimension[0]
        Xmax = Dimension[1]
        Ymin = Dimension[2]
        Ymax = Dimension[3]
        Zmin = Dimension[4]
        Zmax = Dimension[5]
        eps1 = abs(Zmax - Zmin) * 0.0001
        eps2 = abs(Zmax - Zmin) * 0.01
        tolerance = 0.1
        BX = Xmax - Xmin
        BY = Ymax - Ymin
        BZ = Zmax - Zmin
        node = mdb.models[modelName].rootAssembly.instances[instanceName].nodes
        node_E1 = node[1:1]
        node_E2 = node[1:1]
        node_E3 = node[1:1]
        node_E4 = node[1:1]
        node_E5 = node[1:1]
        node_E6 = node[1:1]
        node_E7 = node[1:1]
        node_E8 = node[1:1]
        node_E9 = node[1:1]
        node_E10 = node[1:1]
        node_E11 = node[1:1]
        node_E12 = node[1:1]
        node_FXP = node[1:1]
        node_FXN = node[1:1]
        node_FYP = node[1:1]
        node_FYN = node[1:1]
        node_FZP = node[1:1]
        node_FZN = node[1:1]
        for i in range(len(node)):
            x = node[i].coordinates[0]
            y = node[i].coordinates[1]
            z = node[i].coordinates[2]
            if abs(x - Xmin) < eps1 and abs(y - Ymin) < eps1 and abs(z - Zmin) < eps1:
                node_v1 = node[i:i + 1]
                continue
            if abs(x - Xmax) < eps1 and abs(y - Ymin) < eps1 and abs(z - Zmin) < eps1:
                node_v2 = node[i:i + 1]
                continue
            if abs(x - Xmax) < eps1 and abs(y - Ymax) < eps1 and abs(z - Zmin) < eps1:
                node_v3 = node[i:i + 1]
                continue
            if abs(x - Xmin) < eps1 and abs(y - Ymax) < eps1 and abs(z - Zmin) < eps1:
                node_v4 = node[i:i + 1]
                continue
            if abs(x - Xmin) < eps1 and abs(y - Ymin) < eps1 and abs(z - Zmax) < eps1:
                node_v5 = node[i:i + 1]
                continue
            if abs(x - Xmax) < eps1 and abs(y - Ymin) < eps1 and abs(z - Zmax) < eps1:
                node_v6 = node[i:i + 1]
                continue
            if abs(x - Xmax) < eps1 and abs(y - Ymax) < eps1 and abs(z - Zmax) < eps1:
                node_v7 = node[i:i + 1]
                continue
            if abs(x - Xmin) < eps1 and abs(y - Ymax) < eps1 and abs(z - Zmax) < eps1:
                node_v8 = node[i:i + 1]
                continue
            if abs(x - Xmin) < eps1 and abs(y - Ymin) < eps1 and abs(z - Zmin) > eps2 and abs(z - Zmax) > eps2:
                node_E1 = node_E1 + node[i:i + 1]
                continue
            if abs(x - Xmax) < eps1 and abs(y - Ymin) < eps1 and abs(z - Zmin) > eps2 and abs(z - Zmax) > eps2:
                node_E2 = node_E2 + node[i:i + 1]
                continue
            if abs(x - Xmax) < eps1 and abs(y - Ymax) < eps1 and abs(z - Zmin) > eps2 and abs(z - Zmax) > eps2:
                node_E3 = node_E3 + node[i:i + 1]
                continue
            if abs(x - Xmin) < eps1 and abs(y - Ymax) < eps1 and abs(z - Zmin) > eps2 and abs(z - Zmax) > eps2:
                node_E4 = node_E4 + node[i:i + 1]
                continue
            if abs(x - Xmin) < eps1 and abs(z - Zmin) < eps1 and abs(y - Ymin) > eps2 and abs(y - Ymax) > eps2:
                node_E5 = node_E5 + node[i:i + 1]
                continue
            if abs(x - Xmax) < eps1 and abs(z - Zmin) < eps1 and abs(y - Ymin) > eps2 and abs(y - Ymax) > eps2:
                node_E6 = node_E6 + node[i:i + 1]
                continue
            if abs(x - Xmax) < eps1 and abs(z - Zmax) < eps1 and abs(y - Ymin) > eps2 and abs(y - Ymax) > eps2:
                node_E7 = node_E7 + node[i:i + 1]
                continue
            if abs(x - Xmin) < eps1 and abs(z - Zmax) < eps1 and abs(y - Ymin) > eps2 and abs(y - Ymax) > eps2:
                node_E8 = node_E8 + node[i:i + 1]
                continue
            if abs(y - Ymin) < eps1 and abs(z - Zmin) < eps1 and abs(x - Xmin) > eps2 and abs(x - Xmax) > eps2:
                node_E9 = node_E9 + node[i:i + 1]
                continue
            if abs(y - Ymax) < eps1 and abs(z - Zmin) < eps1 and abs(x - Xmin) > eps2 and abs(x - Xmax) > eps2:
                node_E10 = node_E10 + node[i:i + 1]
                continue
            if abs(y - Ymax) < eps1 and abs(z - Zmax) < eps1 and abs(x - Xmin) > eps2 and abs(x - Xmax) > eps2:
                node_E11 = node_E11 + node[i:i + 1]
                continue
            if abs(y - Ymin) < eps1 and abs(z - Zmax) < eps1 and abs(x - Xmin) > eps2 and abs(x - Xmax) > eps2:
                node_E12 = node_E12 + node[i:i + 1]
                continue
            if abs(x - Xmax) < eps1 and abs(y - Ymin) > eps2 and abs(y - Ymax) > eps2 and abs(z - Zmin) > eps2 and abs(z - Zmax) > eps2:
                node_FXP = node_FXP + node[i:i + 1]
                continue
            if abs(x - Xmin) < eps1 and abs(y - Ymin) > eps2 and abs(y - Ymax) > eps2 and abs(z - Zmin) > eps2 and abs(z - Zmax) > eps2:
                node_FXN = node_FXN + node[i:i + 1]
                continue
            if abs(y - Ymax) < eps1 and abs(x - Xmin) > eps2 and abs(x - Xmax) > eps2 and abs(z - Zmin) > eps2 and abs(z - Zmax) > eps2:
                node_FYP = node_FYP + node[i:i + 1]
                continue
            if abs(y - Ymin) < eps1 and abs(x - Xmin) > eps2 and abs(x - Xmax) > eps2 and abs(z - Zmin) > eps2 and abs(z - Zmax) > eps2:
                node_FYN = node_FYN + node[i:i + 1]
                continue
            if abs(z - Zmax) < eps1 and abs(y - Ymin) > eps2 and abs(y - Ymax) > eps2 and abs(x - Xmin) > eps2 and abs(x - Xmax) > eps2:
                node_FZP = node_FZP + node[i:i + 1]
                continue
            if abs(z - Zmin) < eps1 and abs(y - Ymin) > eps2 and abs(y - Ymax) > eps2 and abs(x - Xmin) > eps2 and abs(x - Xmax) > eps2:
                node_FZN = node_FZN + node[i:i + 1]
                continue
    
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-1', nodes = node_v1)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-2', nodes = node_v2)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-3', nodes = node_v3)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-4', nodes = node_v4)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-5', nodes = node_v5)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-6', nodes = node_v6)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-7', nodes = node_v7)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-8', nodes = node_v8)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-edge1', nodes = node_E1)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-edge2', nodes = node_E2)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-edge3', nodes = node_E3)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-edge4', nodes = node_E4)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-edge5', nodes = node_E5)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-edge6', nodes = node_E6)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-edge7', nodes = node_E7)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-edge8', nodes = node_E8)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-edge9', nodes = node_E9)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-edge10', nodes = node_E10)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-edge11', nodes = node_E11)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-edge12', nodes = node_E12)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Face-XP', nodes = node_FXP)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Face-XN', nodes = node_FXN)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Face-YP', nodes = node_FYP)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Face-YN', nodes = node_FYN)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Face-ZP', nodes = node_FZP)
        mdb.models[modelName].rootAssembly.Set(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Face-ZN', nodes = node_FZN)
        nodesI = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge1'].nodes
        nodesII = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge2'].nodes
        nodesIII = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge3'].nodes
        nodesIV = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge4'].nodes
        nodesV = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge5'].nodes
        nodesVI = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge6'].nodes
        nodesVII = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge7'].nodes
        nodesVIII = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge8'].nodes
        nodesIX = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge9'].nodes
        nodesX = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge10'].nodes
        nodesXI = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge11'].nodes
        nodesXII = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge12'].nodes
        maxIndex = len(nodesI)
        countII = 0
        countIII = 0
        countIV = 0
        for i in range(maxIndex):
            x_1 = nodesI[i].coordinates[0]
            y_1 = nodesI[i].coordinates[1]
            z_1 = nodesI[i].coordinates[2]
            nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge1-')
            mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesI[i:i + 1])
            maxIndexTarget = len(nodesII)
            for j in range(maxIndexTarget):
                x_2 = nodesII[j].coordinates[0]
                y_2 = nodesII[j].coordinates[1]
                z_2 = nodesII[j].coordinates[2]
                distance = abs((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2 - BX ** 2)
                if distance <= tolerance:
                    countII = countII + 1
                    nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge2-')
                    mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesII[j:j + 1])
                    indexNumberToRemove = j
                    break
                    continue
        
            nodesII = nodesII[:indexNumberToRemove] + nodesII[indexNumberToRemove + 1:]
            maxIndexTarget = len(nodesIII)
            for j in range(maxIndexTarget):
                x_2 = nodesIII[j].coordinates[0]
                y_2 = nodesIII[j].coordinates[1]
                z_2 = nodesIII[j].coordinates[2]
                distance = abs((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2 - BX ** 2 - BY ** 2)
                if distance <= tolerance:
                    countIII = countIII + 1
                    nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge3-')
                    mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesIII[j:j + 1])
                    indexNumberToRemove = j
                    break
                    continue
        
            nodesIII = nodesIII[:indexNumberToRemove] + nodesIII[indexNumberToRemove + 1:]
            maxIndexTarget = len(nodesIV)
            for j in range(maxIndexTarget):
                x_2 = nodesIV[j].coordinates[0]
                y_2 = nodesIV[j].coordinates[1]
                z_2 = nodesIV[j].coordinates[2]
                distance = abs((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2 - BY ** 2)
                if distance <= tolerance:
                    countIV = countIV + 1
                    nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge4-')
                    mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesIV[j:j + 1])
                    indexNumberToRemove = j
                    break
                    continue
        
            nodesIV = nodesIV[:indexNumberToRemove] + nodesIV[indexNumberToRemove + 1:]
    
        if countII == maxIndex:
            print 'Node matching edge I-II successful'
        else:
            print 'Node matching edge I-II has failed'
        if countIII == maxIndex:
            print 'Node matching edge I-III successful'
        else:
            print 'Node matching edge I-III has failed'
        if countIV == maxIndex:
            print 'Node matching edge I-IV successful'
        else:
            print 'Node matching edge I-IV has failed'
        maxIndex = len(nodesV)
        countVI = 0
        countVII = 0
        countVIII = 0
        for i in range(maxIndex):
            x_1 = nodesV[i].coordinates[0]
            y_1 = nodesV[i].coordinates[1]
            z_1 = nodesV[i].coordinates[2]
            nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge5-')
            mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesV[i:i + 1])
            maxIndexTarget = len(nodesVI)
            for j in range(maxIndexTarget):
                x_2 = nodesVI[j].coordinates[0]
                y_2 = nodesVI[j].coordinates[1]
                z_2 = nodesVI[j].coordinates[2]
                distance = abs((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2 - BX ** 2)
                if distance <= tolerance:
                    countVI = countVI + 1
                    nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge6-')
                    mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesVI[j:j + 1])
                    indexNumberToRemove = j
                    break
                    continue
        
            nodesVI = nodesVI[:indexNumberToRemove] + nodesVI[indexNumberToRemove + 1:]
            maxIndexTarget = len(nodesVII)
            for j in range(maxIndexTarget):
                x_2 = nodesVII[j].coordinates[0]
                y_2 = nodesVII[j].coordinates[1]
                z_2 = nodesVII[j].coordinates[2]
                distance = abs((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2 - BX ** 2 - BZ ** 2)
                if distance <= tolerance:
                    countVII = countVII + 1
                    nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge7-')
                    mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesVII[j:j + 1])
                    indexNumberToRemove = j
                    break
                    continue
        
            nodesVII = nodesVII[:indexNumberToRemove] + nodesVII[indexNumberToRemove + 1:]
            maxIndexTarget = len(nodesVIII)
            for j in range(maxIndexTarget):
                x_2 = nodesVIII[j].coordinates[0]
                y_2 = nodesVIII[j].coordinates[1]
                z_2 = nodesVIII[j].coordinates[2]
                distance = abs((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2 - BZ ** 2)
                if distance <= tolerance:
                    countVIII = countVIII + 1
                    nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge8-')
                    mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesVIII[j:j + 1])
                    indexNumberToRemove = j
                    break
                    continue
        
            nodesVIII = nodesVIII[:indexNumberToRemove] + nodesVIII[indexNumberToRemove + 1:]
    
        if countVI == maxIndex:
            print 'Node matching edge V-VI successful'
        else:
            print 'Node matching edge V-VI has failed'
        if countVII == maxIndex:
            print 'Node matching edge V-VII successful'
        else:
            print 'Node matching edge V-VII has failed'
        if countVIII == maxIndex:
            print 'Node matching edge V-VIII successful'
        else:
            print 'Node matching edge V-VIII has failed'
        maxIndex = len(nodesIX)
        countX = 0
        countXI = 0
        countXII = 0
        for i in range(maxIndex):
            x_1 = nodesIX[i].coordinates[0]
            y_1 = nodesIX[i].coordinates[1]
            z_1 = nodesIX[i].coordinates[2]
            nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge9-')
            mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesIX[i:i + 1])
            maxIndexTarget = len(nodesX)
            for j in range(maxIndexTarget):
                x_2 = nodesX[j].coordinates[0]
                y_2 = nodesX[j].coordinates[1]
                z_2 = nodesX[j].coordinates[2]
                distance = abs((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2 - BY ** 2)
                if distance <= tolerance:
                    countX = countX + 1
                    nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge10-')
                    mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesX[j:j + 1])
                    indexNumberToRemove = j
                    break
                    continue
        
            nodesX = nodesX[:indexNumberToRemove] + nodesX[indexNumberToRemove + 1:]
            maxIndexTarget = len(nodesXI)
            for j in range(maxIndexTarget):
                x_2 = nodesXI[j].coordinates[0]
                y_2 = nodesXI[j].coordinates[1]
                z_2 = nodesXI[j].coordinates[2]
                distance = abs((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2 - BY ** 2 - BZ ** 2)
                if distance <= tolerance:
                    countXI = countXI + 1
                    nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge11-')
                    mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesXI[j:j + 1])
                    indexNumberToRemove = j
                    break
                    continue
        
            nodesXI = nodesXI[:indexNumberToRemove] + nodesXI[indexNumberToRemove + 1:]
            maxIndexTarget = len(nodesXII)
            for j in range(maxIndexTarget):
                x_2 = nodesXII[j].coordinates[0]
                y_2 = nodesXII[j].coordinates[1]
                z_2 = nodesXII[j].coordinates[2]
                distance = abs((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2 - BZ ** 2)
                if distance <= tolerance:
                    countXII = countXII + 1
                    nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge12-')
                    mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesXII[j:j + 1])
                    indexNumberToRemove = j
                    break
                    continue
        
            nodesXII = nodesXII[:indexNumberToRemove] + nodesXII[indexNumberToRemove + 1:]
    
        if countX == maxIndex:
            print 'Node matching edge IX-X successful'
        else:
            print 'Node matching edge IX-X has failed'
        if countXI == maxIndex:
            print 'Node matching edge IX-XI successful'
        else:
            print 'Node matching edge IX-XI has failed'
        if countXII == maxIndex:
            print 'Node matching edge IX-XII successful'
        else:
            print 'Node matching edge IX-XII has failed'
        nodesFaceXN = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Face-XN'].nodes
        nodesFaceXP = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Face-XP'].nodes
        maxIndex = len(nodesFaceXN)
        print 'Number of nodes Face XN:', maxIndex
        count = 0
        for i in range(maxIndex):
            x_1 = nodesFaceXN[i].coordinates[0]
            y_1 = nodesFaceXN[i].coordinates[1]
            z_1 = nodesFaceXN[i].coordinates[2]
            nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceXN-')
            mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesFaceXN[i:i + 1])
            maxIndexTarget = len(nodesFaceXP)
            for j in range(maxIndexTarget):
                x_2 = nodesFaceXP[j].coordinates[0]
                y_2 = nodesFaceXP[j].coordinates[1]
                z_2 = nodesFaceXP[j].coordinates[2]
                distance = abs((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2 - BX ** 2)
                if distance <= tolerance:
                    count = count + 1
                    nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceXP-')
                    mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesFaceXP[j:j + 1])
                    indexNumberToRemove = j
                    break
                    continue
        
            nodesFaceXP = nodesFaceXP[:indexNumberToRemove] + nodesFaceXP[indexNumberToRemove + 1:]
    
        if count == maxIndex:
            print 'Node matching XN-XP successful'
        else:
            print 'Node matching XN-XP has failed: ', maxIndex - count
        nodesFaceYN = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Face-YN'].nodes
        nodesFaceYP = mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Face-YP'].nodes
        maxIndex = len(nodesFaceYN)
        print 'Number of nodes Face YN:', maxIndex
        count = 0
        for i in range(maxIndex):
            x_1 = nodesFaceYN[i].coordinates[0]
            y_1 = nodesFaceYN[i].coordinates[1]
            z_1 = nodesFaceYN[i].coordinates[2]
            nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceYN-')
            mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesFaceYN[i:i + 1])
            maxIndexTarget = len(nodesFaceYP)
            for j in range(maxIndexTarget):
                x_2 = nodesFaceYP[j].coordinates[0]
                y_2 = nodesFaceYP[j].coordinates[1]
                z_2 = nodesFaceYP[j].coordinates[2]
                distance = abs((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2 - BY ** 2)
                if distance <= tolerance:
                    count = count + 1
                    nodeLabel = LabelName(i, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceYP-')
                    mdb.models[modelName].rootAssembly.Set(name = nodeLabel, nodes = nodesFaceYP[j:j + 1])
                    indexNumberToRemove = j
                    break
                    continue
        
            nodesFaceYP = nodesFaceYP[:indexNumberToRemove] + nodesFaceYP[indexNumberToRemove + 1:]
    
        if count == maxIndex:
            print 'Node matching  YN-YP successful'
        else:
            print 'Node matching YN-YP has failed: ', maxIndex - count



        print 'Creat Mesh sets Complete La La La'


####################################################
###Define PBC
####################################################
        P=mdb.models[modelName].rootAssembly.sets['Part-1-ele-'+str(ele)+'-GP-'+str(jj)+'.M00']
        N1,N2,N3,N4,NX1,NX2,NX3,NX4,NY1,NY2,NY3,NY4=shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
        mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Center-node-M00-on-w', 
        terms = ((1, 'Part-1-ele-'+str(ele)+'-GP-'+str(jj)+'.M00', 3), 
        (-N1, 'N'+str(nodes[0])+'-RP', 3), (-NX1, 'N'+str(nodes[0])+'-RP', 4),(-NY1, 'N'+str(nodes[0])+'-RP', 5),
        (-N2, 'N'+str(nodes[1])+'-RP', 3), (-NX2, 'N'+str(nodes[1])+'-RP', 4),(-NY2, 'N'+str(nodes[1])+'-RP', 5),
        (-N3, 'N'+str(nodes[2])+'-RP', 3), (-NX3, 'N'+str(nodes[2])+'-RP', 4),(-NY3, 'N'+str(nodes[2])+'-RP', 5),
        (-N4, 'N'+str(nodes[3])+'-RP', 3), (-NX4, 'N'+str(nodes[3])+'-RP', 4),(-NY4, 'N'+str(nodes[3])+'-RP', 5)))


        N1_H,N2_H,N3_H,N4_H,N1_H_dx,N2_H_dx,N3_H_dx,N4_H_dx,N1_H_dy,N2_H_dy,N3_H_dy,N4_H_dy=Horizontal_shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
        mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Center-node-M00-on-u', 
        terms = ((1, 'Part-1-ele-'+str(ele)+'-GP-'+str(jj)+'.M00', 1), 
        (-N1_H, 'N'+str(nodes[0])+'-RP', 1),
        (-N2_H, 'N'+str(nodes[1])+'-RP', 1),
        (-N3_H, 'N'+str(nodes[2])+'-RP', 1),
        (-N4_H, 'N'+str(nodes[3])+'-RP', 1)))

        mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Center-node-M00-on-v', 
        terms = ((1, 'Part-1-ele-'+str(ele)+'-GP-'+str(jj)+'.M00', 2),
        (-N1_H, 'N'+str(nodes[0])+'-RP', 2),
        (-N2_H, 'N'+str(nodes[1])+'-RP', 2),
        (-N3_H, 'N'+str(nodes[2])+'-RP', 2),
        (-N4_H, 'N'+str(nodes[3])+'-RP', 2)))






        P=mdb.models[modelName].rootAssembly.sets['Part-1-ele-'+str(ele)+'-GP-'+str(jj)+'.M0']
        N1,N2,N3,N4,NX1,NX2,NX3,NX4,NY1,NY2,NY3,NY4=shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
        N1_d2x,N2_d2x,N3_d2x,N4_d2x,NX1_d2x,NX2_d2x,NX3_d2x,NX4_d2x,NY1_d2x,NY2_d2x,NY3_d2x,NY4_d2x=d2w_dx2_shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
        N1_d2y,N2_d2y,N3_d2y,N4_d2y,NX1_d2y,NX2_d2y,NX3_d2y,NX4_d2y,NY1_d2y,NY2_d2y,NY3_d2y,NY4_d2y=d2w_dy2_shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
        N1_H,N2_H,N3_H,N4_H,N1_H_dx,N2_H_dx,N3_H_dx,N4_H_dx,N1_H_dy,N2_H_dy,N3_H_dy,N4_H_dy=Horizontal_shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
        N1_dxdy,N2_dxdy,N3_dxdy,N4_dxdy,NX1_dxdy,NX2_dxdy,NX3_dxdy,NX4_dxdy,NY1_dxdy,NY2_dxdy,NY3_dxdy,NY4_dxdy=d2w_dxdy_shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4)


        N1_dx,N2_dx,N3_dx,N4_dx,NX1_dx,NX2_dx,NX3_dx,NX4_dx,NY1_dx,NY2_dx,NY3_dx,NY4_dx=minus_dw_dx_shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
        N1_dy,N2_dy,N3_dy,N4_dy,NX1_dy,NX2_dy,NX3_dy,NX4_dy,NY1_dy,NY2_dy,NY3_dy,NY4_dy=dw_dy_shape_fn(P,X1,X2,X3,X4,Y1,Y2,Y3,Y4)

        # print N1_dx,N2_dx,N3_dx,N4_dx,NX1_dx,NX2_dx,NX3_dx,NX4_dx,NY1_dx,NY2_dx,NY3_dx,NY4_dx
        # print N1_dy,N2_dy,N3_dy,N4_dy,NX1_dy,NX2_dy,NX3_dy,NX4_dy,NY1_dy,NY2_dy,NY3_dy,NY4_dy

        #remove Rigid body displacement 
        #M2-M1 vertical 
        P1=mdb.models[modelName].rootAssembly.sets['Part-1-ele-'+str(ele)+'-GP-'+str(jj)+'.M1']
        P2=mdb.models[modelName].rootAssembly.sets['Part-1-ele-'+str(ele)+'-GP-'+str(jj)+'.M2']
        Xcoord_P1=P1.nodes[0].coordinates[0]
        Ycoord_P1=P1.nodes[0].coordinates[1]
        Zcoord_P1=P1.nodes[0].coordinates[2]
        Xcoord_P2=P2.nodes[0].coordinates[0]
        Ycoord_P2=P2.nodes[0].coordinates[1]
        Zcoord_P2=P2.nodes[0].coordinates[2]
        Length_X=float(abs(Xcoord_P2-Xcoord_P1))
        Length_Y=float(abs(Ycoord_P2-Ycoord_P1))
        # print Length_X, Length_Y
        mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Center-node-M2-M1-on-w', 
        terms = ((1, 'Part-1-ele-'+str(ele)+'-GP-'+str(jj)+'.M2', 3), (-1, 'Part-1-ele-'+str(ele)+'-GP-'+str(jj)+'.M1', 3),
        (-Length_X*N1_dx, 'N'+str(nodes[0])+'-RP', 3), (-Length_X*NX1_dx, 'N'+str(nodes[0])+'-RP', 4),(-Length_X*NY1_dx, 'N'+str(nodes[0])+'-RP', 5),
        (-Length_X*N2_dx, 'N'+str(nodes[1])+'-RP', 3), (-Length_X*NX2_dx, 'N'+str(nodes[1])+'-RP', 4),(-Length_X*NY2_dx, 'N'+str(nodes[1])+'-RP', 5),
        (-Length_X*N3_dx, 'N'+str(nodes[2])+'-RP', 3), (-Length_X*NX3_dx, 'N'+str(nodes[2])+'-RP', 4),(-Length_X*NY3_dx, 'N'+str(nodes[2])+'-RP', 5),
        (-Length_X*N4_dx, 'N'+str(nodes[3])+'-RP', 3), (-Length_X*NX4_dx, 'N'+str(nodes[3])+'-RP', 4),(-Length_X*NY4_dx, 'N'+str(nodes[3])+'-RP', 5)))

        #M4-M1 vertical 
        P1=mdb.models[modelName].rootAssembly.sets['Part-1-ele-'+str(ele)+'-GP-'+str(jj)+'.M1']
        P2=mdb.models[modelName].rootAssembly.sets['Part-1-ele-'+str(ele)+'-GP-'+str(jj)+'.M4']
        Xcoord_P1=P1.nodes[0].coordinates[0]
        Ycoord_P1=P1.nodes[0].coordinates[1]
        Zcoord_P1=P1.nodes[0].coordinates[2]
        Xcoord_P2=P2.nodes[0].coordinates[0]
        Ycoord_P2=P2.nodes[0].coordinates[1]
        Zcoord_P2=P2.nodes[0].coordinates[2]
        Length_X=-float(abs(Xcoord_P2-Xcoord_P1))
        Length_Y=float(abs(Ycoord_P2-Ycoord_P1))
        mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Center-node-M4-M1-on-w', 
        terms = ((1, 'Part-1-ele-'+str(ele)+'-GP-'+str(jj)+'.M4', 3), (-1, 'Part-1-ele-'+str(ele)+'-GP-'+str(jj)+'.M1', 3),
        (-Length_Y*N1_dy, 'N'+str(nodes[0])+'-RP', 3), (-Length_Y*NX1_dy, 'N'+str(nodes[0])+'-RP', 4),(-Length_Y*NY1_dy, 'N'+str(nodes[0])+'-RP', 5),
        (-Length_Y*N2_dy, 'N'+str(nodes[1])+'-RP', 3), (-Length_Y*NX2_dy, 'N'+str(nodes[1])+'-RP', 4),(-Length_Y*NY2_dy, 'N'+str(nodes[1])+'-RP', 5),
        (-Length_Y*N3_dy, 'N'+str(nodes[2])+'-RP', 3), (-Length_Y*NX3_dy, 'N'+str(nodes[2])+'-RP', 4),(-Length_Y*NY3_dy, 'N'+str(nodes[2])+'-RP', 5),
        (-Length_Y*N4_dy, 'N'+str(nodes[3])+'-RP', 3), (-Length_Y*NX4_dy, 'N'+str(nodes[3])+'-RP', 4),(-Length_Y*NY4_dy, 'N'+str(nodes[3])+'-RP', 5)))









        A=1
        B=1
        C=1

        #master node 2-1,3-1,4-1
        for k in range(2,5):  
            P1=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-'+str(1)]
            P2=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-'+str(k)]
            Xcoord_P1=P1.nodes[0].coordinates[0]
            Ycoord_P1=P1.nodes[0].coordinates[1]
            Zcoord_P1=P1.nodes[0].coordinates[2]
            Xcoord_P2=P2.nodes[0].coordinates[0]
            Ycoord_P2=P2.nodes[0].coordinates[1]
            Zcoord_P2=P2.nodes[0].coordinates[2]
            Length_X=float(abs(Xcoord_P2-Xcoord_P1))
            Length_Y=float(abs(Ycoord_P2-Ycoord_P1))
            Zcoord=(Zcoord_P1+Zcoord_P2)/2*Scale
            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-node-'+str(k)+'-1'+'-on-u', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-'+str(k), 1), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-'+str(1), 1),
            # plane
            (Length_X*-N1_H_dx+A*Length_Y*-N1_H_dy, 'N'+str(nodes[0])+'-RP', 1),
            (Length_X*-N2_H_dx+A*Length_Y*-N2_H_dy, 'N'+str(nodes[1])+'-RP', 1),
            (Length_X*-N3_H_dx+A*Length_Y*-N3_H_dy, 'N'+str(nodes[2])+'-RP', 1),
            (Length_X*-N4_H_dx+A*Length_Y*-N4_H_dy, 'N'+str(nodes[3])+'-RP', 1),
            # bending
            (Length_X*Zcoord*N1_d2x+C*Length_Y*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_X*Zcoord*NX1_d2x+C*Length_Y*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_X*Zcoord*NY1_d2x+C*Length_Y*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_X*Zcoord*N2_d2x+C*Length_Y*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_X*Zcoord*NX2_d2x+C*Length_Y*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_X*Zcoord*NY2_d2x+C*Length_Y*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_X*Zcoord*N3_d2x+C*Length_Y*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_X*Zcoord*NX3_d2x+C*Length_Y*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_X*Zcoord*NY3_d2x+C*Length_Y*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_X*Zcoord*N4_d2x+C*Length_Y*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_X*Zcoord*NX4_d2x+C*Length_Y*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_X*Zcoord*NY4_d2x+C*Length_Y*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))







            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-node-'+str(k)+'-1'+'-on-v', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-'+str(k), 2), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-'+str(1), 2),
            # plane
            (Length_Y*-N1_H_dy+A*Length_X*-N1_H_dx, 'N'+str(nodes[0])+'-RP', 2),
            (Length_Y*-N2_H_dy+A*Length_X*-N2_H_dx, 'N'+str(nodes[1])+'-RP', 2),
            (Length_Y*-N3_H_dy+A*Length_X*-N3_H_dx, 'N'+str(nodes[2])+'-RP', 2),
            (Length_Y*-N4_H_dy+A*Length_X*-N4_H_dx, 'N'+str(nodes[3])+'-RP', 2),
            # bending
            (Length_Y*Zcoord*N1_d2y+B*Length_X*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_Y*Zcoord*NX1_d2y+B*Length_X*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_Y*Zcoord*NY1_d2y+B*Length_X*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_Y*Zcoord*N2_d2y+B*Length_X*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_Y*Zcoord*NX2_d2y+B*Length_X*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_Y*Zcoord*NY2_d2y+B*Length_X*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_Y*Zcoord*N3_d2y+B*Length_X*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_Y*Zcoord*NX3_d2y+B*Length_X*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_Y*Zcoord*NY3_d2y+B*Length_X*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_Y*Zcoord*N4_d2y+B*Length_X*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_Y*Zcoord*NX4_d2y+B*Length_X*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_Y*Zcoord*NY4_d2y+B*Length_X*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))












        #master node 6-5,7-5,8-5
        for k in range(6,9):  
            P1=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-'+str(5)]
            P2=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-'+str(k)]
            Xcoord_P1=P1.nodes[0].coordinates[0]
            Ycoord_P1=P1.nodes[0].coordinates[1]
            Zcoord_P1=P1.nodes[0].coordinates[2]
            Xcoord_P2=P2.nodes[0].coordinates[0]
            Ycoord_P2=P2.nodes[0].coordinates[1]
            Zcoord_P2=P2.nodes[0].coordinates[2]
            Length_X=float(abs(Xcoord_P2-Xcoord_P1))
            Length_Y=float(abs(Ycoord_P2-Ycoord_P1))
            Zcoord=(Zcoord_P1+Zcoord_P2)/2*Scale
            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-node-'+str(k)+'-5'+'-on-u', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-'+str(k), 1), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-'+str(5), 1),
            # plane
            (Length_X*-N1_H_dx+A*Length_Y*-N1_H_dy, 'N'+str(nodes[0])+'-RP', 1),
            (Length_X*-N2_H_dx+A*Length_Y*-N2_H_dy, 'N'+str(nodes[1])+'-RP', 1),
            (Length_X*-N3_H_dx+A*Length_Y*-N3_H_dy, 'N'+str(nodes[2])+'-RP', 1),
            (Length_X*-N4_H_dx+A*Length_Y*-N4_H_dy, 'N'+str(nodes[3])+'-RP', 1),
            # bending
            (Length_X*Zcoord*N1_d2x+C*Length_Y*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_X*Zcoord*NX1_d2x+C*Length_Y*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_X*Zcoord*NY1_d2x+C*Length_Y*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_X*Zcoord*N2_d2x+C*Length_Y*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_X*Zcoord*NX2_d2x+C*Length_Y*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_X*Zcoord*NY2_d2x+C*Length_Y*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_X*Zcoord*N3_d2x+C*Length_Y*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_X*Zcoord*NX3_d2x+C*Length_Y*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_X*Zcoord*NY3_d2x+C*Length_Y*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_X*Zcoord*N4_d2x+C*Length_Y*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_X*Zcoord*NX4_d2x+C*Length_Y*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_X*Zcoord*NY4_d2x+C*Length_Y*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))


            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-node-'+str(k)+'-5'+'-on-v', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-'+str(k), 2), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Master-Node-'+str(5), 2),
            # plane
            (Length_Y*-N1_H_dy+A*Length_X*-N1_H_dx, 'N'+str(nodes[0])+'-RP', 2),
            (Length_Y*-N2_H_dy+A*Length_X*-N2_H_dx, 'N'+str(nodes[1])+'-RP', 2),
            (Length_Y*-N3_H_dy+A*Length_X*-N3_H_dx, 'N'+str(nodes[2])+'-RP', 2),
            (Length_Y*-N4_H_dy+A*Length_X*-N4_H_dx, 'N'+str(nodes[3])+'-RP', 2),
            # bending
            (Length_Y*Zcoord*N1_d2y+B*Length_X*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_Y*Zcoord*NX1_d2y+B*Length_X*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_Y*Zcoord*NY1_d2y+B*Length_X*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_Y*Zcoord*N2_d2y+B*Length_X*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_Y*Zcoord*NX2_d2y+B*Length_X*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_Y*Zcoord*NY2_d2y+B*Length_X*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_Y*Zcoord*N3_d2y+B*Length_X*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_Y*Zcoord*NX3_d2y+B*Length_X*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_Y*Zcoord*NY3_d2y+B*Length_X*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_Y*Zcoord*N4_d2y+B*Length_X*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_Y*Zcoord*NX4_d2y+B*Length_X*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_Y*Zcoord*NY4_d2y+B*Length_X*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))







        # edges 2-1,3-1,4-1
        for k in range(2,5): 
            maxIndex = len(mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge'+str(1)].nodes)
            for kk in range(maxIndex): 
                P1=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(1)+'-'+str(kk)]
                P2=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(k)+'-'+str(kk)]
                Xcoord_P1=P1.nodes[0].coordinates[0]
                Ycoord_P1=P1.nodes[0].coordinates[1]
                Zcoord_P1=P1.nodes[0].coordinates[2]
                Xcoord_P2=P2.nodes[0].coordinates[0]
                Ycoord_P2=P2.nodes[0].coordinates[1]
                Zcoord_P2=P2.nodes[0].coordinates[2]
                Length_X=float(abs(Xcoord_P2-Xcoord_P1))
                Length_Y=float(abs(Ycoord_P2-Ycoord_P1))
                Zcoord=(Zcoord_P1+Zcoord_P2)/2*Scale
                mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge-'+str(k)+'-1-'+str(kk)+'-on-u', 
                terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(k)+'-'+str(kk), 1), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(1)+'-'+str(kk), 1),
                # plane
                (Length_X*-N1_H_dx+A*Length_Y*-N1_H_dy, 'N'+str(nodes[0])+'-RP', 1),
                (Length_X*-N2_H_dx+A*Length_Y*-N2_H_dy, 'N'+str(nodes[1])+'-RP', 1),
                (Length_X*-N3_H_dx+A*Length_Y*-N3_H_dy, 'N'+str(nodes[2])+'-RP', 1),
                (Length_X*-N4_H_dx+A*Length_Y*-N4_H_dy, 'N'+str(nodes[3])+'-RP', 1),
                # bending
                (Length_X*Zcoord*N1_d2x+C*Length_Y*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
                (Length_X*Zcoord*NX1_d2x+C*Length_Y*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
                (Length_X*Zcoord*NY1_d2x+C*Length_Y*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
                (Length_X*Zcoord*N2_d2x+C*Length_Y*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
                (Length_X*Zcoord*NX2_d2x+C*Length_Y*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
                (Length_X*Zcoord*NY2_d2x+C*Length_Y*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
                (Length_X*Zcoord*N3_d2x+C*Length_Y*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
                (Length_X*Zcoord*NX3_d2x+C*Length_Y*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
                (Length_X*Zcoord*NY3_d2x+C*Length_Y*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
                (Length_X*Zcoord*N4_d2x+C*Length_Y*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
                (Length_X*Zcoord*NX4_d2x+C*Length_Y*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
                (Length_X*Zcoord*NY4_d2x+C*Length_Y*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))


                mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge-'+str(k)+'-1-'+str(kk)+'-on-v', 
                terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(k)+'-'+str(kk), 2), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(1)+'-'+str(kk), 2),
                # plane
                (Length_Y*-N1_H_dy+A*Length_X*-N1_H_dx, 'N'+str(nodes[0])+'-RP', 2),
                (Length_Y*-N2_H_dy+A*Length_X*-N2_H_dx, 'N'+str(nodes[1])+'-RP', 2),
                (Length_Y*-N3_H_dy+A*Length_X*-N3_H_dx, 'N'+str(nodes[2])+'-RP', 2),
                (Length_Y*-N4_H_dy+A*Length_X*-N4_H_dx, 'N'+str(nodes[3])+'-RP', 2),
                # bending
                (Length_Y*Zcoord*N1_d2y+B*Length_X*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
                (Length_Y*Zcoord*NX1_d2y+B*Length_X*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
                (Length_Y*Zcoord*NY1_d2y+B*Length_X*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
                (Length_Y*Zcoord*N2_d2y+B*Length_X*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
                (Length_Y*Zcoord*NX2_d2y+B*Length_X*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
                (Length_Y*Zcoord*NY2_d2y+B*Length_X*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
                (Length_Y*Zcoord*N3_d2y+B*Length_X*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
                (Length_Y*Zcoord*NX3_d2y+B*Length_X*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
                (Length_Y*Zcoord*NY3_d2y+B*Length_X*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
                (Length_Y*Zcoord*N4_d2y+B*Length_X*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
                (Length_Y*Zcoord*NX4_d2y+B*Length_X*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
                (Length_Y*Zcoord*NY4_d2y+B*Length_X*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))








        #edge 6-5
        maxIndex = len(mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge'+str(5)].nodes)
        for kk in range(maxIndex):
            P1=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(5)+'-'+str(kk)]
            P2=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(6)+'-'+str(kk)]
            Xcoord_P1=P1.nodes[0].coordinates[0]
            Ycoord_P1=P1.nodes[0].coordinates[1]
            Zcoord_P1=P1.nodes[0].coordinates[2]
            Xcoord_P2=P2.nodes[0].coordinates[0]
            Ycoord_P2=P2.nodes[0].coordinates[1]
            Zcoord_P2=P2.nodes[0].coordinates[2]
            Length_X=float(abs(Xcoord_P2-Xcoord_P1))
            Length_Y=float(abs(Ycoord_P2-Ycoord_P1))
            Zcoord=(Zcoord_P1+Zcoord_P2)/2*Scale
            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge-6-5-'+str(kk)+'-on-u', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(6)+'-'+str(kk), 1), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(5)+'-'+str(kk), 1),
            # plane
            (Length_X*-N1_H_dx+A*Length_Y*-N1_H_dy, 'N'+str(nodes[0])+'-RP', 1),
            (Length_X*-N2_H_dx+A*Length_Y*-N2_H_dy, 'N'+str(nodes[1])+'-RP', 1),
            (Length_X*-N3_H_dx+A*Length_Y*-N3_H_dy, 'N'+str(nodes[2])+'-RP', 1),
            (Length_X*-N4_H_dx+A*Length_Y*-N4_H_dy, 'N'+str(nodes[3])+'-RP', 1),
            # bending
            (Length_X*Zcoord*N1_d2x+C*Length_Y*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_X*Zcoord*NX1_d2x+C*Length_Y*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_X*Zcoord*NY1_d2x+C*Length_Y*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_X*Zcoord*N2_d2x+C*Length_Y*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_X*Zcoord*NX2_d2x+C*Length_Y*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_X*Zcoord*NY2_d2x+C*Length_Y*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_X*Zcoord*N3_d2x+C*Length_Y*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_X*Zcoord*NX3_d2x+C*Length_Y*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_X*Zcoord*NY3_d2x+C*Length_Y*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_X*Zcoord*N4_d2x+C*Length_Y*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_X*Zcoord*NX4_d2x+C*Length_Y*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_X*Zcoord*NY4_d2x+C*Length_Y*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))

            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge-6-5-'+str(kk)+'-on-v', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(6)+'-'+str(kk), 2), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(5)+'-'+str(kk), 2),
            # plane
            (Length_Y*-N1_H_dy+A*Length_X*-N1_H_dx, 'N'+str(nodes[0])+'-RP', 2),
            (Length_Y*-N2_H_dy+A*Length_X*-N2_H_dx, 'N'+str(nodes[1])+'-RP', 2),
            (Length_Y*-N3_H_dy+A*Length_X*-N3_H_dx, 'N'+str(nodes[2])+'-RP', 2),
            (Length_Y*-N4_H_dy+A*Length_X*-N4_H_dx, 'N'+str(nodes[3])+'-RP', 2),
            # bending
            (Length_Y*Zcoord*N1_d2y+B*Length_X*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_Y*Zcoord*NX1_d2y+B*Length_X*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_Y*Zcoord*NY1_d2y+B*Length_X*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_Y*Zcoord*N2_d2y+B*Length_X*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_Y*Zcoord*NX2_d2y+B*Length_X*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_Y*Zcoord*NY2_d2y+B*Length_X*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_Y*Zcoord*N3_d2y+B*Length_X*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_Y*Zcoord*NX3_d2y+B*Length_X*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_Y*Zcoord*NY3_d2y+B*Length_X*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_Y*Zcoord*N4_d2y+B*Length_X*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_Y*Zcoord*NX4_d2y+B*Length_X*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_Y*Zcoord*NY4_d2y+B*Length_X*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))








        #edge 7-8
        maxIndex = len(mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge'+str(8)].nodes)
        for kk in range(maxIndex):
            P1=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(8)+'-'+str(kk)]
            P2=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(7)+'-'+str(kk)]
            Xcoord_P1=P1.nodes[0].coordinates[0]
            Ycoord_P1=P1.nodes[0].coordinates[1]
            Zcoord_P1=P1.nodes[0].coordinates[2]
            Xcoord_P2=P2.nodes[0].coordinates[0]
            Ycoord_P2=P2.nodes[0].coordinates[1]
            Zcoord_P2=P2.nodes[0].coordinates[2]
            Length_X=float(abs(Xcoord_P2-Xcoord_P1))
            Length_Y=float(abs(Ycoord_P2-Ycoord_P1))
            Zcoord=(Zcoord_P1+Zcoord_P2)/2*Scale
            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge-7-8-'+str(kk)+'-on-u', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(7)+'-'+str(kk), 1), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(8)+'-'+str(kk), 1),
            # plane
            (Length_X*-N1_H_dx+A*Length_Y*-N1_H_dy, 'N'+str(nodes[0])+'-RP', 1),
            (Length_X*-N2_H_dx+A*Length_Y*-N2_H_dy, 'N'+str(nodes[1])+'-RP', 1),
            (Length_X*-N3_H_dx+A*Length_Y*-N3_H_dy, 'N'+str(nodes[2])+'-RP', 1),
            (Length_X*-N4_H_dx+A*Length_Y*-N4_H_dy, 'N'+str(nodes[3])+'-RP', 1),
            # bending
            (Length_X*Zcoord*N1_d2x+C*Length_Y*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_X*Zcoord*NX1_d2x+C*Length_Y*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_X*Zcoord*NY1_d2x+C*Length_Y*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_X*Zcoord*N2_d2x+C*Length_Y*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_X*Zcoord*NX2_d2x+C*Length_Y*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_X*Zcoord*NY2_d2x+C*Length_Y*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_X*Zcoord*N3_d2x+C*Length_Y*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_X*Zcoord*NX3_d2x+C*Length_Y*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_X*Zcoord*NY3_d2x+C*Length_Y*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_X*Zcoord*N4_d2x+C*Length_Y*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_X*Zcoord*NX4_d2x+C*Length_Y*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_X*Zcoord*NY4_d2x+C*Length_Y*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))


            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge-7-8-'+str(kk)+'-on-v', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(7)+'-'+str(kk), 2), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(8)+'-'+str(kk), 2),
            # plane
            (Length_Y*-N1_H_dy+A*Length_X*-N1_H_dx, 'N'+str(nodes[0])+'-RP', 2),
            (Length_Y*-N2_H_dy+A*Length_X*-N2_H_dx, 'N'+str(nodes[1])+'-RP', 2),
            (Length_Y*-N3_H_dy+A*Length_X*-N3_H_dx, 'N'+str(nodes[2])+'-RP', 2),
            (Length_Y*-N4_H_dy+A*Length_X*-N4_H_dx, 'N'+str(nodes[3])+'-RP', 2),
            # bending
            (Length_Y*Zcoord*N1_d2y+B*Length_X*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_Y*Zcoord*NX1_d2y+B*Length_X*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_Y*Zcoord*NY1_d2y+B*Length_X*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_Y*Zcoord*N2_d2y+B*Length_X*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_Y*Zcoord*NX2_d2y+B*Length_X*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_Y*Zcoord*NY2_d2y+B*Length_X*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_Y*Zcoord*N3_d2y+B*Length_X*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_Y*Zcoord*NX3_d2y+B*Length_X*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_Y*Zcoord*NY3_d2y+B*Length_X*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_Y*Zcoord*N4_d2y+B*Length_X*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_Y*Zcoord*NX4_d2y+B*Length_X*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_Y*Zcoord*NY4_d2y+B*Length_X*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))







        #edge 10-9
        maxIndex = len(mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge'+str(9)].nodes)
        for kk in range(maxIndex):
            P1=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(9)+'-'+str(kk)]
            P2=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(10)+'-'+str(kk)]
            Xcoord_P1=P1.nodes[0].coordinates[0]
            Ycoord_P1=P1.nodes[0].coordinates[1]
            Zcoord_P1=P1.nodes[0].coordinates[2]
            Xcoord_P2=P2.nodes[0].coordinates[0]
            Ycoord_P2=P2.nodes[0].coordinates[1]
            Zcoord_P2=P2.nodes[0].coordinates[2]
            Length_X=float(abs(Xcoord_P2-Xcoord_P1))
            Length_Y=float(abs(Ycoord_P2-Ycoord_P1))
            Zcoord=(Zcoord_P1+Zcoord_P2)/2*Scale
            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge-10-9-'+str(kk)+'-on-u', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(10)+'-'+str(kk), 1), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(9)+'-'+str(kk), 1),
            # plane
            (Length_X*-N1_H_dx+A*Length_Y*-N1_H_dy, 'N'+str(nodes[0])+'-RP', 1),
            (Length_X*-N2_H_dx+A*Length_Y*-N2_H_dy, 'N'+str(nodes[1])+'-RP', 1),
            (Length_X*-N3_H_dx+A*Length_Y*-N3_H_dy, 'N'+str(nodes[2])+'-RP', 1),
            (Length_X*-N4_H_dx+A*Length_Y*-N4_H_dy, 'N'+str(nodes[3])+'-RP', 1),
            # bending
            (Length_X*Zcoord*N1_d2x+C*Length_Y*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_X*Zcoord*NX1_d2x+C*Length_Y*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_X*Zcoord*NY1_d2x+C*Length_Y*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_X*Zcoord*N2_d2x+C*Length_Y*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_X*Zcoord*NX2_d2x+C*Length_Y*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_X*Zcoord*NY2_d2x+C*Length_Y*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_X*Zcoord*N3_d2x+C*Length_Y*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_X*Zcoord*NX3_d2x+C*Length_Y*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_X*Zcoord*NY3_d2x+C*Length_Y*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_X*Zcoord*N4_d2x+C*Length_Y*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_X*Zcoord*NX4_d2x+C*Length_Y*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_X*Zcoord*NY4_d2x+C*Length_Y*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))


            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge-10-9-'+str(kk)+'-on-v', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(10)+'-'+str(kk), 2), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(9)+'-'+str(kk), 2),
            # plane
            (Length_Y*-N1_H_dy+A*Length_X*-N1_H_dx, 'N'+str(nodes[0])+'-RP', 2),
            (Length_Y*-N2_H_dy+A*Length_X*-N2_H_dx, 'N'+str(nodes[1])+'-RP', 2),
            (Length_Y*-N3_H_dy+A*Length_X*-N3_H_dx, 'N'+str(nodes[2])+'-RP', 2),
            (Length_Y*-N4_H_dy+A*Length_X*-N4_H_dx, 'N'+str(nodes[3])+'-RP', 2),
            # bending
            (Length_Y*Zcoord*N1_d2y+B*Length_X*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_Y*Zcoord*NX1_d2y+B*Length_X*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_Y*Zcoord*NY1_d2y+B*Length_X*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_Y*Zcoord*N2_d2y+B*Length_X*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_Y*Zcoord*NX2_d2y+B*Length_X*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_Y*Zcoord*NY2_d2y+B*Length_X*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_Y*Zcoord*N3_d2y+B*Length_X*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_Y*Zcoord*NX3_d2y+B*Length_X*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_Y*Zcoord*NY3_d2y+B*Length_X*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_Y*Zcoord*N4_d2y+B*Length_X*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_Y*Zcoord*NX4_d2y+B*Length_X*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_Y*Zcoord*NY4_d2y+B*Length_X*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))







        #edge 11-12
        maxIndex = len(mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-edge'+str(12)].nodes)
        for kk in range(maxIndex):
            P1=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(12)+'-'+str(kk)]
            P2=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(11)+'-'+str(kk)]
            Xcoord_P1=P1.nodes[0].coordinates[0]
            Ycoord_P1=P1.nodes[0].coordinates[1]
            Zcoord_P1=P1.nodes[0].coordinates[2]
            Xcoord_P2=P2.nodes[0].coordinates[0]
            Ycoord_P2=P2.nodes[0].coordinates[1]
            Zcoord_P2=P2.nodes[0].coordinates[2]
            Length_X=float(abs(Xcoord_P2-Xcoord_P1))
            Length_Y=float(abs(Ycoord_P2-Ycoord_P1))
            Zcoord=(Zcoord_P1+Zcoord_P2)/2*Scale
            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge-11-12-'+str(kk)+'-on-u', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(11)+'-'+str(kk), 1), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(12)+'-'+str(kk), 1),
            # plane
            (Length_X*-N1_H_dx+A*Length_Y*-N1_H_dy, 'N'+str(nodes[0])+'-RP', 1),
            (Length_X*-N2_H_dx+A*Length_Y*-N2_H_dy, 'N'+str(nodes[1])+'-RP', 1),
            (Length_X*-N3_H_dx+A*Length_Y*-N3_H_dy, 'N'+str(nodes[2])+'-RP', 1),
            (Length_X*-N4_H_dx+A*Length_Y*-N4_H_dy, 'N'+str(nodes[3])+'-RP', 1),
            # bending
            (Length_X*Zcoord*N1_d2x+C*Length_Y*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_X*Zcoord*NX1_d2x+C*Length_Y*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_X*Zcoord*NY1_d2x+C*Length_Y*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_X*Zcoord*N2_d2x+C*Length_Y*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_X*Zcoord*NX2_d2x+C*Length_Y*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_X*Zcoord*NY2_d2x+C*Length_Y*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_X*Zcoord*N3_d2x+C*Length_Y*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_X*Zcoord*NX3_d2x+C*Length_Y*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_X*Zcoord*NY3_d2x+C*Length_Y*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_X*Zcoord*N4_d2x+C*Length_Y*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_X*Zcoord*NX4_d2x+C*Length_Y*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_X*Zcoord*NY4_d2x+C*Length_Y*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))


            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge-11-12-'+str(kk)+'-on-v', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(11)+'-'+str(kk), 2), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-Edge'+str(12)+'-'+str(kk), 2),
            # plane
            (Length_Y*-N1_H_dy+A*Length_X*-N1_H_dx, 'N'+str(nodes[0])+'-RP', 2),
            (Length_Y*-N2_H_dy+A*Length_X*-N2_H_dx, 'N'+str(nodes[1])+'-RP', 2),
            (Length_Y*-N3_H_dy+A*Length_X*-N3_H_dx, 'N'+str(nodes[2])+'-RP', 2),
            (Length_Y*-N4_H_dy+A*Length_X*-N4_H_dx, 'N'+str(nodes[3])+'-RP', 2),
            # bending
            (Length_Y*Zcoord*N1_d2y+B*Length_X*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_Y*Zcoord*NX1_d2y+B*Length_X*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_Y*Zcoord*NY1_d2y+B*Length_X*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_Y*Zcoord*N2_d2y+B*Length_X*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_Y*Zcoord*NX2_d2y+B*Length_X*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_Y*Zcoord*NY2_d2y+B*Length_X*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_Y*Zcoord*N3_d2y+B*Length_X*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_Y*Zcoord*NX3_d2y+B*Length_X*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_Y*Zcoord*NY3_d2y+B*Length_X*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_Y*Zcoord*N4_d2y+B*Length_X*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_Y*Zcoord*NX4_d2y+B*Length_X*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_Y*Zcoord*NY4_d2y+B*Length_X*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))








        #Face XP-XN
        maxIndex = len(mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Face-XN'].nodes)
        for k in range(maxIndex):
            P1=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceXN-'+str(k)]
            P2=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceXP-'+str(k)]
            Xcoord_P1=P1.nodes[0].coordinates[0]
            Ycoord_P1=P1.nodes[0].coordinates[1]
            Zcoord_P1=P1.nodes[0].coordinates[2]
            Xcoord_P2=P2.nodes[0].coordinates[0]
            Ycoord_P2=P2.nodes[0].coordinates[1]
            Zcoord_P2=P2.nodes[0].coordinates[2]
            Length_X=float(abs(Xcoord_P2-Xcoord_P1))
            Length_Y=float(abs(Ycoord_P2-Ycoord_P1))
            Zcoord=(Zcoord_P1+Zcoord_P2)/2*Scale
            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceXP-XN-'+str(k)+'-on-u', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceXP-'+str(k), 1), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceXN-'+str(k), 1),
            # plane
            (Length_X*-N1_H_dx+A*Length_Y*-N1_H_dy, 'N'+str(nodes[0])+'-RP', 1),
            (Length_X*-N2_H_dx+A*Length_Y*-N2_H_dy, 'N'+str(nodes[1])+'-RP', 1),
            (Length_X*-N3_H_dx+A*Length_Y*-N3_H_dy, 'N'+str(nodes[2])+'-RP', 1),
            (Length_X*-N4_H_dx+A*Length_Y*-N4_H_dy, 'N'+str(nodes[3])+'-RP', 1),
            # bending
            (Length_X*Zcoord*N1_d2x+C*Length_Y*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_X*Zcoord*NX1_d2x+C*Length_Y*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_X*Zcoord*NY1_d2x+C*Length_Y*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_X*Zcoord*N2_d2x+C*Length_Y*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_X*Zcoord*NX2_d2x+C*Length_Y*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_X*Zcoord*NY2_d2x+C*Length_Y*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_X*Zcoord*N3_d2x+C*Length_Y*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_X*Zcoord*NX3_d2x+C*Length_Y*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_X*Zcoord*NY3_d2x+C*Length_Y*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_X*Zcoord*N4_d2x+C*Length_Y*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_X*Zcoord*NX4_d2x+C*Length_Y*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_X*Zcoord*NY4_d2x+C*Length_Y*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))

            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceXP-XN-'+str(k)+'-on-v', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceXP-'+str(k), 2), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceXN-'+str(k), 2),
            # plane
            (Length_Y*-N1_H_dy+A*Length_X*-N1_H_dx, 'N'+str(nodes[0])+'-RP', 2),
            (Length_Y*-N2_H_dy+A*Length_X*-N2_H_dx, 'N'+str(nodes[1])+'-RP', 2),
            (Length_Y*-N3_H_dy+A*Length_X*-N3_H_dx, 'N'+str(nodes[2])+'-RP', 2),
            (Length_Y*-N4_H_dy+A*Length_X*-N4_H_dx, 'N'+str(nodes[3])+'-RP', 2),
            # bending
            (Length_Y*Zcoord*N1_d2y+B*Length_X*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_Y*Zcoord*NX1_d2y+B*Length_X*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_Y*Zcoord*NY1_d2y+B*Length_X*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_Y*Zcoord*N2_d2y+B*Length_X*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_Y*Zcoord*NX2_d2y+B*Length_X*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_Y*Zcoord*NY2_d2y+B*Length_X*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_Y*Zcoord*N3_d2y+B*Length_X*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_Y*Zcoord*NX3_d2y+B*Length_X*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_Y*Zcoord*NY3_d2y+B*Length_X*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_Y*Zcoord*N4_d2y+B*Length_X*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_Y*Zcoord*NX4_d2y+B*Length_X*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_Y*Zcoord*NY4_d2y+B*Length_X*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))







        #Face YP-YN
        maxIndex = len(mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-Face-YN'].nodes)
        for k in range(maxIndex):
            P1=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceYN-'+str(k)]
            P2=mdb.models[modelName].rootAssembly.sets['Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceYP-'+str(k)]
            Xcoord_P1=P1.nodes[0].coordinates[0]
            Ycoord_P1=P1.nodes[0].coordinates[1]
            Zcoord_P1=P1.nodes[0].coordinates[2]
            Xcoord_P2=P2.nodes[0].coordinates[0]
            Ycoord_P2=P2.nodes[0].coordinates[1]
            Zcoord_P2=P2.nodes[0].coordinates[2]
            Length_X=float(abs(Xcoord_P2-Xcoord_P1))
            Length_Y=float(abs(Ycoord_P2-Ycoord_P1))
            Zcoord=(Zcoord_P1+Zcoord_P2)/2*Scale
            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceYP-YN-'+str(k)+'-on-u', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceYP-'+str(k), 1), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceYN-'+str(k), 1),
            # plane
            (Length_X*-N1_H_dx+A*Length_Y*-N1_H_dy, 'N'+str(nodes[0])+'-RP', 1),
            (Length_X*-N2_H_dx+A*Length_Y*-N2_H_dy, 'N'+str(nodes[1])+'-RP', 1),
            (Length_X*-N3_H_dx+A*Length_Y*-N3_H_dy, 'N'+str(nodes[2])+'-RP', 1),
            (Length_X*-N4_H_dx+A*Length_Y*-N4_H_dy, 'N'+str(nodes[3])+'-RP', 1),
            # bending
            (Length_X*Zcoord*N1_d2x+C*Length_Y*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_X*Zcoord*NX1_d2x+C*Length_Y*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_X*Zcoord*NY1_d2x+C*Length_Y*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_X*Zcoord*N2_d2x+C*Length_Y*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_X*Zcoord*NX2_d2x+C*Length_Y*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_X*Zcoord*NY2_d2x+C*Length_Y*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_X*Zcoord*N3_d2x+C*Length_Y*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_X*Zcoord*NX3_d2x+C*Length_Y*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_X*Zcoord*NY3_d2x+C*Length_Y*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_X*Zcoord*N4_d2x+C*Length_Y*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_X*Zcoord*NX4_d2x+C*Length_Y*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_X*Zcoord*NY4_d2x+C*Length_Y*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))


            mdb.models[modelName].Equation(name = 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceYP-YN-'+str(k)+'-on-v', 
            terms = ((1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceYP-'+str(k), 2), (-1, 'Ele-'+str(ele)+'-GP-'+str(jj)+'-FaceYN-'+str(k), 2),
            # plane
            (Length_Y*-N1_H_dy+A*Length_X*-N1_H_dx, 'N'+str(nodes[0])+'-RP', 2),
            (Length_Y*-N2_H_dy+A*Length_X*-N2_H_dx, 'N'+str(nodes[1])+'-RP', 2),
            (Length_Y*-N3_H_dy+A*Length_X*-N3_H_dx, 'N'+str(nodes[2])+'-RP', 2),
            (Length_Y*-N4_H_dy+A*Length_X*-N4_H_dx, 'N'+str(nodes[3])+'-RP', 2),
            # bending
            (Length_Y*Zcoord*N1_d2y+B*Length_X*Zcoord*N1_dxdy, 'N'+str(nodes[0])+'-RP', 3), 
            (Length_Y*Zcoord*NX1_d2y+B*Length_X*Zcoord*NX1_dxdy, 'N'+str(nodes[0])+'-RP', 4),
            (Length_Y*Zcoord*NY1_d2y+B*Length_X*Zcoord*NY1_dxdy, 'N'+str(nodes[0])+'-RP', 5),
            (Length_Y*Zcoord*N2_d2y+B*Length_X*Zcoord*N2_dxdy, 'N'+str(nodes[1])+'-RP', 3), 
            (Length_Y*Zcoord*NX2_d2y+B*Length_X*Zcoord*NX2_dxdy, 'N'+str(nodes[1])+'-RP', 4),
            (Length_Y*Zcoord*NY2_d2y+B*Length_X*Zcoord*NY2_dxdy, 'N'+str(nodes[1])+'-RP', 5),
            (Length_Y*Zcoord*N3_d2y+B*Length_X*Zcoord*N3_dxdy, 'N'+str(nodes[2])+'-RP', 3), 
            (Length_Y*Zcoord*NX3_d2y+B*Length_X*Zcoord*NX3_dxdy, 'N'+str(nodes[2])+'-RP', 4),
            (Length_Y*Zcoord*NY3_d2y+B*Length_X*Zcoord*NY3_dxdy, 'N'+str(nodes[2])+'-RP', 5),
            (Length_Y*Zcoord*N4_d2y+B*Length_X*Zcoord*N4_dxdy, 'N'+str(nodes[3])+'-RP', 3), 
            (Length_Y*Zcoord*NX4_d2y+B*Length_X*Zcoord*NX4_dxdy, 'N'+str(nodes[3])+'-RP', 4),
            (Length_Y*Zcoord*NY4_d2y+B*Length_X*Zcoord*NY4_dxdy, 'N'+str(nodes[3])+'-RP', 5)))








        print 'Ele-',str(ele),'-GP-',str(jj),'Finish lalalalalallalalala'
