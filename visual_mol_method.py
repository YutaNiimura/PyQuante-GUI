import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pylab as p
from mpl_toolkits.mplot3d import Axes3D
from PyQuante.Molecule import Molecule
from PyQuante.hartree_fock import rhf
from PyQuante.Ints import getbasis
import numpy as np
import math
from visual import *
from PyQuante import configure_output
import csv

def default(self):
    self.coor1_x.SetValue("0.70000000")
    self.coor1_y.SetValue("0.00000000")
    self.coor1_z.SetValue("0.00000000")
    self.coor2_x.SetValue("-0.70000000")
    self.coor2_y.SetValue("0.00000000")
    self.coor2_z.SetValue("0.00000000")
    self.mol_box.SetValue("H2")
    self.Bf_box.SetValue("sto-3g")
    #self.orb_bound_box.SetValue("boundary")
    self.z_pos.SetValue("0.0")

def calculate(self):
    
    if self.log == 1:
        f = open('calc.log','w')
        f.close()
        configure_output("calc.log")
    x_1 = float(self.coor1_x.GetValue())
    y_1 = float(self.coor1_y.GetValue())
    z_1 = float(self.coor1_z.GetValue())
    x_2 = float(self.coor2_x.GetValue())
    y_2 = float(self.coor2_y.GetValue())
    z_2 = float(self.coor2_z.GetValue())
 
    mol = Molecule('h2',[(self.atom_no_1,(x_1,y_1,z_1)),
                          (self.atom_no_2,(x_2,y_2,z_2))])
    print mol
    bf_name = self.Bf_box.GetValue()
    
        #Why?
    if bf_name == "sto-3g":
        bf_name = "sto-3g"
    elif bf_name == "6-31g":
        bf_name = "6-31g"
    elif bf_name == "3-21g":
        bf_name = "3-21g"
    elif bf_name == "6-31g**":
        bf_name = "6-31g**"
    elif bf_name == "6-31g**++":
        bf_name = "6-31g**++"
    elif bf_name == "6-311g**":
        bf_name = "6-311g**"
    elif bf_name == "6-311g++(2d,2p)":
        bf_name = "6-311g++(2d,2p)"
    elif bf_name == "6-311g++(3d,3p)":
        bf_name = "6-311g++(3d,3p)"
    elif bf_name == "6-311g++(3df,3pd)":
        bf_name = "6-311g++(3df,3pd)"
    elif bf_name == "sto-6g":
        bf_name = "sto-6g"
    elif bf_name == "lacvp":
        bf_name = "lacvp"
    elif bf_name == "cc-pvdz":
        bf_name = "cc-pvdz"
    elif bf_name == "cc-pvtz":
        bf_name = "cc-pvtz"
    elif bf_name == "dzvp":
        bf_name = "dzvp"
        
    en,orbe,orbs = rhf(mol,basis_data=bf_name)
    if self.log == 1:
        with open('calc.log','a') as f_handle:
            f_handle.write("\norbs is\n")
            np.savetxt(f_handle,orbs)
            f_handle.write("\norbe is\n")
            np.savetxt(f_handle,orbe)
    
    self.console_box.SetValue('energy=%f' % en)

def visualize(self):
    if self.log == 1:
        f = open('calc.log','w')
        f.close()
        configure_output("calc.log")
    x_1 = float(self.coor1_x.GetValue())
    y_1 = float(self.coor1_y.GetValue())
    z_1 = float(self.coor1_z.GetValue())
    x_2 = float(self.coor2_x.GetValue())
    y_2 = float(self.coor2_y.GetValue())
    z_2 = float(self.coor2_z.GetValue())
    z_pos = float(self.z_pos.GetValue())
    mol = Molecule('h2',[(self.atom_no_1,(y_1,x_1,z_1)),
                          (self.atom_no_2,(y_2,x_2,z_2))],
                   units='Angstrom')
    bf_name = self.Bf_box.GetValue()
    
        #Why?
    if bf_name == "sto-3g":
        bf_name = "sto-3g"
    elif bf_name == "6-31g":
        bf_name = "6-31g"
    elif bf_name == "3-21g":
        bf_name = "3-21g"
    elif bf_name == "6-31g**":
        bf_name = "6-31g**"
    elif bf_name == "6-31g**++":
        bf_name = "6-31g**++"
    elif bf_name == "6-311g**":
        bf_name = "6-311g**"
    elif bf_name == "6-311g++(2d,2p)":
        bf_name = "6-311g++(2d,2p)"
    elif bf_name == "6-311g++(3d,3p)":
        bf_name = "6-311g++(3d,3p)"
    elif bf_name == "6-311g++(3df,3pd)":
        bf_name = "6-311g++(3df,3pd)"
    elif bf_name == "sto-6g":
        bf_name = "sto-6g"
    elif bf_name == "lacvp":
        bf_name = "lacvp"
    elif bf_name == "cc-pvdz":
        bf_name = "cc-pvdz"
    elif bf_name == "cc-pvtz":
        bf_name = "cc-pvtz"
    elif bf_name == "dzvp":
        bf_name = "dzvp"
        
    en,orbe,orbs = rhf(mol,basis_data=bf_name)
    self.console_box.SetValue('energy=%f' % en)

    bfs = getbasis(mol,bf_name)

    delta = 0.1
    c_range = 5.0
#setting for visualize
    x = np.arange(-1*c_range,c_range,delta)
    y = np.arange(-1*c_range,c_range,delta)
    X,Y = p.meshgrid(x,y)
    Z = np.zeros((len(X),len(Y)))
    Z_2 = np.zeros((len(X),len(Y)))
    
#print C matrix
#print "C=",orbs

    #bounding = self.orb_bound_box.GetValue()
        
#calculate wave function
    for k,bf in enumerate(bfs.bfs):
        for i,x1 in enumerate(x):
            for j,y1 in enumerate(y):
            #calculate wave function
            #basic function multiply
                #if bounding == "boundary":
                Z[i,j] += bf.amp(x1,y1,z_pos) * orbs[k,0]
                #else:
                #    Z[i,j] += bf.amp(x1,y1,z_pos) * orbs[k,1]

    if self.log == 1:
        with open('calc.log','a') as f_handle:
            f_handle.write("\norbs is\n")
            np.savetxt(f_handle,orbs)
            f_handle.write("\norbe is\n")
            np.savetxt(f_handle,orbe)
                    
#####  visualize #####
    fig = p.figure()
    ax = Axes3D(fig)
    #ax.plot_surface(X,Y,Z)
    ax.plot_wireframe(X,Y,Z,color = 'b')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    p.show()

def setmolecule(self):
    mol = self.mol_box.GetValue()
    if mol == "H2":
        self.atom_no_1 = 1
        self.atom_no_2 = 1
    elif mol == "He2":
        self.atom_no_1 = 2
        self.atom_no_2 = 2
    elif mol == "OH":
        self.atom_no_1 = 8
        self.atom_no_2 = 1
