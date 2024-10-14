#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import os
from calculix import Calculix,Mesh

# get the legnth scale if provided
parser = argparse.ArgumentParser(description='Generates the Finite-element mesh..')
parser.add_argument('-L','--lengthscale', help='Length scale of the body')
parser.add_argument('-N','--npoints', help='number of points along the length')
parser.add_argument('-O','--order', help='Order of the finite-element mesh')
args = parser.parse_args()

# set default variables
L = 1. if(not args.lengthscale) else float(args.lengthscale)
o = 1 if(not args.order) else int(args.order)
N = 16 if(not args.npoints) else int(args.npoints)

# # generate temporary mesh file
replace("template.geo","LIMO.geo",{"NX": str(N)})
os.system("gmsh LIMO.geo -order %d -2 -o init.inp" % o)

def roll(x,L=10):
    r = L/(2*np.pi)
    return [r*np.cos(2*np.pi*x[0]/L),x[1],r*np.sin(2*np.pi*x[0]/L)]

# mesh generation
mesh = Mesh("init.inp")
# duplicate inside the pouch
mesh.duplicate_srf([1],[2,3])
# duplicate the complete mesh, except the outer edges
# mesh.mirror_srf([0,1,2],[0,1],"z",[0,0,1e-3])
# select BC nodes
# mesh.make_nodeset([3,4,5],lambda x: x[1]==0.0,"EDGE_2")
mesh.save_lines([0,1],"OUTER","outer.nam")
mesh.save_lines([4,5,6],"EDGE","edge.nam")
mesh.save("geom.inp")

# material properties
E = 10000
nu = 0.0
density = 800

# generate calculix.inp file
calculix = Calculix()
calculix.include(["outer.nam","edge.nam"])
calculix.set_material(E,nu,density,alpha=0.65,beta=0.)
calculix.set_bc("EDGE", 2, 2, 0)
calculix.write("*AMPLITUDE, NAME=A1\n 0.,0., 20., 1. \n")
calculix.set_thickness(1e-4)
calculix.write("*STEP, NLGEOM, INC=1000000\n*DYNAMIC\n 0.05, 20., 0.0001, 0.1\n")
calculix.write("*DLOAD, AMPLITUDE=A1\n SRF_1, P1, 0.1\n SRF_2, P1, 0.1\n SRF_3, P1, 0.1\n")
calculix.close("\n*EL FILE\n E, ME, S")
#
