#
# builder for calculix input file: calculix.inp
#
import numpy as np

EOL = "\n"

def rotate(a,center,alpha,axis):
    x = a-center; b=a
    j=axis%3; k=(j+1)%3
    b[j] = np.cos(alpha)*x[j]-np.sin(alpha)*x[k]
    b[k] = np.cos(alpha)*x[k]+np.sin(alpha)*x[j]
    return b + center

def replace(template,target,dic):
    """
    Write target file by replacing the dic on the template
    """
    f1 = open(template,'r')
    f2 = open(target,'w')
    for text in f1:
        for i, j in dic.items():
            text = text.replace(i, j)
        f2.write(text)
    f1.close()
    f2.close()
    return None

class Calculix(object):

    def __init__(self):
        self.file = open("calculix.inp", "w")
        self.write("*INCLUDE, INPUT=geom.inp\n")

    def write(self, args):
        self.file.write(args)

    def include(self, fname):
        for name in fname:
            self.write("*INCLUDE, INPUT="+name+"\n")
        self.write(EOL)

    def set_bc(self, fname, bc_1, bc_2, value, func=""):
        txt = "*BOUNDARY"
        if(func != ""): txt += ",AMPLITUDE="+func
        self.write(txt+"\n")
        self.set_constraints(fname, bc_1, bc_2, value)

    def set_constraints(self, fname, bc_1, bc_2, value):
        self.write(" %s, %s, %s, %5.6f" % (fname, bc_1, bc_2, value))
        self.write(EOL)

    def set_amplitude(self, *uamplitude):
        for amplitude in uamplitude:
            self.write("*AMPLITUDE,NAME="+amplitude+",USER\n")
        self.write(EOL)

    def set_material(self, E, nu, density, alpha=0.0, beta=0.0):
        self.write("*MATERIAL,NAME=MEMBRANE\n")
        self.write("*ELASTIC\n %.1f, %.2f\n" % (E, nu))
        self.write("*DENSITY\n %.2f\n" % density)
        self.write("*DAMPING, ALPHA=%.2f, BETA=%.2f\n" % (alpha, beta))
        self.write(EOL)

    def set_thickness(self, thickness, nodal=False, offset=0.0):
        tmp = "*SHELL SECTION, MATERIAL=MEMBRANE, ELSET=PLATE"
        if(offset!=0.0): tmp +=", OFFSET=%.4f" % offset
        if(nodal):
            tmp+=",NODAL THICKNESS\n 0.1\n*INCLUDE, INPUT=thickness.nam\n"
        else:
            tmp += "\n %.10f\n" % thickness
        self.write(tmp)
        self.write(EOL)

    def add_node(self, nodes, nset="Nodes"):
        self.write("*NODE,NSET="+nset+"\n")
        for node in nodes:
            tmp = "%d" % node[0]
            for i in range(1,len(node)): tmp += ", %.4f" % node[i] 
            self.write(tmp)
        self.write(EOL)

    def add_surface(self, name, nodes):
        self.write("*SURFACE,NAME="+name+",TYPE=NODE\n")
        for node in nodes:
            self.write(" %d,\n" % node)
        self.write(EOL)

    def add_orientation(self, A, B, system="CYLINDRICAL", name="OR1"):
        self.write("*ORIENTATION,NAME="+name+",SYSTEM="+system+"\n")
        tmp = ""
        for i in range(3): tmp += "%.4f," % A[i]
        for i in range(3): tmp += "%.4f," % B[i]
        self.write(tmp[:-1])
        self.write(EOL+EOL)
    
    def add_coupling(self, ref, surface, orientation, name, kinematic=[]):
        self.write("*COUPLING,REF NODE="+str(ref)+",SURFACE="+surface)
        self.write(",ORIENTATION="+orientation+",CONSTRAINT NAME="+name+"\n ")
        if len(kinematic)>0:
            for i in range(len(kinematic)): self.write("%d," % kinematic[i])
            self.write(EOL)

    def set_solver(self, dt, T=1e5, Ninc=1e10, alpha=-0.05):
        self.write("*STEP, NLGEOM, INC=%d\n" % Ninc)
        self.write("*DYNAMIC, ALPHA=%-3f\n" % alpha)
        self.write(" %.4f, %.1f\n" % (dt, T))
        self.write("*CONTROLS,PARAMETERS=TIME INCREMENTATION\n ,,,,,,,10,,\n\n")
        self.write("*CLOAD\n")
        self.set_constraints("Ninterface", 1, 3, 0.0)

    def close(self, *args, frequency=1):
        self.write("*NODE FILE,OUTPUT=3D,FREQUENCY=%d\n" % frequency)
        out = " U"
        for arg in args:
            out += ", "+arg
        self.write(out+"\n")
        self.write("*END STEP\n")
        self.file.close()



class Mesh:
    def __init__(self, filename):
        self.nodes = []
        self.lines = []
        self.surface = []

        # open the mesh file and remove first useless lines
        lines = open("init.inp", "r").readlines()[3:]

        # check where the good stuff is
        eof = len(lines); elem_id=[]; line_id=[]
        for k,line in enumerate(lines):
            if("**" in line): i = k
            if("type=T3D2" in line): line_id.append(k)
            if("type=CPS" in line): elem_id.append(k)
            if("*ELSET,ELSET=" in line): eof = k; break
        line_id.append(elem_id[0]); elem_id.append(eof)
        # read the mesh and populate the nodes and elements
        for l in range(i):
            # get floats
            ids,x,y,z = to_float(lines[l].strip().replace(",","").split())
            self.nodes.append([ids+1,x,y,z])
        # read each of the line/element and populate accordingly
        for i in range(len(line_id)-1):
            node_list = []
            for l in range(line_id[i]+1,line_id[i+1]):
                out = to_int(lines[l].strip().replace(",","").split())[1:]
                node_list += [ids+1 for ids in out]
            self.lines.append(np.unique(node_list))
        # read each of the element and populate accordingly
        for i in range(len(elem_id)-1):
            node_list = []
            for l in range(elem_id[i]+1,elem_id[i+1]):
                out = to_int(lines[l].strip().replace(",","").split())
                node_list.append([out[0]]+[ids+1 for ids in out[1:]])
            self.surface.append(node_list)
            nelem = sum([len(el) for el in self.surface])

        self.k = self.nodes[-1][0]
        self.last_srf = self.surface[-1][-1][0]+1
        print("Done reading mesh: %d nodes, %d surfaces, %d elements" % (len(self.nodes),len(self.surface),nelem) )

    def add_node(self, node_id, new_nodes):
        for node in self.nodes: 
            if node[0] == node_id: # get the node that has this ID
                return self.check_existing(node, new_nodes) # check if the node has already been duplicated
        return None
    
    def check_existing(self, node, new_nodes):
        all_nodes = new_nodes # gather all the nodes together
        for new_node in all_nodes:
            if np.allclose(node[1:],new_node[1:]): # the node has already been duplicated/exists
                return new_node[0] # return existing duplicate node id  
        new_nodes.append([self.k+1]+node[1:]); self.k+=1
        return self.k # duplicate node and change index
    
    # mirrors a surface with boundary curves (not duplicated) in a `direction` around a point `origin`
    def mirror_srf(self, surface_id, line_id, direction, origin):
        line_nodes = [] if line_id==None else np.concatenate([self.lines[i] for i in line_id]) if len(line_id)>0 else self.lines[line_id]
        idx = int(self.k) # how many nodes we have before duplication
        unique_nodes = []
        for srf in surface_id:
            new_srf = []
            for el in self.surface[srf]:
                connectivity = []
                for nds in el[1:]: # for each node in the element
                    if nds in line_nodes: # if the node is in the outer line, we don't duplicate
                        connectivity.append(int(nds)) 
                    else: # duploicate it
                        unique_nodes.append(int(nds))
                        connectivity.append(int(nds)+idx)
                new_srf.append([self.last_srf]+[cx for cx in connectivity])
                self.last_srf+=1
            self.surface.append(new_srf)
        unique_nodes = np.unique(unique_nodes) # unique node ids
        # duplicate each unique node
        new_nodes = []
        for node_id in unique_nodes:
            for node in self.nodes:
                if node[0] == node_id:
                    new_nodes.append([node[0]+idx]+self.reflect(node[1:],direction,origin))
        self.nodes+=new_nodes
        return None
    
    def duplicate_srf(self, surface_id, line_id):
        for srf in surface_id:
            new_nodes,new_srf = self._duplicate_srf(self.surface[srf], line_id)
            # append the new nodes and elements to the list
            self.nodes += new_nodes
            self.surface.append(new_srf)
        return None

    def _duplicate_srf(self, surface, line_id):
        new_nodes = []; new_surface = []
        line_nodes = [] if line_id==None else np.concatenate([self.lines[i] for i in line_id]) if len(line_id)>0 else self.lines[line_id]
        # duplicate all the elements
        for el in surface:
            connectivity=[] # connectivity of the duplicated element
            for nd in el[1:]: # loop over the nodes of the element we duplicate
                # if the node is on the seams, we should duplicate, but reuse
                if nd in line_nodes:
                    connectivity.append(int(nd))
                else: # if the node needs duplication, check that it is not already duplicated
                    new_id = self.add_node(nd,new_nodes)
                    connectivity.append(int(new_id))
            # generate the new elements connectivity
            new_surface.append([int(self.last_srf)]+[int(cx) for cx in connectivity])
            self.last_srf+=1
        return new_nodes,new_surface
    
    def reflect(self,nodes,dir,o):
        if(dir=="x"): return [-(nodes[0]-o[0]),nodes[1],nodes[2]]
        elif(dir=="y"): return [nodes[0],-(nodes[1]-o[1]),nodes[2]]
        else: return [nodes[0],nodes[1],-(nodes[2]-o[2])]

    def save_lines(self, line_id, name, filename="lines.nam"):
        output = open(filename, "w")
        output.write("*NSET, NSET=%s\n" % name)
        unique = np.unique(np.concatenate([self.lines[i] for i in line_id]) if len(line_id)>0 else self.lines[line_id])
        output.write("".join(["%d,\n" % i for i in unique]))
        output.close()
        return None
    
    """
        map the nodes using a function f
    """
    def map_nodes(self, f):
        for i in range(len(self.nodes)):
            self.nodes[i][1:] = f(self.nodes[i][1:])
        return None
    
    def get_node_from_id(self, node_id):
        for node in self.nodes:
            if node[0] == node_id: return node
        return None
        
    def make_nodeset(self, surface_id, func, NSET="EDGE"):
        output = open(NSET+".nam", "w")
        output.write("*NSET, NSET=%s\n" % NSET)
        done = []
        for srf in surface_id:
            for el in self.surface[srf]:
                for nd in el[1:]:
                    if nd not in done:
                        if func(self.get_node_from_id(nd)[1:])==True:
                            output.write("%d,\n" % nd); done.append(nd)
        output.close()

    def save(self, filename="geometry.inp"):
        output = open(filename, "w")
        output.write("*NODE, NSET=Ninterface\n")
        # write all the nodes
        for ids,x,y,z in self.nodes:
            output.write("%d, %.8f, %.8f, %.8f\n" % (ids,x,y,z))
        output.write("*ELEMENT, type=S%d, ELSET=PLATE\n" % (3*o) )
        # for each surface, write the elements
        for srf in self.surface:
            for el in srf:
                line =  "%d" % el[0]
                for i in range(1,len(el)): line += ", %d" % (el[i])
                output.write(line+"\n")
        # write ELSETS for the different surfaces
        for i,srf in enumerate(self.surface):
            output.write("*ELSET, ELSET=SRF_%d\n" % (i+1))
            for el in srf:
                output.write("%d,\n" % el[0])
        output.close()
        return None

def to_float(ls):
    return [float(ls[i]) for i in range(len(ls))]
def to_int(ls):
    return [int(ls[i]) for i in range(len(ls))]