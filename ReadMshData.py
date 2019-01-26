from __future__ import print_function
import numpy as np
import numpy
import struct

class GmshMesh(object):
    """This is a class for storing nodes and elements. Based on Gmsh.py

    Members:
    nodes -- A dict of the form { nodeID: [ xcoord, ycoord, zcoord] }
    elements -- A dict of the form { elemID: (type, [tags], [nodeIDs]) }

    Methods:
    read([file]) -- Parse a Gmsh version mesh file
    """

    def __init__(self, filename=None):
        """Initialise Gmsh data structure"""
        self.nodes = {}
        self.elements = {}
        self.filename = filename
        if self.filename:
            self.read()

    def reset(self):
        """Reinitialise Gmsh data structure"""
        self.nodes = {}
        self.elements = {}

    def read(self, mshfile=None):
        """Read a Gmsh .msh file.
        
        Reads Gmsh format 1.0 and 2.0 mesh files, storing the nodes and
        elements in the appropriate dicts.
        """

        if not mshfile:
            mshfile = open(self.filename,'r')

        readmode = 0
        print('Reading %s'%mshfile.name)
        line='a'
        while line:
            line=mshfile.readline()
            line = line.strip()
            if line.startswith('$'):
                if line == '$Nodes':
                    readmode = 1
                elif line == '$Elements':
                    readmode = 3
                elif line == '$MeshFormat':
                    readmode = 4
                else:
                    readmode = 0
            elif readmode:
                columns = line.split()
                if readmode == 4:
                    if len(columns)==3:
                        vno,ftype,dsize=(float(columns[0]),
                                         int(columns[1]),
                                         int(columns[2]))
                        print(('ASCII','Binary')[ftype]+' format')
                    else:
                        endian=struct.unpack('i',columns[0])
                if readmode == 1:
                    # Version 1.0 or 2.0 Nodes
                    try:
                        if ftype==0 and len(columns)==4:
                            self.nodes[int(columns[0])] = map(float,
							                                  columns[1:])
                        elif ftype==1:
                            nnods=int(columns[0])
                            for N in range(nnods):
                                data=mshfile.read(4+3*dsize)
                                i,x,y,z=struct.unpack('=i3d',data)
                                self.nodes[i]=(x,y,z)
                            mshfile.read(1)
                    except ValueError:
                        print('Node format error: '+line, ERROR)
                        readmode = 0
                elif ftype==0 and  readmode > 1 and len(columns) > 5:
                    # Version 1.0 or 2.0 Elements 
                    try:
                        columns = map(int, columns)
                    except ValueError:
                        print('Element format error: '+line, ERROR)
                        readmode = 0
                    else:
                        (id, type) = columns[0:2]
                        if readmode == 2:
                            # Version 1.0 Elements
                            tags = columns[2:4]
                            nodes = columns[5:]
                        else:
                            # Version 2.0 Elements
                            ntags = columns[2]
                            tags = columns[3:3+ntags]
                            nodes = columns[3+ntags:]
                        self.elements[id] = (type, tags, nodes)
                elif readmode == 3 and ftype==1:
                    tdict={1:2,2:3,3:4,4:4,5:5,6:6,7:5,8:3,9:6,10:9,11:10}
                    try:
                        neles=int(columns[0])
                        k=0
                        while k<neles:
                            etype,ntype,ntags=struct.unpack('=3i',
								                            mshfile.read(3*4))
                            k+=ntype
                            for j in range(ntype):
                                mysize=1+ntags+tdict[etype]
                                data=struct.unpack('=%di'%mysize,
                                                   mshfile.read(4*mysize))
                                self.elements[data[0]]=(etype,
                                                        data[1:1+ntags],
                                                        data[1+ntags:])
                    except:
                        raise
                    mshfile.read(1)
                            
        print('  %d Nodes'%len(self.nodes))
        print('  %d Elements'%len(self.elements))

        mshfile.close()
