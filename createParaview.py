# -*- coding: utf-8 -*-
import numpy as np


import time

def createVtk(u_sol, nodes, triangles, alpha, k):
	"""Create paraview file format UnstructuredGrid"""
	fichier="visualisation_alpha="+str(round(alpha,2))+"_k="+str(round(k,2))+".vtu"
	file = open(fichier,"w")
	file.write('<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n')
	file.write("<UnstructuredGrid>\n")
	file.write('<Piece NumberOfPoints="' + str(len(nodes))+'" NumberOfCells="'+str(len(triangles))+'">\n')
	file.write("<Points>\n")
	file.write('<DataArray NumberOfComponents="3" type="Float64">\n')
	for i in range(1,len(nodes)+1):
		for j in nodes[i]:
			file.write(str(j) + ' ')
		file.write('\n')

	file.write('</DataArray>\n</Points>\n<Cells>\n<DataArray type="Int32" Name="connectivity">\n')

	for i in triangles.values():
		for j in i[2]:
			file.write(str(j-1) + ' ')
		file.write('\n')

	file.write('</DataArray>\n<DataArray type="Int32" Name="offsets">\n')

	for i in range(len(triangles)):
		file.write(str((i+1)*3) + '\n')

	file.write('</DataArray>\n<DataArray type="UInt8" Name="types">\n')

	for i in range(len(triangles)):
		file.write('5\n')

	file.write('</DataArray>\n</Cells>\n<PointData Scalars="solution">\n<DataArray type="Float64" Name="Real part" format="ascii">\n')
	for i in u_sol:
		file.write(str(np.real(i)) + "\n")
	file.write('</DataArray>\n')
	file.write('<DataArray type="Float64" Name="Imag part" format="ascii">\n')
	for i in u_sol:
		file.write(str(np.imag(i)) + "\n")
	file.write('</DataArray>\n')

	file.write('</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n')
	file.close()