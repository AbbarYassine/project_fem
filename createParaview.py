# -*- coding: utf-8 -*-
import numpy as np

from paraview.simple import *
import time



def pngImage(filename):
	reader = XMLUnstructuredGridReader(FileName=filename)


	Show(reader)
	Render()
	ai = reader.PointData[1]
	dp = GetDisplayProperties(reader)
	print(ai.GetRange())
	# To color the representation by an array, we need to first create
	# a lookup table.  We use the range of the Elevation array
	dp.LookupTable = MakeBlueToRedLT(ai.GetRange()[0], ai.GetRange()[1])
	dp.ColorArrayName = ("POINTS", "Real part")
	Render()

	#save screenshot
	print(filename)
	WriteImage(filename+"_img.png")

def createVtkfile(u_sol, nodes, elements,alpha): 
	nomfichier="OutputVtu/out_alpha="+str(alpha)+".vtu"
	file = open(nomfichier,"w")
	file.write('<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n')
	file.write("<UnstructuredGrid>\n")
	file.write('<Piece NumberOfPoints="' + str(len(nodes))+'" NumberOfCells="'+str(len(elements))+'">\n')
	file.write("<Points>\n")
	file.write('<DataArray NumberOfComponents="3" type="Float64">\n')
	for i in nodes:
		for j in i:
			file.write(str(j) + ' ')
		file.write('\n')

	file.write('</DataArray>\n</Points>\n<Cells>\n<DataArray type="Int32" Name="connectivity">\n')

	for i in elements:
		for j in i:
			file.write(str(j) + ' ')
		file.write('\n')

	file.write('</DataArray>\n<DataArray type="Int32" Name="offsets">\n')

	for i in range(len(elements)):
		file.write(str((i+1)*3) + '\n')

	file.write('</DataArray>\n<DataArray type="UInt8" Name="types">\n')

	for i in range(len(elements)):
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
	pngImage(nomfichier)