#! /usr/bin/python
import fileinput, string, sys ,os ,shutil

#####################################
def readPARAM(var_Name,fileQuery):
    flag = 0
    var_Name=var_Name.strip()
    var_Value="###"
    ##################################################################
    if flag == 0:
	    file1 = open(fileQuery, "r") 
    	    text1 = file1.readlines() 
    	    file1.close() 
	    #print var_Name
            for line1 in text1:
		#print var_Name
		word1 = line1.split()
		if len(word1) > 1:
			name1 = word1[0].strip()   
			var_Value=line1
			var_Value=var_Value.replace(name1,"")
			var_Value=var_Value.strip() 
			var_Value=var_Value.strip(";") 
			if var_Name==name1:
				return var_Value
	    return "###"
##################################################################
def updatePARAM(var_Name, var_Value,fileQuery):
        var_Name=var_Name.strip()
        file4 = open(fileQuery, "r")
        text4 = file4.readlines()
        file4.close()
        file4 = open(fileQuery, "w")
        #print var_Name
        for line4 in text4:
                #print var_Name
                word4 = line4.split()
                if len(word4) > 1:
                        name4 = word4[0].strip()
                        if var_Name.find(name4) != -1:
                                if var_Name==name4:
                                        flag = 1
                                        line4 = var_Name+"\t\t"+var_Value+";\n"
                                        #os.system("touch t 0/*")
                                        #print "0"
                file4.write(line4)
        file4.close()
#############################
######################################################################

os.system("rm -rf constant/extendedFeatureEdgeMesh")
os.system("mv constant/polyMesh/blockMeshDict ../../")
os.system("rm -rf constant/polyMesh/*")
os.system("mv blockMeshDict constant/polyMesh/")
os.system("rm -rf constant/triSurface/*")
os.system("rm -f system/snappyHexMeshGeometryRefinement")
os.system("cp -ar system/surfaceFeatureExtract_bak system/surfaceFeatureExtract")

######################################################################

########################### blockMesh ################################
cellsize=readPARAM("cellsize","input")
min_coord_x=readPARAM("min_coord_x","input")
min_coord_y=readPARAM("min_coord_y","input")
min_coord_z=readPARAM("min_coord_z","input")
max_coord_x=readPARAM("max_coord_x","input")
max_coord_y=readPARAM("max_coord_y","input")
max_coord_z=readPARAM("max_coord_z","input")

ncells_x = int((float(max_coord_x)-float(min_coord_x))/int(cellsize))
ncells_y = int((float(max_coord_y)-float(min_coord_y))/int(cellsize))
ncells_z = int((float(max_coord_z)-float(min_coord_z))/int(cellsize))

file6 = open('blocktmp','w')
	file6.write("ncells_x	")
	file6.write(ncells_x)
	file6.write(";")
	file6.write("ncells_y	")
	file6.write(ncells_y)
	file6.write(";")
	file6.write("ncells_z	")
	file6.write(ncells_z)
	file6.write(";")
file6.close()


###################### refinement region ################################
    
stlfilestl=readPARAM("stlfile","input")
stlfile=readPARAM("stlfile","input")
ref_box1_coord_min=readPARAM("ref_box1_coord_min","input")
ref_box1_coord_max=readPARAM("ref_box1_coord_max","input")
ref_box2_coord_min=readPARAM("ref_box1_coord_min","input")
ref_box2_coord_max=readPARAM("ref_box1_coord_max","input")

os.system("cp " +stlfile " constant/triSurface/")

stlfile=stlfile.strip(".stl")

updatePARAM("name",stlfile,"system/snappyHexMeshGeometry")

file1 = open('snappyHexMeshGeometryRefinement','w')
  file1.write('refinementBox1\n{\ntype searchableBox;\nmin (')
  file1.write(ref_box1_coord_min)
  file1.write(');\n')
  file1.write('max (')
  file1.write(ref_box1_coord_max)
  file1.write(');\n}\n')
  file1.write('refinementBox2\n{\ntype searchableBox;\nmin (')
  file1.write(ref_box2_coord_min)
  file1.write(');\n')
  file1.write('max (')
  file1.write(ref_box2_coord_max)
  file1.write(');\n}\n')
file1.close()

ref_box1_level=readPARAM("ref_box1_level","input")
ref_box2_level=readPARAM("ref_box2_level","input")

file2 = open('snappyHexMeshRefinement','w')
  file2.write('refinementBox1\n{\nmode inside;\nlevels ((1E15 ')
  file2.write(ref_box1_level)
  file1.write('));\n}\n')
  file2.write('refinementBox2\n{\nmode inside;\nlevels ((1E15 ')
  file2.write(ref_box2_level)
  file1.write('));\n}\n')

file2.close()
##################### surfaceFeatureExtract #######################
text1 = "test.stl"

#print ("Text to replace it with:",textToReplace)

fileToSearch=open('system/surfaceFeatureExtractDict','r')
text2 = fileToSearch.readlines()

fileToSearch1=open('system/surfaceFeatureExtractDict','w')

for line in text2:
    fileToSearch1.write( line.replace( text1, stlfilestl ) )
fileToSearch1.close()
fileToSearch.close()
#######################################################################

######################## Feature edge refinement ######################

file3 = open('snappyHexMeshFeature','w')
  file3.write("{\nfile \"")
  file3.write(stlfile)
  file3.write(".eMesh\";\nlevel 0;\n}")

######################################################################

##################### refinement surface #############################
file4=open('input','r')
textall = file4.readlines()
file4.close()
file5=open('snappyHexMeshRefSurface','w')
switch = 0
for ref1 in textall:
	if "refinement_surface_end" in ref1:
		switch = 0 
	if switch == 1:
		ref2 = ref1.split()
		ref2[1] =ref2[1].strip(';')
		file5.write(ref2[0])
		file5.write("\n{\n level (")
		file5.write(ref2[1])
		file5.write(" ")
		file5.write(ref2[1])
		file5.write(");\n}\n")
	if "refinement_surface_start" in ref1:
		switch = 1
file5.close()
os.system("mv snappyHexMeshRefSurface system/")

######################################################################

################################# Run Mesh ###########################
No_of_Processors=readPARAM("No_of_Processors","input")

os.system("surfaceFeatureExtract")
os.system("blockMesh")
os.system("decomposePar")
snappy="mpirun -np "+No_of_Processors+" snappyHexMesh -parallel >& log &"
os.system(snappy)
######################################################################







