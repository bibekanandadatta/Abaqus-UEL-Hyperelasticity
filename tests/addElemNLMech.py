# this python file changes the definition of user element and 
# add the dummy element set based overlaid on the original 
# elements at the end of the original input file. 
# this code is not general-purpose and right now only 
# adds fully-integrated elements. But it can be easily 
# modified to add reduced integration elements. 
# users also have to modify the input file manually in case
# they want to add definition of surfaces based on dummy elements
# to use this script, generate the mesh input file without boundary
# conditions and step or any other additional features.
# run the file from terminal using "python addElemMech.py"
# enter the input file name with .inp extension and enter the element type

import sys

input_fname = input("ABAQUS input file: ")
jtype       = int(input("Enter element type: "))

offset      = 100000

if jtype == 1:
    oldElemStr  = '*Element, type=C3D4'
    newElemStr  = '*User Element,Type=U1,Nodes=4,Coordinates=3,Properties=3,Iproperties=3\n'\
                  '1,2,3\n'\
                  '*Element, type=U1'
elif jtype == 2:
    oldElemStr  = '*Element, type=C3D10'
    newElemStr  = '*User Element,Type=U2,Nodes=10,Coordinates=3,Properties=3,Iproperties=3\n'\
                  '1,2,3\n'\
                  '*Element, type=U2'
elif jtype == 3:
    oldElemStr  = '*Element, type=C3D8'
    newElemStr  = '*User Element,Type=U3,Nodes=8,Coordinates=3,Properties=3,Iproperties=3\n'\
                  '1,2,3\n'\
                  '*Element, type=U3'
elif jtype == 4:
    oldElemStr  = '*Element, type=C3D20'
    newElemStr  = '*User Element,Type=U4,Nodes=8,Coordinates=3,Properties=3,Iproperties=3\n'\
                  '1,2,3\n'\
                  '*Element, type=U4'
elif jtype == 5:
    oldElemStr  = '*Element, type=CPE3'
    newElemStr  = '*User Element,Type=U5,Nodes=3,Coordinates=2,Properties=3,Iproperties=3\n'\
                  '1,2\n'\
                  '*Element, type=U5'
elif jtype == 6:
    oldElemStr  = '*Element, type=CPE6'
    newElemStr  = '*User Element,Type=U6,Nodes=3,Coordinates=2,Properties=3,Iproperties=3\n'\
                  '1,2\n'\
                  '*Element, type=U6'
elif jtype == 7:
    oldElemStr  = '*Element, type=CPE4'
    newElemStr  = '*User Element,Type=U7,Nodes=4,Coordinates=2,Properties=3,Iproperties=3\n'\
                  '1,2\n'\
                  '*Element, type=U7'
elif jtype == 8:
    oldElemStr  = '*Element, type=CPE8'
    newElemStr  = '*User Element,Type=U8,Nodes=4,Coordinates=2,Properties=3,Iproperties=3\n'\
                  '1,2\n'\
                  '*Element, type=U8'
else:
    print("element is not supported")
    sys.exit(1)
    

elem_info   = list()
elem_record = False

#replacing the original element type definition with new definition
with open(input_fname, "r") as fin:
    new_file_content = ""
    for line in fin:
        stripped_line = line.strip()
        new_line = stripped_line.replace(oldElemStr,newElemStr)
        new_file_content += new_line +"\n"
    fin.close()

with open(input_fname, "w") as fout:
    fout.write(new_file_content)
    fout.close()

#reading original element connectivity for dummy elements
with open(input_fname, 'r') as fin:
    lines = fin.readlines()
    for line in lines:
        line_splt = line.strip().split(',')
        
        if len(line_splt) and line_splt[0] == '*Element':
            elem_record = True
            continue
        
        if len(line_splt) and line_splt[0] == '*System':
            elem_record = False
        
        if elem_record is True:
            elem_list = [int(x) for x in line_splt]
            elem_list[0] += offset
            elem_info.append(elem_list)
    fin.close()

# open the input file again for modification
with open(input_fname,'a+') as fout:

    #add dummy elements to the end of the file
    if jtype == 1:
        fout.write('*Element, type=C3D4')
    elif jtype == 2:
        fout.write('*Element, type=C3D10')
    elif jtype == 3:
        fout.write('*Element, type=C3D8')
    elif jtype == 4:
        fout.write('*Element, type=C3D20')
    elif jtype == 5:
        fout.write('*Element, type=CPE4')
    elif jtype == 6:
        fout.write('*Element, type=CPE6')
    elif jtype == 7:
        fout.write('*Element, type=CPE4')
    elif jtype == 8:
        fout.write('*Element, type=CPE8')
    else:
        print("element is not supported")
        sys.exit(1)

    for elem_list in elem_info:
        fout.write('\n')
        fout.write(str(elem_list)[1:-1])
    fout.write('\n*Elset, elset=ElDummy, generate\n')
    fout.write(str(elem_info[0][0]) + "," + str(elem_info[-1][0]) + "," + str(1))

    fout.close()
    
print("dummy elements added to the input file")