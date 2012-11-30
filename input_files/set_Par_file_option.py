import sys

def check_MODEL(value):
    # sys.argv[1:][0] is the model option
    # available options are:
    #   default (model parameters described by mesh properties)
    # 1D models available are:
    #   1d_prem,1d_socal,1d_cascadia
    # 3D models available are:
    #   aniso,external,gll,salton_trough,tomo
    allowed_MODEL_types = ['default','1d_prem','1d_socal','1d_cascadia','aniso','external','gll','salton_trough','tomo']
    if not (value in allowed_MODEL_types):
        print "ERROR:"
        print "Requested MODEL type is " + value + ", but only supported values are " + str(allowed_MODEL_types)
        sys.exit(1)

def check_SAVE_MESH_FILES(value):   
    allowed_inputs = ['.true.','.false.']
    if not (value.lower() in allowed_inputs):
        print "ERROR:"
        print "Requested value is " + value + ", but only supported values are " + str(allowed_inputs)
        sys.exit(1)

# assumes all parameters are in upper-case
input_param = sys.argv[1:][0].upper() 
param_value = sys.argv[1:][1]
# assumes all (non-numeric) parameter values are lower-case
if type(param_value)==str:
    param_value = param_value.lower()
# convert boolean data types to proper format
if param_value in ['true','false']:
    param_value = '.'+param_value+'.'

# Parameters that this script supports so far
allowed_parameters = ['MODEL','SAVE_MESH_FILES']
if input_param not in allowed_parameters:
    print "ERROR: requested input parameter is " + input_param + ", but this script only supports " + str(allowed_parameters)
    sys.exit(1)
else:
    if input_param=='MODEL':
        check_MODEL(param_value)
    elif input_param=='SAVE_MESH_FILES':
        check_SAVE_MESH_FILES(param_value)

print "Setting parameter " + input_param + " to value " + param_value


# Read the Par_file
fin = open('Par_file','r')
fl=fin.readlines()
fin.close()

# Output a new Par_file
fout = open('Par_file','w')
for line in fl:
    if len(line.split())==0:
        fout.write(line)
    elif line.split()[0]==input_param:
        fout.write(input_param+' = ' + param_value +'\n')
    else:
        fout.write(line)
