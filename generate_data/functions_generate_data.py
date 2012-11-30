def initialize_folder():
    import directory_parameters, os
    print "Initializing folder"
    for fname in ['bin','MESH','true_model_seismograms',directory_parameters.in_data_files_path,directory_parameters.OUTPUT_FILES_path,directory_parameters.DATABASES_MPI_path]:
        print "Removing " + fname
        os.system('rm -rf ' + fname)

def link_bin_MESH():
    import os
    print  "Linking binaries"
    os.system('mkdir bin')
    os.chdir('bin/')
    os.system('ln -s ../../binary_files/* .')
    os.chdir('../')
    print "Linking MESH"
    os.system('ln -s ../meshes/MESH .')

def setup_input_data_files():
    import directory_parameters, os
    print "Setting up input data files"
    os.system('mkdir ' + directory_parameters.in_data_files_path)
    os.chdir(directory_parameters.in_data_files_path)
    os.system('cp ../../input_files/Par_file .')
    os.system('ln -s ../../input_files/STATIONS .')
    os.system('ln -s ../../input_files/set_Par_file_option.py .')
    # Set options in Par_file
    os.system('python set_Par_file_option.py MODEL default')
    os.system('python set_Par_file_option.py SAVE_MESH_FILES false')
    print "Linking true-model tomgraphy (.xyz) file"
    os.system('ln -s ../../models/tomography_model_true.xyz tomography_model.xyz')
    # unimportant for now, but required for running xgenerate_databases
    os.system('ln -s ../../input_files/CMTSOLUTION_files/event0 CMTSOLUTION')
    os.chdir('../')

def link_Database():
    import directory_parameters, os
    print "Linking mesh partitioning *Database files"
    os.system('mkdir -p ' + directory_parameters.DATABASES_MPI_path)
    os.system('mkdir -p ' + directory_parameters.OUTPUT_FILES_path)
    root_dir=os.getcwd()+'/'
    os.chdir(directory_parameters.DATABASES_MPI_path)
    os.system('ln -s ' + root_dir+'../models/true_model/*Database .')
    os.chdir(root_dir)

def generate_external_mesh_bin(NPROC):
    import os
    print "Generating true-model databases"
    os.chdir('bin/')
    os.system('mpirun.openmpi -np ' + str(NPROC) + ' xgenerate_databases')
    os.chdir('../')

def simulate_event(enum,NPROC,output_data_directory):
    import os, directory_parameters
    event_dir = 'temp_event'+str(enum)
    root_dir = os.getcwd()+'/'
    os.system('rm -rf ' + event_dir)
    os.system('mkdir ' + event_dir)
    os.chdir(event_dir)
    os.system('mkdir -p ' + directory_parameters.OUTPUT_FILES_path)
    # Link files from root directory
    for dname in ['bin',directory_parameters.DATABASES_MPI_path,directory_parameters.OUTPUT_FILES_path,directory_parameters.in_data_files_path]:
        os.system('mkdir -p ' + dname)
        os.system('ln -s ' + root_dir + dname + '/* ' + dname + '/')

    # Copy CMTSOLUTION file
    os.system('rm ' + directory_parameters.in_data_files_path + 'CMTSOLUTION')
    os.system('cp ../../input_files/CMTSOLUTION_files/event'+str(enum)+' '+directory_parameters.in_data_files_path+'CMTSOLUTION')
    os.chdir('bin')
    os.system('mpirun.openmpi -np ' + str(NPROC) + ' xspecfem3D')
    os.chdir('../')
    dtemp = root_dir + output_data_directory + '/event'+str(enum)+'/'
    os.system('mkdir -p ' + dtemp)
    os.system('mv ' + directory_parameters.OUTPUT_FILES_path+'*XZ.sema ' + dtemp)
    os.chdir(root_dir)
    os.system('rm -rf ' + event_dir)
    print "Event " + str(enum) + " -- done"
