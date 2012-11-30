def initialize_folder():
    import os, directory_parameters
    print "Initializing folder"
    for fname in ['bin','MESH','true_model','initial_model','tomography_model_true.xyz','tomography_model_starting.xyz',directory_parameters.in_data_files_path,directory_parameters.OUTPUT_FILES_path,directory_parameters.DATABASES_MPI_path]:
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
    import os, directory_parameters
    print "Creating xyz files"
    execfile('create_xyz_starting.py')
    execfile('create_xyz_true.py')
    os.system('mkdir ' + directory_parameters.in_data_files_path)
    os.chdir(directory_parameters.in_data_files_path)
    os.system('cp ../../input_files/Par_file .')
    for fname in ['STATIONS','set_Par_file_option.py']:
        os.system('ln -s ../../input_files/'+fname+' .')
    # unimportant for now, but required for running xgenerate_databases
    os.system('ln -s ../../input_files/CMTSOLUTION_files/event0 CMTSOLUTION')
    # set model type in Par_file
    os.system('python set_Par_file_option.py MODEL default')
    # must have SAVE_MESH_FILES = true to generate model (vp/vs/rho) binaries
    os.system('python set_Par_file_option.py SAVE_MESH_FILES true')
    os.chdir('../')

def partition_mesh(num_processors):
    import os, directory_parameters
    print "Running mesh partitioning"
    os.system('mkdir -p ' + directory_parameters.DATABASES_MPI_path)
    os.system('mkdir -p ' + directory_parameters.OUTPUT_FILES_path)
    os.system('./bin/xdecompose_mesh_SCOTCH ' + str(num_processors) + ' MESH/ ' + directory_parameters.DATABASES_MPI_path)

def generate_tomographic_model(xyz_filename,num_processors):
    import os, directory_parameters
    print "Generating model from tomo model: " + xyz_filename
    os.system('rm '+directory_parameters.in_data_files_path+'tomography_model.xyz')
    os.system('cp ' + xyz_filename + ' ' +directory_parameters.in_data_files_path+'tomography_model.xyz')
    os.chdir('bin/')
    os.system('mpirun.openmpi -np ' + str(num_processors) + ' xgenerate_databases')
    os.chdir('../')

def move_model_files(output_dir):
    import os, directory_parameters
    print "Moving model files to directory " + output_dir
    os.system('mkdir ' + output_dir)
    for ftype in ['vp','rho','vs']:
        os.system('mv ' + directory_parameters.DATABASES_MPI_path+'*'+ftype+'.bin ' + output_dir)

def copy_Database_files(output_dir):
    import os, directory_parameters
    print "Copying Database files from " + directory_parameters.DATABASES_MPI_path + " to " + output_dir
    os.system('cp ' + directory_parameters.DATABASES_MPI_path+'*Database ' + output_dir)

