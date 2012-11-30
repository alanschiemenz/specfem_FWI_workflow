import os, functions_models

print "Setting up input files (linking from ../input_files)"
os.system('ln -s ../input_files/directory_parameters.py . ')
import directory_parameters

NPROC=4
print "Running with " + str(NPROC) + " processors"

functions_models.initialize_folder()
functions_models.link_bin_MESH()
functions_models.setup_input_data_files()
functions_models.partition_mesh(NPROC)
functions_models.generate_tomographic_model('tomography_model_true.xyz',NPROC)
functions_models.move_model_files('true_model/')
functions_models.copy_Database_files('true_model/')
functions_models.generate_tomographic_model('tomography_model_starting.xyz',NPROC)
functions_models.move_model_files('initial_model/')
functions_models.copy_Database_files('initial_model/')
