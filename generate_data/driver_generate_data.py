from glob import glob
import functions_generate_data
import os, pprocess

NPROC = 4 # number of processors for 1 specfem simulation
total_number_processors = 36 # total number of processors available
num_events = len(glob('../input_files/CMTSOLUTION_files/*'))
output_data_directory = 'true_model_seismograms'

print "Linking python files"
os.system('ln -s ../input_files/directory_parameters.py . ')
import directory_parameters

functions_generate_data.initialize_folder()
functions_generate_data.link_bin_MESH()
functions_generate_data.setup_input_data_files()
functions_generate_data.link_Database()
functions_generate_data.generate_external_mesh_bin(NPROC)

print "Simulating data for " + str(num_events) + " events"
# Creates a temporary directory for each data simulation
# This way, events can be simulated in parallel with no problem
results = pprocess.Map(limit=total_number_processors/NPROC,reuse=1)
parallel_function = results.manage(pprocess.MakeReusable(functions_generate_data.simulate_event))
for i in range(num_events):
    parallel_function(i,NPROC,output_data_directory)
results.finish()
