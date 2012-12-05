from glob import glob

NPROCS=4 # number of processors used in 1 SPECFEM3D solve
total_processors = 36 # total number of processors available 

number_of_events=len(glob('../input_files/CMTSOLUTION_files/*'))
zero_line = -1.0 # kernels are zeroed-out for z > zero_line (i.e. no model update there)

# Where the exact-model data is stored
data_path = '/export/data/schiemenz/workflow_examples/toy_cube/generate_data/true_model_seismograms/'

# Where source files ('CMTSOLUTION') for each event are stored
CMTSOLUTION_path = '/export/data/schiemenz/workflow_examples/toy_cube/input_files/CMTSOLUTION_files/'

# Where the individual events (both forward and adjoint) are simulated
# Note that most of the disk space for the inversion is consumed here
event_scratch_directory='/export/data/schiemenz/workflow_examples/toy_cube/inversion/event_scratch_directory/'

# Currently this workflow is specialized only for acoustic-model, VP inversion
kernels_to_save = ['alpha_acoustic_kernel']

# path where raw kernels computed from adjoint solve are saved
unsmoothed_kernel_path=event_scratch_directory+'kernels_unsmoothed/'

summed_kernel_directory = 'summed_kernels/'

step_length=0.025

# Decide whether to save individual traces for each iteration and simulation
# These are placed in <scratch directory event path>/save_traces/it<iteration number>
save_all_synthetics = False

# Decide whether to save summed event kernels each iteration
# These are placed in ./kernelsave/it<iteration number>
save_summed_kernels = False

