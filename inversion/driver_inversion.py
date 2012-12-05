import os, sys
import functions_inversion
from glob import glob

print "Linking input directory parameters"
os.system('ln -s ../input_files/directory_parameters.py . ')
import directory_parameters

functions_inversion.initialize_folder()
functions_inversion.link_setup()
functions_inversion.copy_initial_model()
functions_inversion.setup_in_data_files()
functions_inversion.setup_event_scratch_directory()
functions_inversion.setup_SEM_list()
functions_inversion.create_STATIONS_ADJOINT()

sys.exit()
for iteration in range(10):
    current_model_directory = os.getcwd()+'/models_inversion/m'+str(iteration)+'/'

    # Generate proc*external_mesh.bin files for current model
    # We This assumes that only the vp-model is changing (vs, rho are fixed at initial model)
    # Par_file has been modified so that input model is read from vp/rho/vs binary files ("MODEL = gll")
    functions_inversion.generate_external_mesh_dot_bin(current_model_directory)
    if iteration==0: functions_inversion.write_grid_gll() # schiemenz, not elegant
    
    print "Doing forward simulations"
    save_forward=True
    simtype='forward'
    functions_inversion.simulation_loop(iteration,simtype,save_forward)

    print "Computing misfit"
    functions_inversion.compute_misfit_write_SEM(iteration)

    print "Doing adjoint simulations"
    save_forward=False
    simtype='adjoint'
    functions_inversion.simulation_loop(iteration,simtype,save_forward)

    print "Summing event kernels"
    functions_inversion.sum_smooth_kernels(iteration)

    print "Updating vp model"
    next_model_directory = os.getcwd()+'/models_inversion/m'+str(iteration+1)+'/'
    functions_inversion.update_model(iteration,current_model_directory,next_model_directory,'alpha_acoustic_kernel','vp')

