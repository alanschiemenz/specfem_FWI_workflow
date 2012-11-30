def initialize_folder():
    import directory_parameters, input_parameters, os
    print "Initializing folder"
    for fname in ['bin','MESH','linked_data','linked_CMTSOLUTION_files',directory_parameters.in_data_files_path,directory_parameters.OUTPUT_FILES_path,directory_parameters.DATABASES_MPI_path,input_parameters.event_scratch_directory,'L2_misfit','models_inversion']:
        print "Removing " + fname
        os.system('rm -rf ' + fname)

def link_setup():
    import os, input_parameters, directory_parameters
    print "Linking MESH"
    os.system('ln -s ../meshes/MESH .')
    print "Linking binary executables"
    os.system('mkdir bin')
    os.system('ln -s ' + os.getcwd()+'/../binary_files/* bin/')
    print "Linking data and CMTSOLUTION_files"
    os.system('mkdir linked_CMTSOLUTION_files linked_data')
    for enum in range(input_parameters.number_of_events):
        os.chdir('linked_data')
        os.system('ln -s '+input_parameters.data_path+'event'+str(enum)+' .')
        os.chdir('../')
        os.chdir('linked_CMTSOLUTION_files')
        os.system('ln -s ' + input_parameters.CMTSOLUTION_path + 'event'+str(enum) + ' .')
        os.chdir('../')
    print "Linking initial model files and (constant) vs, rho, Database files"
    os.system('ln -s ../models_inversion/initial_model . ')
    os.system('mkdir -p ' + directory_parameters.OUTPUT_FILES_path)
    os.system('mkdir -p ' + directory_parameters.DATABASES_MPI_path)
    for fname in ['vs.bin','rho.bin','Database']:
        os.system('ln -s ' + os.getcwd()+'/initial_model/*' + fname + ' ' + directory_parameters.DATABASES_MPI_path)

def copy_initial_model():
    import os
    print "Copying initial vp-model to models_inversion/m0"
    os.system('mkdir -p models_inversion/m0')
    os.system('cp initial_model/*vp.bin models_inversion/m0')


def read_bin_file(filename, data_types):
    import array
    num_objects = len(data_types)
    f = open(filename,'rb')
    file_out=list()
    for n in range(num_objects):
        binlength = array.array('i')
        binlength.read(f,1)
        binvalues = array.array(data_types[n])
        if (data_types[n]=='d'):
            binvalues.read(f,binlength[0]/8) # double precision
        else:
            binvalues.read(f,binlength[0]/4) # 4 for single precision
        binlength.read(f,1)
        file_out.append(binvalues)
    f.close()
    return file_out

def write_bin_file(filename, file_in, data_types, num_objects):
    import array
    # file_in is a LIST
    fileobj = open(filename, mode='wb')        
    for n in range(num_objects):
        fheader = array.array('i')
        fheader.fromlist([len(file_in[n])*4]) # 4 bytes per length of input
        outvalues = array.array(data_types[n])
        outvalues.fromlist(file_in[n])#.tolist())
        fheader.tofile(fileobj)
        outvalues.tofile(fileobj)
        fheader.tofile(fileobj)
    fileobj.close()
    

def read_bin_file_noheader(filename, data_type): 
    f = open(filename,'rb')
    file_out=list()
    num_bytes = os.path.getsize(filename)
    binvalues = array.array(data_type)        
    binvalues.read(f,num_bytes/4) # 4 for single precision
    file_out.append(binvalues)
    f.close()
    return file_out

def compute_misfit_write_SEM(iteration,compute_scalar_only=False):
    import input_parameters, os, pprocess, sys
    import numpy as np
    print "Computing adjoint source"
    print "Assuming the default SPECFEM3D format for data: 2 columns of ascii text, < time value, recorded waveform >" 
    
    def process_adj_file(line):#syn_file,data_file,adj_file):
        syn_file = line.split()[0]
        data_file = line.split()[1]
        adj_file = line.split()[2]
        syn_trace = np.loadtxt(syn_file)
        data_trace = np.loadtxt(data_file)

        if len(syn_trace)>len(data_trace):
            print "ERROR: synthetic traces have more time samples thatn data traces!"
            sys.exit(1)
        elif len(syn_trace)<len(data_trace):
            # Assuming that synthetics and data are sampled at equivalent intervals
            # Synthetics can be shorter time series than data.  In this case we just read the first <length of synthetic trace> samples of data
            data_trace=data_trace[:len(syn_trace),:]        
        
        pressure_difference = syn_trace[:,1]-data_trace[:,1]    
        adj_source = np.zeros(len(syn_trace))
        # for acoustic FWI, L2 adjoint source is the second derivative of pressure difference
        # See Peter et al. 2011, GJI, eq. (A8)
        for i in range(1,len(adj_source)-1):
            # do a simple central finite difference
            adj_source[i] = pressure_difference[i+1]-2*pressure_difference[i]+pressure_difference[i-1]
        dt = syn_trace[1,0]-syn_trace[0,0]
        adj_source *= 1.0/(dt**2)
        return_data_misfit = 0.5 * sum(adj_source**2) * dt
        if compute_scalar_only == False:
            ftemp = open(adj_file,'w')
            for i in range(len(adj_source)):
                ftemp.write('%e %e\n' % (syn_trace[i,0],adj_source[i]))
            ftemp.close()
            if iteration==0: # schiemenz, gotta be a better way to do this
                os.system('ln -s ' + adj_file + ' ' + adj_file[:-5] + 'Y.adj')
                os.system('ln -s ' + adj_file + ' ' + adj_file[:-5] + 'Z.adj')                      
        return return_data_misfit

    results_misfit = pprocess.Map(limit=input_parameters.total_processors,reuse=1)
    parfun = results_misfit.manage(pprocess.MakeReusable(process_adj_file))
    f=open('SEM_file_list','r')
    SEM_file_list = f.readlines()
    f.close()
    for line in SEM_file_list:
        parfun(line)
    results_misfit.finish()

    total_misfit = 0.0
    for i in range(len(SEM_file_list)):
        total_misfit += results_misfit[i]
        
    fmisfit=open('L2_misfit','a')
    if iteration==0:
        fmisfit.write("Iteration, L2-misfit\n")
    fmisfit.write('%i, %e\n' % (iteration,total_misfit))
    fmisfit.close()
 
def setup_in_data_files():
    import os, directory_parameters
    print "Setting up input data files in directory " + directory_parameters.in_data_files_path
    print "By default, this is copying Par_file and STATIONS from directory ../input_files/"
    root_dir = os.getcwd()+'/'
    os.system('mkdir -p ' + directory_parameters.in_data_files_path)
    os.system('ln -s ' + root_dir + '../input_files/STATIONS ' +  directory_parameters.in_data_files_path)
    os.system('ln -s ' + root_dir + '../input_files/set_Par_file_option.py ' +  directory_parameters.in_data_files_path)
    os.system('cp ../input_files/Par_file ' + directory_parameters.in_data_files_path)
    os.chdir(directory_parameters.in_data_files_path)
    # set options in Par_file
    os.system('python set_Par_file_option.py MODEL gll')
    os.system('python set_Par_file_option.py SAVE_MESH_FILES false')
    os.chdir(root_dir)


def generate_external_mesh_dot_bin(vp_model_dir): 
    import os, directory_parameters, input_parameters
    print "Linking vp model files from " + vp_model_dir + " to " + directory_parameters.DATABASES_MPI_path
    os.system('rm ' + directory_parameters.DATABASES_MPI_path + '/*vp.bin')
    os.system('ln -s ' + vp_model_dir + '/*vp.bin ' + directory_parameters.DATABASES_MPI_path)
    # xgenerate_databases requires a CMTSOLUTION file to be placed in the input data folder
    os.system('cp linked_CMTSOLUTION_files/event0 ' +  directory_parameters.in_data_files_path + 'CMTSOLUTION')
    # With our current MESH/nummateterial_velocity_file settings, a tomography_model.xyz file is needed for xgenerate_databasese
    # However, this file is not actually read when Par_file -> MODEL type is set to gll
    os.system('cp ../models/tomography_model_starting.xyz ' + directory_parameters.in_data_files_path + 'tomography_model.xyz')
    print "Generating *external_mesh.bin databases"
    os.chdir('bin/')
    os.system('mpirun.openmpi -np ' + str(input_parameters.NPROCS) + ' xgenerate_databases')
    os.chdir('../')

def setup_event_scratch_directory():
    import directory_parameters, input_parameters, os
    root_dir=os.getcwd()+'/'
    for enum in range(input_parameters.number_of_events):
        event_dir = input_parameters.event_scratch_directory+'event'+str(enum)+'/'
        os.system('mkdir -p ' + event_dir)
        os.system('mkdir -p ' + event_dir + '/bin')
        os.system('mkdir -p ' + event_dir + directory_parameters.in_data_files_path)
        os.system('mkdir -p ' + event_dir + directory_parameters.OUTPUT_FILES_path)
        os.system('mkdir -p ' + event_dir + directory_parameters.DATABASES_MPI_path)
        os.system('mkdir -p ' + event_dir + directory_parameters.SEM_path)
        os.system('mkdir -p ' + event_dir + 'traces_forward')
        os.system('ln -s ' + root_dir + '/bin/* ' + event_dir + '/bin/')
        os.system('ln -s ' + root_dir + directory_parameters.in_data_files_path + 'Par_file ' + event_dir + directory_parameters.in_data_files_path)
        os.system('ln -s ' + root_dir + directory_parameters.in_data_files_path + 'STATIONS ' + event_dir + directory_parameters.in_data_files_path)
        os.system('ln -s ' + root_dir + 'linked_CMTSOLUTION_files/event' + str(enum) + ' ' + event_dir + directory_parameters.in_data_files_path + 'CMTSOLUTION')
        os.system('ln -s ' + root_dir + directory_parameters.OUTPUT_FILES_path+'values_from_mesher.h ' + event_dir + directory_parameters.OUTPUT_FILES_path)

        # Link *external_databases.bin into event directories, from root inversion directory
        # At first these are not defined, but will be once the inversion loop gets going
        for proc in range(input_parameters.NPROCS):
            fn = get_proc_prefix(proc) # returns string "proc<6-digit prefix>_"
            os.system('ln -s ' + root_dir + directory_parameters.DATABASES_MPI_path + fn + 'external_mesh.bin ' + event_dir + directory_parameters.DATABASES_MPI_path)

def get_proc_prefix(p):
    # Most SPECFEM3D filenames follow the convention "proc<6-digit prefix>_<name>.bin"
    # This extracts the prefix
    if p<10: return 'proc00000'+str(p)+'_'
    elif p<100: return 'proc0000'+str(p)+'_'
    elif p<1000: return 'proc000'+str(p)+'_'
    elif p<10000: return 'proc00'+str(p)+'_'
    elif p<100000: return 'proc0'+str(p)+'_'
    else: return 'proc'+str(p)+'_'
    
def simulation_loop(iteration,simtype='forward',save_forward=True):
    import input_parameters, os, pprocess, sys, directory_parameters
    def do_event(event):
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print "     (ITERATION, EVENT, SIMULATION TYPE)=  : (" + str(iteration) + ' '+str(event)+' ' + simtype+')'
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" 
        inversion_home_directory=os.getcwd()
        event_dir = input_parameters.event_scratch_directory+'event'+str(event)+'/'
        os.chdir(event_dir+'bin')
        os.system('mpirun.openmpi -np ' + str(input_parameters.NPROCS) + ' xspecfem3D')        
        os.chdir(inversion_home_directory) # back to root directory
        if simtype=='forward':
            if input_parameters.save_all_synthetics == True: # save traces and link them into directory to be read for adjoint source construction
                save_trace_dir = event_dir + '/save_traces/it' + str(iteration) + '/'
                os.system('mkdir -p ' + save_trace_dir)                
                os.system('mv ' + event_dir + directory_parameters.OUTPUT_FILES_path + '*XZ.sema ' + save_trace_dir)
                # schiemenz: gotta be a better way to do this
                os.system('rm ' + event_dir + 'traces_forward/*')
                os.system('ln -s ' + save_trace_dir + '*' + ' ' + event_dir + 'traces_forward/')
            else: # just move the traces into a directory to be read for adjoint source construction
                os.system('mv ' + event_dir + directory_parameters.OUTPUT_FILES_path + '*XZ.sema ' + event_dir + 'traces_forward/')
        elif simtype=='adjoint':
            for kerneltype in input_parameters.kernels_to_save:
                kernel_dir = input_parameters.unsmoothed_kernel_path+'/it'+str(iteration)+'/event'+str(event)+'/'
                os.system('mkdir -p ' + kernel_dir)
                os.system('mv ' + event_dir + directory_parameters.DATABASES_MPI_path + '*'+ kerneltype + '.bin ' + kernel_dir)
        
    print "Changing Par_file for simtype = " + simtype + ', save_forward = ' + str(save_forward)
    if simtype not in ['forward','adjoint']:
        print 'ERROR: simulation type must be forward or adjoint'
        sys.exit(1)
    elif simtype=='forward':
        if save_forward==True: os.system('./change_simulation_type.pl -F')
        else: os.system('./change_simulation_type.pl -f')
    elif simtype=='adjoint':
        os.system('./change_simulation_type.pl -b')            

    # Simulate events, an embarrasingly-parallel loop
    results = pprocess.Map(limit=input_parameters.total_processors/input_parameters.NPROCS,reuse=1)
    parfun = results.manage(pprocess.MakeReusable(do_event))
    for i in range(input_parameters.number_of_events):
        parfun(i)
    results.finish()
        
def setup_SEM_list():
    # Constructs a file with 3 columns: <synthetic file name, data file name, adjoint source file name>
    # This allows for easier implementation later when computing adjoint sources
    
    from glob import glob
    import input_parameters, os, directory_parameters    
    ftemp=open('SEM_file_list','w')
    for eventnum in range(input_parameters.number_of_events):
        event_dir=input_parameters.event_scratch_directory+'event'+str(eventnum)+'/'
        data_fnames = glob('linked_data/event'+str(eventnum)+'/*XZ.sema')
        for f in data_fnames:
            filename = f[f.rfind('/')+1:] # gives file name without path
            synthetic_name = event_dir+'/traces_forward/'+filename
            # replace "*XZ.sema" with "*XX.adj", per SPECFEM3D conventions
            adj_name = event_dir+directory_parameters.SEM_path+filename[:-7]+'XX.adj'
            ftemp.write("%s %s %s\n" % (synthetic_name, f, adj_name))
    ftemp.close()
                             
def create_STATIONS_ADJOINT():
    import numpy as np
    import directory_parameters, input_parameters, os
    root_dir = os.getcwd()
    for eventnum in range(input_parameters.number_of_events):
        os.chdir(input_parameters.event_scratch_directory+'event'+str(eventnum))
        statfile = open(directory_parameters.in_data_files_path+'/STATIONS',"r")
        slines = statfile.readlines()
        statfile.close()
        SAfile = open(directory_parameters.in_data_files_path+'/STATIONS_ADJOINT','w')
        for i in range(len(slines)):
            yr = float(slines[i].split()[2])
            xr = float(slines[i].split()[3])
            SAfile.write("%s" % (slines[i]))  
        SAfile.close()
        os.chdir(root_dir)


def sum_smooth_kernels(iteration):
    import os, input_parameters, pprocess, directory_parameters
    from glob import glob
    import numpy as np

    print "Summing kernels: " + str(input_parameters.kernels_to_save)
    print "These are placed in " + input_parameters.summed_kernel_directory
    os.system('rm -rf '+input_parameters.summed_kernel_directory)
    os.system('mkdir -p '+input_parameters.summed_kernel_directory)
    
    def sum_over_kernels(proc):
        proc_prefix = get_proc_prefix(proc) # returns string "proc<6-digit prefix>_"
        for kerneltype in input_parameters.kernels_to_save:
            for event in range(input_parameters.number_of_events):
                # Read file proc<proc index>_<kernel type>.bin, e.g. proc000000_alpha_acoustic_kernel.bin
                fname_kernel = input_parameters.unsmoothed_kernel_path+'/it'+str(iteration)+'/event'+str(event)+'/'+proc_prefix+kerneltype+'.bin'
                fileout = read_bin_file(fname_kernel,['f'])[0] # read first item from list-output
                data_single = np.array(fileout,dtype='f')
                if event==0:
                    data_summed = data_single
                else:
                    data_summed += data_single
            filename_out = input_parameters.summed_kernel_directory+'/'+proc_prefix+kerneltype+'.bin'
            write_bin_file(filename_out,[data_summed.tolist()],['f'],1)
            print "Finished summation over events: " + proc_prefix+kerneltype+'.bin'

    results = pprocess.Map(limit=input_parameters.total_processors,reuse=1)
    parfun = results.manage(pprocess.MakeReusable(sum_over_kernels))
    for i in range(input_parameters.NPROCS):
        parfun(i)
    results.finish()

    def zero_out_kernels(proc):
        proc_prefix = get_proc_prefix(proc) # returns string "proc<6-digit prefix>_"
        for kerneltype in input_parameters.kernels_to_save:
            kernel_fname = input_parameters.summed_kernel_directory+'/'+proc_prefix+kerneltype+'.bin'            
            kernelvalues = read_bin_file(kernel_fname,['f'])[0]
            kernelvalues = np.array(kernelvalues,dtype='f')
            zlocal = read_bin_file(directory_parameters.DATABASES_MPI_path+proc_prefix+'zlocal.bin',['d'])[0]
            for gp in range(len(zlocal)):
                if (zlocal[gp] > input_parameters.zero_line):
                    kernelvalues[gp] = 0.0
            write_bin_file(kernel_fname,[kernelvalues.tolist()],['f'],1)
    print "Zeroing out summed kernels for z > " + str(input_parameters.zero_line)
    results = pprocess.Map(limit=input_parameters.NPROCS,reuse=1)
    calc = results.manage(pprocess.MakeReusable(zero_out_kernels))
    for proc in range(input_parameters.NPROCS):
        calc(proc)
    results.finish()

    # Optional feature: save kernels
    if input_parameters.save_summed_kernels==True:
        def savekernels(procnum):
            proc_prefix = get_proc_prefix(procnum) # returns string "proc<6-digit prefix>_"
            tmpfile = input_parameters.summed_kernel_directory+'/'+proc_prefix+'alpha_acoustic_kernel'+'.bin'
            os.system('cp '+tmpfile+' '+ksave_dir)
        ksave_dir = 'kernelsave/it'+str(iteration)
        os.system("mkdir -p "+ksave_dir)
        results = pprocess.Map(limit=input_parameters.NPROCS,reuse=1)
        parallel_function = results.manage(pprocess.MakeReusable(savekernels))
        for ppp in range(input_parameters.NPROCS):
            parallel_function(ppp)
        results.finish() # wait for things to finish before moving on


def update_model(iteration,model_in_dir,model_out_dir,kerneltype,modeltype):#,update_files,step_length):
    # e.g. kerneltype = 'alpha_acoustic_kernel'
    # e.g. modeltype = 'vp'
    import input_parameters, os, pprocess, sys
    import numpy as np
    
    def update_proc(current_proc,maxkernel):
        proc_prefix = get_proc_prefix(current_proc) # returns string "proc<6-digit prefix>_"        
        modelvalues = read_bin_file(model_in_dir+proc_prefix+modeltype+'.bin',['f'])[0]  # current model values
        modelvalues = np.array(modelvalues,dtype='f')
        search_direction = read_bin_file(input_parameters.summed_kernel_directory+proc_prefix+kerneltype+'.bin',['f'])[0]  
        search_direction = np.array(search_direction,dtype='f')
        # use steepest descent algorithm with fixed step length
        # normalize update direction by maximum absolute gradient
        newmodel = modelvalues * (1 + input_parameters.step_length * search_direction / maxkernel)
        file_out = model_out_dir+'/'+proc_prefix+modeltype+'.bin'
        write_bin_file(file_out,[newmodel.tolist()],['f'],1)

    def get_maxabs_kernel(proc):
        filename = input_parameters.summed_kernel_directory+get_proc_prefix(proc)+kerneltype+'.bin'
        u = read_bin_file(filename,['f'])[0]  
        u = np.array(u,dtype='f')
        vp=read_bin_file(model_in_dir+get_proc_prefix(proc)+modeltype+'.bin',['f'])[0]
        maxval = max(abs(u))
        return maxval
    results = pprocess.Map(limit=input_parameters.NPROCS,reuse=1)
    calc = results.manage(pprocess.MakeReusable(get_maxabs_kernel))
    for proc in range(input_parameters.NPROCS):
        calc(proc)
    results.finish()
    
    max_kernel_value = 0.0
    for i in range(input_parameters.NPROCS):
        max_kernel_value = max(max_kernel_value,results[i])

    print "Initializing directory for updated model : " + model_out_dir
    os.system('mkdir -p '+model_out_dir)
    
    results = pprocess.Map(limit=input_parameters.NPROCS,reuse=1)
    parfun = results.manage(pprocess.MakeReusable(update_proc))
    for proc in range(input_parameters.NPROCS):        
        parfun(proc,max_kernel_value)
    results.finish()

def write_model_misfit():
    import array
    import input_parameters
    import numpy as np
    from glob import glob
    import os

    it = num_models(copy_flag=False)-1
    exact_model_files = glob('EXACT_MODEL/*vp.bin')
    exact_model_files.sort()
    model_files = glob('models_inversion/m'+str(it)+'/*vp.bin')
    model_files.sort()
    mmisfit = 0.0
    for proc in range(input_parameters.NPROCS):
        u1=read_bin_file(exact_model_files[proc],['f'])[0]
        u2=read_bin_file(model_files[proc],['f'])[0]
        mmisfit += np.sum(np.abs(np.array(u1)-np.array(u2)))        
    if ((it == 0) or (os.path.isfile('model_misfit')==False)):
        modelmisfit = open('model_misfit','w')
        modelmisfit.write("%f 1.0 0\n" % (mmisfit))
        modelmisfit.close()
    else:
        mtemp = open('model_misfit','r')
        mmisfit_initial = np.float(mtemp.readline().split()[0])
        mtemp.close()
        modelmisfit = open('model_misfit','a')
        modelmisfit.write("%f %f %i\n" % (mmisfit,mmisfit/mmisfit_initial,it))
        modelmisfit.close()
