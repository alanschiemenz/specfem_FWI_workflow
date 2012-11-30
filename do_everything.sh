echo `date` 

# Assumes MESH has been created and placed in folder ./meshes/
# Assumes that all specfem3D binaries have been placed in folder ./binary_files/

# Create source and receiver files
cd input_files
python write_CMTSOLUTION_files.py # modify this for different sources
python write_STATIONS.py          # modify this for difference receivers
cd ../

# Create true and initial models
# For different models, modify the file create_xyz_<starting/true>.py
cd models
python driver_models.py
cd ../

# Generate "data" from true model
cd generate_data
python driver_generate_data.py
cd ../

# Performs acoustic FWI with initial model from above
# L2 misfit is applied
# steepest descent with constant step length is taken
# Inverts only for Vp (vs=0, rho is assumed known)
cd inversion
python driver_inversion.py
cd ../

echo `date`
