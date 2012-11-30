#!/bin/bash
#
# script runs decomposition,database generation and solver
# using this example setup
#
# prior to running this script, you must create the mesh files
# in directory MESH/
#

#./go_solver $1 $2
#$1 = simulation name (synthetics or adjoint)
#$2 = number of processors

###################################################

# number of processes
NPROC=$2

##################################################

echo "running example: `date`"
currentdir=`pwd`

THISHOST=$(hostname) 


simname=$1
echo "Making directory $simname"
mkdir $simname
mkdir -p $simname/in_data_files
cp in_data_files/Par_file in_data_files/CMTSOLUTION $simname/in_data_files/

echo
echo

# runs simulation
echo
echo "  running solver..."
echo
cd bin/
if [ "$THISHOST" == "tethys" ]
then
    mpirun.openmpi -np $NPROC -nolocal -hostfile ~/TETHYS.machines.32 ./xspecfem3D
elif [ "$THISHOST" == "segfault"  -o "$THISHOST" == "coredump" -o "$THISHOST" == "walchen" ]
then
    mpirun.openmpi -np $NPROC  ./xspecfem3D
else
    mpiexec -n $NPROC ./xspecfem3D
    #mpiexec.hydra  -bootstrap ll -n $NPROC -hostfile hostfile$3 ./xspecfem3D
fi



cd ../


if [ "$1" == "synthetics" ]
 then 
    echo "moving results to directory: $simname"
    mkdir -p $simname/in_out_files/OUTPUT_FILES
    mv in_out_files/OUTPUT_FILES/* $simname/in_out_files/OUTPUT_FILES/
fi

#cp xconvolve_source_timefunction convolve_source_timefunction.csh convolve_script $simname/
#halfdur=0.142
#cd $simname
#./convolve_script #$halfdur
#cd ../

echo

echo
echo "done"
echo `date`

