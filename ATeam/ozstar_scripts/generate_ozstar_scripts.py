#!/usr/bin/env python
#generate scripts to lauch on ozstar for data processing

def initiate_script(filename):
   header_text = ("
                 #!/bin/bash -l \n
                 #SBATCH --nodes=1 \n
                 #SBATCH --cpus-per-task=8 \n
                 #SBATCH --mem=100000 \n
                 #SBATCH --account=oz048 \n
                 #SBATCH --time=12:00:00 \n
                 #SBATCH --gres=gpu:1 \n
                 #SBATCH --partition=skylake-gpu \n
                 \n
                 module load gcc/6.4.0 openmpi/3.0.0 \n
                 module load fftw/3.3.7 \n
                 module load gsl/2.4 \n
                 module load cfitsio/3.420 \n
                 module load boost/1.67.0-python-2.7.14 \n
                 module load hdf5/1.10.1 \n
                 module load openblas/0.2.20 \n
                 module load cuda/9.0.176 \n
                 \n
                 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fred/oz048/MWA/CODE/lib/
                 ")
                 
   print(header_text)
                
   #with open(filename,'w') as f:
   

def generate_download(obsid_list,dest_dir,timeres=4,freqres=40,ms=True):
   print('generating download script') 
   
def generate_model_cal():
   print('generating model calibration script')

def generate_wsclean_image():
   print('generating wsclean script')
   
def generate_selfcal():
   print('generating selfcal script')
   
def generate_final_sbatch_script():
   print('generating final sbatch script')
   

initiate_script('test.sh')