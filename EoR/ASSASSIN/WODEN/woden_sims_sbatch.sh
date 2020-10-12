#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=2,3

module swap gcc gcc/5.5.0
module use /pawsey/mwa/software/python3/modulefiles
module load erfa/1.7.0
module load json-c/0.14
module load hdf5/1.10.5
module load cfitsio/3.48
module load cmake/3.15.0
module load cuda/10.2
module load pal/0.9.8

module load python/3.8.2
module load astropy/4.0.1.post1
#export PYTHONPATH=$PYTHONPATH:/astro/mwaeor/jline/software/local_python
source /astro/mwaeor/jline/software/WODEN_EDA2/build/init_WODEN.sh
export LD_LIBRARY_PATH=$ERFA_LIB:/pawsey/mwa/software/python3/json-c/0.14-20200419/lib64:$LD_LIBRARY_PATH

latitude = -26.70331940 (#this is actually MWA)

#example sbatch script lives here: /astro/mwaeor/jline/test_WODEN/EDA2_sims/sbatch_run_woden.sh
#cd /astro/mwaeor/jline/test_WODEN/EDA2_sims
#time python /astro/mwaeor/jline/software/WODEN_EDA2/build/run_woden.py \
#    --ra0=60.0 --dec0=-30.0 \
#    --num_freq_channels=16 --num_time_steps=14 \
#    --freq_res=80e+3 --time_res=8.0 \
#    --cat_filename=/astro/mwaeor/jline/test_WODEN/pumav3_foregrounds/woden-srclist_pumav3.txt \
#    --metafits_filename=/astro/mwaeor/jline/test_WODEN/pumav3_foregrounds/1132846440_metafits_ppds.fits \
#    --output_uvfits_prepend=/astro/mwaeor/jline/test_WODEN/EDA2_sims/data/EDA2_EoR1_pumav3 \
#    --sky_crop_components \
#    --EDA2_sim \
#    --array_layout=/astro/mwaeor/jline/test_WODEN/EDA2_sims/AAVS1_loc_uvgen_255.ant \
#    --band_nums=$SLURM_ARRAY_TASK_ID \
#    --chunking_size=5000 \
#    # --no_tidy

#python /astro/mwaeor/jline/software/WODEN_EDA2/build/run_woden.py \
#    --ra0=metafits_LST --dec0=latitude \
#    --num_freq_channels=1 --num_time_steps=1 \
#    --freq_res=80e+3 --time_res=0.28 \
#    --cat_filename=some_catalogue_frequency_scaled_$SLURM_ARRAY_TASK_ID  \
#    --metafits_filename=/astro/mwaeor/jline/test_WODEN/haslam_map/1132846440_metafits_ppds.fits \
#    --output_uvfits_prepend=your_directory_and_name_here \
#    --sky_crop_components \
#    --EDA2_sim \
#    --telescope_name=EDA2 \
#    --array_layout=AAVS1_loc_uvgen.ant \
#    --band_nums=$SLURM_ARRAY_TASK_ID \
#    --chunking_size=2500
    
    
cd /astro/mwaeor/jline/test_WODEN/EDA2_sims
time python /astro/mwaeor/jline/software/WODEN_EDA2/build/run_woden.py \
    --ra0=58.13223745343605 --dec0=-26.70 \
    --num_freq_channels=1 --num_time_steps=1 \
    --freq_res=10e+3 --time_res=0.28 \
    --cat_filename=/astro/mwaeor/bmckinley/EoR/ASSASSIN/WODEN/hpx_fits/woden_$SLURM_ARRAY_TASK_ID.txt \
    --metafits_filename=/astro/mwaeor/bmckinley/EoR/ASSASSIN/WODEN/1132846440_metafits_ppds.fits \
    --output_uvfits_prepend=/astro/mwaeor/bmckinley/EoR/ASSASSIN/WODEN/uniform_sky_test \
    --sky_crop_components \
    --EDA2_sim \
    --array_layout=/astro/mwaeor/jline/test_WODEN/EDA2_sims/AAVS1_loc_uvgen_255.ant \
    --band_nums=$SLURM_ARRAY_TASK_ID \
    --chunking_size=5000 \
    # --no_tidy