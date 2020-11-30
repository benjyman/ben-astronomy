#



def write_woden_sims_sbatch_file(nbands,daniel=False,time_string='',pol_list=['X','Y']):
   nbands = int(nbands)
   start_freq_MHz = 50.0
   if (daniel):
       coarse_band_width = '1.28e+6'
       print('coarse_band_width daniel is %s ' % (coarse_band_width))       
   else:
       coarse_band_width = '1.00e+6'
       print('coarse_band_width is %s ' % (coarse_band_width)) 
   #array_layout = "/astro/mwaeor/bmckinley/code/ben-astronomy/AAVS-1/AAVS1_loc_uvgen_match_daniel_255.ant"
   array_layout = "/astro/mwaeor/bmckinley/code/ben-astronomy/EoR/ASSASSIN/eda2_antenna_order_daniel_255.txt"
   year,month,day,hour,min,sec = time_string.split('_')
   time_formatted = '%d-%02d-%02dT%02d:%02d:%02d' % (float(year),float(month),float(day),float(hour),float(min),float(sec))
   type_list = ["gsm","gsm_uniform","EDGES_uniform","unity_uniform","angular"]
   #type_list = ["gsm"]
   #LST_deg = 58.13223745343605    #from template metafits RA=LST
   LST_deg = 60.0
   name_base = "assassin_sims_sbatch_lst_%0.3f" % (LST_deg)  
   sbatch_filename = "%s.sh" % name_base

   with open('%s' % sbatch_filename,'w') as outfile:
      outfile.write("#!/bin/bash --login\n#SBATCH --nodes=1\n#SBATCH --partition=gpuq\n#SBATCH --gres=gpu:1\n")
      outfile.write("#SBATCH --time=00:30:00\n#SBATCH --account=mwaeor\n#SBATCH --nodes=1\n#SBATCH --mem=10gb\n")
      outfile.write("#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n#SBATCH --array=0-%s\n\n" % str(nbands-1))
      #for WODEN
      outfile.write("module use /pawsey/mwa/software/python3/modulefiles\nmodule load erfa/1.7.0\n")
      outfile.write("module load json-c/0.14\nmodule load hdf5/1.10.5\nmodule load cfitsio/3.48\nmodule load cmake/3.15.0\n")
      outfile.write("module load cuda/10.2\nmodule load pal/0.9.8\nmodule load python/3.8.2\nmodule load astropy/4.0.1.post1\n\n")
      #for woden_sourcelists.py
      outfile.write("module load scipy/1.4.1\nmodule load h5py/2.10.0\nmodule load healpy/1.13.0\nmodule load numpy/1.18.2\n\n")
      outfile.write("module load ephem/3.7.7.1\nmodule load matplotlib/3.2.1\n\n")

      outfile.write("source /astro/mwaeor/jline/software/WODEN/build/init_WODEN.sh\n")
      outfile.write("export LD_LIBRARY_PATH=$ERFA_LIB:/pawsey/mwa/software/python3/json-c/0.14-20200419/lib64:$LD_LIBRARY_PATH\n\n")
      if daniel:
         outfile.write("cd /astro/mwaeor/bmckinley/EoR/ASSASSIN/WODEN/noise_coupling\n")
      else:
         outfile.write("cd /astro/mwaeor/bmckinley/EoR/ASSASSIN/WODEN\n")
      
      outfile.write("mkdir -p source_lists\n")
      outfile.write("cd source_lists\n")
      outfile.write("time python /astro/mwaeor/bmckinley/code/ben-astronomy/EoR/ASSASSIN/woden_sourcelists.py --daniel --band=${SLURM_ARRAY_TASK_ID} \n")  
      outfile.write("cd ..\n")
      outfile.write("mkdir -p data\n")
      outfile.write("module swap gcc gcc/5.5.0\n")
      for type in type_list:
         if (type=="gsm" or type=="EDGES_uniform" or type=="unity_uniform"):
            sourcelist_name = "source_lists/woden_map_start_freq_%0.3f_band_${SLURM_ARRAY_TASK_ID}_hpx_%s_sourcelist.txt" % (start_freq_MHz,type)
            output_uvfits_prepend = "data/woden_LST_%0.3f_%s_start_freq_%0.3f" % (LST_deg,type,start_freq_MHz)
            outfile.write("time python /astro/mwaeor/jline/software/WODEN/build/run_woden.py \\\n")
            outfile.write("   --ra0=%0.5f --dec0=-26.70 \\\n" % LST_deg)
            outfile.write("   --num_freq_channels=1 --num_time_steps=1 \\\n")
            outfile.write("   --freq_res=10e+3 --time_res=0.28 \\\n")
            outfile.write("   --lowest_channel_freq=%0.3fe+6 \\\n" % start_freq_MHz)
            outfile.write("   --coarse_band_width=%s \\\n" % coarse_band_width)
            outfile.write("   --cat_filename=%s \\\n" % sourcelist_name)
            #outfile.write("   --metafits_filename=/astro/mwaeor/bmckinley/code/ben-astronomy/EoR/ASSASSIN/WODEN/centre_chan_%03d_metafits_ppds.fits \\\n" % centre_chan)
            outfile.write("   --output_uvfits_prepend=%s \\\n" % output_uvfits_prepend)
            outfile.write("   --sky_crop_components \\\n")
            #outfile.write("   --EDA2_sim \\\n")
            outfile.write("   --primary_beam=EDA2 \\\n")
            outfile.write("   --array_layout=%s \\\n" % array_layout)
            outfile.write("   --band_nums=${SLURM_ARRAY_TASK_ID} \\\n")
            outfile.write("   --date=%s \\\n" % time_formatted)
            outfile.write("   --chunking_size=2500\n\n")
         else:      
            for pol in pol_list:
               if type =="angular":
                  sourcelist_name = "source_lists/woden_map_start_freq_%0.3f_band_${SLURM_ARRAY_TASK_ID}_hpx_gsm_%s_pol_%s_angular_sourcelist.txt" % (start_freq_MHz,time_string,pol)
                  output_uvfits_prepend = "data/woden_LST_%0.3f_gsm_start_freq_%0.3f_pol_%s_angular" % (LST_deg,start_freq_MHz,pol)
               if type =="gsm_uniform":
                  sourcelist_name = "source_lists/woden_map_start_freq_%0.3f_band_${SLURM_ARRAY_TASK_ID}_hpx_gsm_uniform_%s_pol_%s_sourcelist.txt" % (start_freq_MHz,time_string,pol)
                  output_uvfits_prepend = "data/woden_LST_%0.3f_gsm_uniform_start_freq_%0.3f_pol_%s" % (LST_deg,start_freq_MHz,pol)            
               outfile.write("time python /astro/mwaeor/jline/software/WODEN/build/run_woden.py \\\n")
               outfile.write("   --ra0=%0.5f --dec0=-26.70 \\\n" % LST_deg)
               outfile.write("   --num_freq_channels=1 --num_time_steps=1 \\\n")
               outfile.write("   --freq_res=10e+3 --time_res=0.28 \\\n")
               outfile.write("   --lowest_channel_freq=%0.3fe+6 \\\n" % start_freq_MHz)
               outfile.write("   --coarse_band_width=%s \\\n" % coarse_band_width)
               outfile.write("   --cat_filename=%s \\\n" % sourcelist_name)
               outfile.write("   --output_uvfits_prepend=%s \\\n" % output_uvfits_prepend)
               outfile.write("   --sky_crop_components \\\n")
               outfile.write("   --primary_beam=EDA2 \\\n")
               outfile.write("   --array_layout=%s \\\n" % array_layout)
               outfile.write("   --band_nums=${SLURM_ARRAY_TASK_ID} \\\n")
               outfile.write("   --date=%s \\\n" % time_formatted)
               outfile.write("   --chunking_size=2500\n")
                  
   print("wrote %s" % sbatch_filename) 
 
if __name__ == "__main__":
    import argparse
    
    class SmartFormatter(argparse.HelpFormatter):
        def _split_lines(self, text, width):
            if text.startswith('R|'):
                return text[2:].splitlines()
            # this is the RawTextHelpFormatter._split_lines
            return argparse.HelpFormatter._split_lines(self, text, width)
            
            from argparse import RawTextHelpFormatter
            
    parser = argparse.ArgumentParser(description="Run this to make WODEN sky models on \
            supercomputer with array job",formatter_class=SmartFormatter)

    parser.add_argument('--nbands', default='218',
        help='Number of bands as string to go into slurm sbatch file header. Frequencies depend on other arguments such as --daniel e.g. --nbands="218"')
            
    parser.add_argument('--daniel', default=True, action='store_true',
        help='Frequencies correspond to Daniels sims. Band 0 is at 50 MHz and frequency goes up in increments of 1.28 MHz')
            
            
    args = parser.parse_args()
    
    if args.nbands:
       nbands = int(args.nbands)
    #this is for LST 60.0 deg, will do rest of lsts later
    year,month,day,hour,min,sec = 2015,11,29,15,40,29 #LST=60 deg
    time_string = '%d_%02d_%02d_%02d_%02d_%02d' % (year,month,day,hour,min,sec)
    
    write_woden_sims_sbatch_file(nbands=nbands,daniel=args.daniel,time_string=time_string)  
