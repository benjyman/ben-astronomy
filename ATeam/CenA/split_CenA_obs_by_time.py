#!/usr/bin/env python

import os
import datetime

time_to_split_off_s = 16.

year_list = ['2015','2018']
chan_list = ['145','169']

#year = '2015'
#chan = '145'
#obsid_list = ['1117031728','1121334536','1121420968','1121507392','1121593824','1121680256']

#year = '2015'
#chan = '169'
#obsid_list = ['1117031848','1121248288','1121334712','1121421144','1121507576','1121594000']

#year = '2018'
#chan = '145'
#obsid_list = ['1199663088','1200604688','1200691136','1200777584','1200864032','1200950480']

#year = '2018'
#chan = '169'
#obsid_list = ['1200431968','1200518416','1200604864','1200691312','1200777760','1200864208']

for year in year_list:
   for chan in chan_list:
      if year=='2015':
         if chan=='145':
            obsid_list = ['1117031728','1121334536','1121420968','1121507392','1121593824','1121680256']
         elif chan=='169':
            obsid_list = ['1117031848','1121248288','1121334712','1121421144','1121507576','1121594000']
      elif year=='2018':
         if chan=='145':
            obsid_list = ['1199663088','1200604688','1200691136','1200777584','1200864032','1200950480']
         elif chan=='169':
            obsid_list = ['1200431968','1200518416','1200604864','1200691312','1200777760','1200864208']
      
      for obsid in obsid_list:
         input_ms_name = "/md0/ATeam/CenA/%s/%s/%s/%s.ms" % (year,chan,obsid,obsid)
         output_ms_name = "/md0/ATeam/CenA/%s/%s/%s/%s_split.ms" % (year,chan,obsid,obsid)
         
         
         listfile_name = "/md0/ATeam/CenA/%s_info.txt" % obsid
         casa_script_filename = "casa_cmd.sh"
         casa_cmd_string = "listobs(vis='%s',listfile='%s')" % (input_ms_name,listfile_name)
         with open(casa_script_filename,'w') as f:
            f.write(casa_cmd_string)
         
         cmd = "rm -rf %s" % listfile_name
         print cmd
         os.system(cmd)
         
         cmd = "casa -nologger -nogui -c %s" % casa_script_filename
         print cmd
         os.system(cmd)
         
         with open(listfile_name,'r') as f:
            lines = f.readlines()
         
         for line in lines:
            if ("Observed from" in line):
               start_time = line.split('from')[1].split()[0].split('.')[0]
         
         #30-May-2015/14:35:13.0
         start_time_datetime = datetime.datetime.strptime(start_time, "%d-%b-%Y/%H:%M:%S") 
         time_to_split_off_s_timedelta = datetime.timedelta(seconds=time_to_split_off_s)
         end_time_datetime = start_time_datetime + time_to_split_off_s_timedelta
         
         #'YYYY/MM/DD/hh:mm:ss~YYYY/MM/DD/hh:mm:ss'
         start_time_string_casa = start_time_datetime.strftime("%Y/%m/%d/%H:%M:%S")
         end_time_string_casa = end_time_datetime.strftime("%Y/%m/%d/%H:%M:%S")
         
         timerange_string_casa = "%s~%s" % (start_time_string_casa,end_time_string_casa)
         
         print timerange_string_casa
         
         casa_script_filename = "casa_cmd_split.sh"
         casa_cmd_string = "split(vis='%s',outputvis='%s',timerange='%s')" % (input_ms_name,output_ms_name,timerange_string_casa)
         with open(casa_script_filename,'w') as f:
            f.write(casa_cmd_string)
         
         cmd = "rm -rf %s" % output_ms_name
         print cmd
         os.system(cmd)
         
         cmd = "casa -nologger -nogui -c %s" % casa_script_filename
         print cmd
         os.system(cmd)
      