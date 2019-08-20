#!/usr/bin/env python

import os
import datetime

time_to_split_off_s = 8.

year = '2015'
chan = '145'
obsid = '1117031728'

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


