#script to get obsids for on and off moon pointings and also get the RA DEC of the moon
from datetime import datetime, timedelta
import sys, getopt, string, os
from ephem import *
from time import *
import cmd
from astropy.time import Time
import calendar

sidereal_day_sec = 86164

def obsid_list_generator(options):
   
   if (options.start_date!=None):
      start_date_string=options.start_date
   else:
      exit("Must set a start date. Exiting...")
   start=datetime.strptime(start_date_string, '%Y/%m/%d').timetuple()
   start_date_underscore=start_date_string.replace('/','_')
   
   duration_string=options.duration
   duration=float(duration_string)
   duration_seconds=duration*24*60*60

   project=options.project
   
   if project=="G0017":
      if (options.on_moon =="on"):
         on_moon = options.on_moon
      elif (options.on_moon =="off"):
         on_moon = options.on_moon
         if (options.on_moon_date!=None):
            on_moon_date=options.on_moon_date
            on_moon_date_underscore=on_moon_date.replace('/','_')
         else:
            exit("on_moon_date must be set for off_moon obs for project G0017, exiting...")
      else:
         exit("on_moon not set for project G0017, set --on_moon='on' or --on_moon='off' , exiting...")
   
   #convert to GPS time for MWA database queries
   unix_start_time=calendar.timegm(start)
   gps_start_time= int(unix_start_time - 315964782)
   
   #work out end time in gps time
   gps_end_time=int(gps_start_time + duration_seconds)

   obsid_list_all_text_filename="obsid_list_all_%s_%s.txt" % (start_date_underscore,duration_string)

   #look up all obsids for 
   cmd="obsinfo.py --from=%s --to=%s --projectid=%s > %s " % (gps_start_time,gps_end_time,project,obsid_list_all_text_filename)
   #print cmd
   os.system(cmd)
   
   #clear directory of txt files first
   cmd = "rm -rf %s_%s_*.txt" % (start_date_underscore,duration_string)
   #print cmd
   os.system(cmd)
   
   chan_list=[]
   
   #Go through the file line by line and add the obsids to the appropriate lists by channel.
   for line in open(obsid_list_all_text_filename).readlines():
      #print line
      obsid=line.split()[0][0:-1]
      obs_name=line.split()[1].split('(')[0]
      obs_name_prefix=obs_name.split('_')[0]
      chan=obs_name.split('_')[1]
      if (project=="G0017" and obs_name_prefix=="EoRMoon"):
         if chan in chan_list:
            pass
         else:
            chan_list.append(chan)
         if (on_moon == "on"):
            obsid_output_filename="%s_%s_on_moon_obsids_%s.txt" % (start_date_underscore,duration_string,chan)
            #print "appending to file %s " % (obsid_output_filename)
         elif (on_moon=="off"):
            obsid_output_filename="%s_%s_off_moon_obsids_%s.txt" % (start_date_underscore,duration_string,chan)
            #print "appending to file %s " % (obsid_output_filename)
         fileprint = open(obsid_output_filename, 'a')
         print >> fileprint, '%s' %  (obsid)                               
    
   
   #If it is off moon, now go through the two files, match up the on-moon and off-moon obsids and add in the rad dec of the moon
   if (project=="G0017" and on_moon=="off"):
      print chan_list
      for chan in chan_list:
         track_off_moon_filename="track_off_moon_%s_%s_%s.txt" % (on_moon_date_underscore,start_date_underscore,chan)
         #clear previous
         cmd="rm -rf %s" % track_off_moon_filename
         os.system(cmd)
         print track_off_moon_filename
         off_moon_obsid_filename="%s_%s_off_moon_obsids_%s.txt" % (start_date_underscore,duration_string,chan)
         for off_moon_obsid in open(off_moon_obsid_filename).readlines():
            #find the matching obsid in the on-moon file, i.e. where the LST is the same
            on_moon_obsid_filename="%s_%s_on_moon_obsids_%s.txt" % (on_moon_date_underscore,duration_string,chan)
            continue_search=True
            for on_moon_obsid in open(on_moon_obsid_filename).readlines():
               if continue_search:
                  off_moon_obsid=off_moon_obsid.strip()
                  on_moon_obsid=on_moon_obsid.strip()
                  if (len(off_moon_obsid)==10 and len(on_moon_obsid)==10):
                     #one sidereal day is 86164 s
                     LST_difference=float(off_moon_obsid)-float(on_moon_obsid)
                     LST_remainder = LST_difference % sidereal_day_sec
                     if (abs(LST_remainder/sidereal_day_sec) < 0.0001 or abs(1-LST_remainder/sidereal_day_sec) < 0.0001):
                        #print LST_remainder/sidereal_day_sec
                        print "LSTs match for obsids on_moon %s off_moon %s" % (on_moon_obsid,off_moon_obsid)
                        #Find the Moon position  
                        date = Time(int(on_moon_obsid), format='gps')
                        print date
                        new_date = Time(date, format='iso')
                        print new_date
                        #get date in right format for print_src.py eg --date='2015/3/2 12:01:01'
                        #new_date=string.replace(string.replace(date,'-','/'),'T',' ')
                        #print new_date
                        #find position of moon
                        src_file_name='src_file_moon_%s.txt' % on_moon_obsid 
                        cmd="print_src.py --date='"+str(new_date)+"' > "+ src_file_name
                        print cmd
                        os.system(cmd)
                        #read src_file.txt to get Moon ra and dec 
                        with open(src_file_name, "r") as infile:
                           lines=infile.readlines()
                           moon_ra=lines[6].split()[3].split(',')[0]
                           moon_dec=lines[6].split()[6].split(',')[0]
                           #print moon_ra
                           #print moon_dec
                        fileprint = open(track_off_moon_filename, 'a')
                        print >> fileprint, '%s %s %s %s ' %  (on_moon_obsid, moon_ra, moon_dec, off_moon_obsid)
                        continue_search=False
                        continue
                     else:
                        #print "LSTs do not match for obsids on_moon %s off_moon %s" % (on_moon_obsid,off_moon_obsid)
                        pass
                  else:
                     continue
               else:
                  continue
            
            
import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: obsid_list_generator.py [options]'

parser = OptionParser(usage=usage)

#parser.add_option('--track_moon',action='store_true',dest='track_moon',default=False,help='Track the Moon by shifting centre of each image to Moon position on sky [default=%default]')
parser.add_option('--start_date',type='string', dest='start_date',default=None,help='Date to start searching for obs YYYY/MM/DD e.g. --start_date="2017/10/05" [default=%default]')
parser.add_option('--duration',type='string', dest='duration',default="1",help='Duration to search for in days e.g. --duration="10" [default=%default]')
parser.add_option('--project',type='string', dest='project',default="G0017",help='Project ID of obs e.g. --project="G0017" [default=%default]')
parser.add_option('--on_moon',type='string', dest='on_moon',default=None,help='Set to on or off if project G0017 e.g. --on_moon=on [default=%default]')
parser.add_option('--on_moon_date',type='string', dest='on_moon_date',default=None,help='If off moon, set this to indicate the corresponding on-moon date. e.g. --on_moon_date="2017/10/05" [default=%default]')



(options, args) = parser.parse_args()



obsid_list_generator(options)