#!/usr/bin/env python
#1. design sky model with one source above and one below horizon
#2. test with OSKAR test example and wsclean
#3. Design one iteration of Moon ULFA
#4. Repeat tests
#5. write program to generate telescope iterations
#6. Combine iterations and test imaging algorithms
#7. Add noise, test calibration 

import sys,os
import cmd
from datetime import datetime, timedelta, date, time
from astropy.io import fits
import string
import numpy as np
import subprocess

#cd /md0/EoR/OSKAR/oskar-data/example-data/OSKAR-2.6-Example-Data

#set the time of the Observation
start_time_string="2018-10-05 03:00:00"
start_time_datetime=datetime.strptime(start_time_string, '%Y-%m-%d %H:%M:%S')

#At this time (11 am AWST) the sun is up and at approx RA DEC: 13:00:00.00, -5:00:00.0
#A source further north with an opposite RA will be below the horizon  try: RA DEC, 02:00:00.00 30:00:00.0
#example sky model format in sky.osm
#make out own sky model test_sky_01.osm
#edit setup.ini (mv original setup.ini to template_setup.ini ))
