from mwapy import ephem_utils
import mwasources as ms
import sys
import numpy,math,numpy.ma
import os,time,datetime,calendar
import ephem
import logging
from optparse import OptionParser

sec_const = 315964783

# length of individual scan
block_length=112
padding_time = 8
# calibrator scan
cal_block_length=112

# do this source before the target
# might change this to 3C353
pre_calibrator='HerA'
# frequency bands in channels
frequencies=[93,121,145]
project='G0004'
target_name='GP336'
options='--useazel --usegrid='
#creator=os.environ['USER']
creator='jmiller-jones'
starttime=None
command_base='/home/mwa/bin/single_observation.py'

def UTCToGPS(year, month, day, hour, minu, sec) :
    t = calendar.timegm((year, month, day, hour, minu, sec, -1, -1, 0)) 
    return t - sec_const

def gal2eq(l,b):
    """Convert (l,b) -> (ra,dec)
    both in degrees
    """
    gal = ephem.Galactic(numpy.radians(l),numpy.radians(b))
    pos = ephem.Equatorial(gal)
    ra = numpy.degrees(pos.ra)
    dec = numpy.degrees(pos.dec)
    return (ra,dec)

def eq2gal(ra,dec):
    """Convert (ra,dec) -> (l,b)
    both in degrees
    """
    pos = ephem.Equatorial(numpy.radians(ra),numpy.radians(dec))
    gal = ephem.Galactic(pos)
    l = numpy.degrees(gal.lon)
    b = numpy.degrees(gal.lat)
    return (l,b)

def makeSources():
    source_id = 'GP336_{0:d}'
    l_cen=336
    b=0
    sources=[]
    src_count=0
    #100Mhz
    for l in [l_cen + 42*0.5 ,l_cen - 42*0.5]:
        ra,dec = gal2eq(l,b)
        sources.append(ms.Source(name=source_id.format(src_count),ra=ra,dec=dec,freqs=[frequencies[0]],dwelltime_sec=block_length,iscal=False))
        src_count+=1
    #150MHz
    for l in [ l_cen,l_cen+28, l_cen -28]:
        ra,dec = gal2eq(l,b)
        sources.append(ms.Source(name=source_id.format(src_count),ra=ra,dec=dec,freqs=[frequencies[1]],dwelltime_sec=block_length,iscal=False))
        src_count+=1
    #200MHz
    for l in [ l_cen - 21*1.5, l_cen-21*0.5, l_cen+21*0.5, l_cen+21*1.5]:
        ra,dec= gal2eq(l,b)
        sources.append(ms.Source(name=source_id.format(src_count),ra=ra,dec=dec,freqs=[frequencies[2]],dwelltime_sec=block_length,iscal=False))    
        src_count+=1
    return sources

usage="Usage: %prog [options]\n"
usage+="\tSchedules a block of observing\n"
usage+="\tEach block consists of a calibrator scan at each frequency, followed by a compound field\n"
usage+="\tThe pre-calibrator may be omitted\n"
parser = OptionParser(usage='')
parser.add_option('--starttime',dest='starttime',default=None,
                  help='time at start of block (YYYY-MM-DD,hh:mm:ss) [required]')
parser.add_option('--lst',dest='lst',default=False,action="store_true",
                  help='Interpret time as LST? [default=UT]')
parser.add_option('--ut',dest='lst',default=False,action="store_false",
                  help='Interpret time as UT? [default=UT]')
parser.add_option('--source',dest='source',default=None,
                  help='science target [default=%default]')
parser.add_option('--ra',dest='ra',default=None,
                  help='RA of science target [default=%default]')
parser.add_option('--dec',dest='dec',default=None,
                  help='Dec of science target [default=%default]')
parser.add_option('--precal',dest='pre_calibrator',default=pre_calibrator,
                  help='calibrator at start of block [default=%default]')
parser.add_option('--callength',dest='cal_block_length',default=cal_block_length,type='int',
                  help='(s) length of calibrator observation [default=%default]')
parser.add_option('--sourcelength',dest='block_length',default=block_length,type='int',
                  help='(s) length of source observation [default=%default]')
parser.add_option('--frequencies',dest='frequencies',default=frequencies,
                  help='list of center channels [default=%default]')
parser.add_option('--project',dest='project',default=project,
                  help='project ID [default=%default]')
parser.add_option('--targetname',dest='target_name',default=target_name,
                  help='Name for target scan [default=%default]')
parser.add_option('--creator',dest='creator',default=creator,
                  help='creator [default=%default]')
parser.add_option('--options',dest='options',default=options,
                  help='additional options [default=%default]')
parser.add_option('--command',dest='command_base',default=command_base,
                  help='single_observation command [default=%default]')

(options, args) = parser.parse_args()


if options.starttime is None:
    print 'Must supply starttime'
    sys.exit(1)

if isinstance(options.frequencies,str):
    options.frequencies=[int(x) for x in options.frequencies.split(',')]

sources = makeSources()

starttime_datetime0=datetime.datetime.strptime(options.starttime,'%Y-%m-%d,%H:%M:%S')
if options.lst:
    # need to translate the input LST into a UT
    t=ephem_utils.MWATime(datetime=starttime_datetime0)
    LST_target=(t.hour+t.minute/60.0+t.second/3600.0)
    for iter in xrange(2):
        dt=float(t.LST)/15.0-LST_target
        t.datetime-=datetime.timedelta(hours=dt)
    # make sure it is an integer GPStime that is mod 8
    t.gpstime=(int(t.gpstime) / 8)*8
    print "Determined that %s UT = %s LST" % (t.strftime('%H:%M:%S'),t.LST)
    starttime_datetime0=t.datetime
starttime_datetime=starttime_datetime0

time_for_cal=len(options.frequencies)*options.cal_block_length
time_for_source=len(options.frequencies)*options.block_length

commands=[]

totaltime=0
if options.pre_calibrator is not None and len(options.pre_calibrator)>0:
    # do the pre-calibrator
    for frequency in options.frequencies:
        command=options.command_base + ' --starttime=' + starttime_datetime.strftime('%Y-%m-%d,%H:%M:%S') + ' --stoptime=++%ds' % options.cal_block_length
        command+=' --freq=%d,24' % frequency
        command+=' --obsname=%s_%d' % (options.pre_calibrator,frequency)
        command+=' --creator=' + options.creator
        command+=' ' + options.options
        command+=' --project=%s' % options.project
        command+=' --source=%s' % options.pre_calibrator
        command+=' --calibration --calibrator=%s' % options.pre_calibrator
        commands.append(command)
        totaltime+=options.cal_block_length
        totaltime+=padding_time
        starttime_datetime+=datetime.timedelta(seconds=options.cal_block_length)
        starttime_datetime+=datetime.timedelta(seconds=padding_time)

# do the source
for src in sources:
    l,b=eq2gal(src.ra_degs,src.dec_degs)

    command=options.command_base + ' --starttime=' + starttime_datetime.strftime('%Y-%m-%d,%H:%M:%S') + ' --stoptime=++%ds' % options.block_length
    command+=' --freq=%d,24' % src.freqs[0]
    command+=' --obsname=%s_l%3.1f_b0_f%d' % (options.target_name,l,src.freqs[0])
    command+=' --creator=' + options.creator
    command+=' ' + options.options
    command+=' --project=%s' % options.project
    command+=' --ra=%f' % src.ra_degs
    command+=' --dec=%f' % src.dec_degs
    commands.append(command)
    totaltime+=options.block_length
    totaltime+=padding_time
    starttime_datetime+=datetime.timedelta(seconds=options.block_length)
    starttime_datetime+=datetime.timedelta(seconds=padding_time)

print '\n'.join(commands)
#print '\nTotal time for block = %dmin%ds' % (totaltime/60, totaltime%60)
