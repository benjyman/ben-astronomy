import aipy as a
import time
import astropy
import numpy as np
from datetime import datetime
from astropy.time import Time
from scipy import signal
import sys
from SPM import SPM
import struct


spm = SPM(ip='127.0.0.1',firmware='correlate_2020-07-21_1143.fpg')

print("Setting gains prior to quantisation.")
spm.set_gains([16.0])

print("Selecting inputs for correlation.")
spm.set_inputs(0,5)

NFFT = 4096
SITE_LAT = -26.70312
SITE_LON = 116.670575
FREQ_GHz = np.linspace(0.000, 0.250, NFFT);

timestr = time.strftime("%Y%m%d-%H%M%S")
mir_file = str(timestr)+".vis"
# mir_file = "test.vis"
print ("Writing out "+ mir_file)
uv = a.miriad.UV(mir_file, 'new')
uv['history'] = 'test file\n'

uv.add_var('latitud','d')
uv.add_var('npol','i')
uv.add_var('nspect', 'i')
uv.add_var('obsdec', 'd')
uv.add_var('vsource', 'r')
uv.add_var('ischan', 'i')
uv.add_var('operator', 'a')
uv.add_var('nants', 'i')
uv.add_var('baseline', 'r')
uv.add_var('sfreq', 'd')
uv.add_var('inttime', 'r')
uv.add_var('source', 'a')
uv.add_var('epoch', 'r')
uv.add_var('version', 'a')
uv.add_var('ra', 'd')
uv.add_var('restfreq', 'd')
uv.add_var('nschan', 'i')
uv.add_var('sdf', 'd')
uv.add_var('corr', 'r')
uv.add_var('freq', 'd')
uv.add_var('longitu', 'd')
uv.add_var('nchan', 'i')
uv.add_var('tscale', 'r')
uv.add_var('antpos', 'd')
uv.add_var('telescop', 'a')
uv.add_var('pol', 'i')
uv.add_var('coord', 'd')
uv.add_var('veldop', 'r')
uv.add_var('lst', 'd')
uv.add_var('time', 'd')
uv.add_var('dec', 'd')
uv.add_var('obsra', 'd')

uv['latitud'] = SITE_LAT*np.pi/180.0
uv['npol'] = 1
uv['nspect'] = 1
uv['obsdec'] = 0.0 
uv['vsource'] = 0.0
uv['ischan'] = 0
uv['operator'] = 'J'
uv['nants'] = 1
uv['baseline'] = 0.0
uv['sfreq'] = FREQ_GHz[0]
uv['inttime'] = 1.0
uv['source'] = 'SKY'
uv['epoch'] = 2000.0
uv['version'] = 'A'
uv['ra'] = 0.0
uv['restfreq'] = 0.0
uv['nschan'] = NFFT
uv['sdf'] = FREQ_GHz[1]-FREQ_GHz[0]
uv['corr'] = 0.0 
uv['freq'] = 0.0
uv['longitu'] = SITE_LON*np.pi/180.0
uv['nchan'] = NFFT
uv['tscale'] = 0.0
uv['antpos'] = 0.0
uv['telescop'] = 'J'
uv['pol'] = 1
uv['coord'] = 0.0
uv['veldop'] = 0.0
uv['lst'] = 0.0
uv['time'] = 0.0
uv['dec'] = 0.0
uv['obsra'] = 0.0

beam = a.phs.Beam(FREQ_GHz)

ants = []
ants.append(a.phs.Antenna(0,0,0,beam,delay=0,offset=0))
ants.append(a.phs.Antenna(0,3,0,beam,delay=0,offset=0))
aa = a.phs.AntennaArray(ants=ants,location=("-26:42:11.23", "116:40:14.07"))

try:
    while True:
        time.sleep(1.0)
        ut = Time(datetime.utcnow(), scale='utc')
        print (ut)
        data_mask = np.zeros(NFFT)
        aa.set_jultime(ut.jd)
        uv['lst'] = float(aa.sidereal_time())

        signal_cross = spm.get_cross()
        signal_auto_a = spm.get_auto_a()
        signal_auto_b = spm.get_auto_b()

        auto_11 = np.ma.array(signal_auto_a, mask=data_mask, dtype=np.complex64)
        auto_22 = np.ma.array(signal_auto_b, mask=data_mask, dtype=np.complex64)
        cross_12 = np.ma.array(signal_cross, mask=data_mask, dtype=np.complex64)

        uvw_11 = np.array([1,2,3], dtype=np.double)
        preamble = (uvw_11, ut.jd, (0,0)) 
        uv.write(preamble,auto_11)

        uvw_22 = np.array([1,2,3], dtype=np.double)
        preamble = (uvw_22, ut.jd, (1,1)) 
        uv.write(preamble,auto_22)

        uvw_12 = np.array([1,2,3], dtype=np.double)
        preamble = (uvw_12, ut.jd, (0,1)) 
        uv.write(preamble,cross_12)
    
except KeyboardInterrupt:
    print('Exiting')
    del(uv)
    sys.exit

#

