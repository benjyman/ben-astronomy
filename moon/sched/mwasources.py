#!/usr/bin/python

# A module for known sources
# Randall Wayth. April 2011, July 2012.

# NOTE: dwell times should preferably be multiples of 8 seconds

import math
import ephem

class Source:
    def __init__(self,name,ra,dec,freqs=['95','121','145'],dwelltime_sec=112,el_limit=40.0,known=False,gain=None,iscal=False):
        """Constructor requiring name, ra(degs), dec(degs), PFB coarse chan(s), dwell time (s)"""
        self.ra_degs = ra
        self.dec_degs = dec
        self.name = name
        self.freqs = freqs
        self.dwelltime_sec = dwelltime_sec
        self.el_limit = el_limit
        self.gain = gain
        self.known = known
        self.iscal = iscal

# based on Vincenty formula. See: http://en.wikipedia.org/wiki/Great-circle_distance
    def AngularSep(self ,ra_degs,dec_degs):
        """calc the angular distance between the source and an ra/dec location"""
        d2r = math.pi/180.0
        # compute sines and cosines of start "s" and finish "f" positions.
        slat_s = math.sin(self.dec_degs*d2r)
        slat_f = math.sin(dec_degs*d2r)
        clat_s = math.cos(self.dec_degs*d2r)
        clat_f = math.cos(dec_degs*d2r)
        dlon = self.ra_degs - ra_degs
        cdlon = math.cos(dlon*d2r)
        sdlon = math.sin(dlon*d2r)
        numerator_p1 = clat_f*sdlon
        numerator_p2 = clat_s*slat_f - slat_s*clat_f*cdlon
        numerator = math.sqrt(numerator_p1*numerator_p1 + numerator_p2*numerator_p2)
        denom = slat_s*slat_f + clat_s*clat_f*cdlon
        return math.atan2(numerator,denom)

    def __repr__(self):
        return "< "+self.name+": RA(degs): "+str(self.ra_degs)+". DEC: "+str(self.dec_degs)+". Freqs: "+str(self.freqs)+". El_limit: "+str(self.el_limit)+" degs >"
    

def LoadSources():
  sun = ephem.Sun()
  sun.compute()
  allsrcs = []
  # sun pos when the script is run
  s = Source('Sun',sun.ra*180.0/math.pi,sun.dec*180.0/math.pi,freqs=['62;63;69;70;76;77;84;85;93;94;103;104;113;114;125;126;139;140;153;154;169;170;187;188'],known=True,dwelltime_sec=296,gain=14,el_limit=30)
  allsrcs.append(s)
  s = Source('SCP',0.0,-89.9,known=True,dwelltime_sec=296,gain=0,el_limit=20)
  allsrcs.append(s)
  s = Source('EoR0',0.0,-30.0)
  allsrcs.append(s)
  s = Source('3C33',17.2202442,13.3372789,known=True)
  allsrcs.append(s)
  s = Source('PKS0131-36',23.490583,-36.493250)
  allsrcs.append(s)
  s = Source('ForA',50.673825,-37.208227,known=True)
  allsrcs.append(s)
  s = Source('3C098',59.7268408,10.4341758)
  allsrcs.append(s)
  s = Source('PKS0349-27',57.899000,-27.742972) # Gianni's polarised source
  allsrcs.append(s)
  s = Source('EoR1',60.0,-27.0)
  allsrcs.append(s)
  s = Source('PKS0408-65',62.0849171,-65.7525217,iscal=True,known=True)
  allsrcs.append(s)
  s = Source('PKS0410-75',62.2020517,-75.1220352)
  allsrcs.append(s)
  s = Source('J0437-4715',69.3162,-47.2525)
  allsrcs.append(s)
  s = Source('PKS0442-28',71.157083,-28.165139)
  allsrcs.append(s)
  s = Source('PicA',79.9571708,-45.7788278,iscal=True,known=True)
  allsrcs.append(s)
  s = Source('LMC',80.8939,-69.7561,dwelltime_sec=240)
  allsrcs.append(s)
  s = Source('SMC',13.1583333,-72.800278,dwelltime_sec=240)
  allsrcs.append(s)
  s = Source('Crab',83.633212,22.01446,iscal=True,known=True,el_limit=35)
  allsrcs.append(s)
  s = Source('M42',83.8186621,-5.3896789,known=True)
  allsrcs.append(s)
  s = Source('A3376',90.0,-40.0, freqs=['69','93','117','141','165'])
  allsrcs.append(s)
  s = Source('3C161',96.7921548,-5.8847717)
  allsrcs.append(s)
  s = Source('Bullet',105.0,-56.0,freqs=['69','93','117','141','165'])
  allsrcs.append(s)
  s = Source('PupA',126.0307,-42.9963,known=True)
  allsrcs.append(s)
  s = Source('HydA',139.5236198,-12.0955426,iscal=True,known=True)
  allsrcs.append(s)
  s = Source('TN0924',141.0,-19.0,freqs=['175'])
  allsrcs.append(s)
  s = Source('3C227',146.9381792,7.4223178)
  allsrcs.append(s)
  s = Source('EoR2',170.0,-10.0)
  allsrcs.append(s)
  s = Source('VirA',187.7059304, 12.3911231,iscal=True,known=True)
  allsrcs.append(s)
  s = Source('3C278',193.65083,-12.56306)   # this is close to VirA
  allsrcs.append(s)
  s = Source('3C283',197.9171042,-22.2845667)
  allsrcs.append(s)
  s = Source('CenA',201.3650633,-43.0191125,known=True,iscal=True)
  allsrcs.append(s)
  s = Source('3C317',229.185372,7.021626)
  allsrcs.append(s)
  s = Source('CirX1',230.17021,-57.166694, el_limit=45, dwelltime_sec=56)
  allsrcs.append(s)
  s = Source('PKS1610-60',243.7659992,-60.9071658,known=False)
  allsrcs.append(s)
  s = Source('B1642-03',251.258500, -3.299867,dwelltime_sec=240,known=False)
  allsrcs.append(s)
  s = Source('HerA',252.788,4.99278,iscal=True,known=True)
  allsrcs.append(s)
  s = Source('3C353',260.1173258,-0.9796169,known=True)
  allsrcs.append(s)
  s = Source('SgrA',266.417,-29.0078,gain=5,known=True)
  allsrcs.append(s)
  s = Source('PKS1814-63',274.8958429,-63.7633858)
  allsrcs.append(s)
  s = Source('SS433',287.95685, 4.9830, el_limit=45)
  allsrcs.append(s)
  s = Source('GRS1915',288.79813,10.945778, el_limit=45, dwelltime_sec=56)
  allsrcs.append(s)
  #s = Source('PKS1932-46',293.985649,-46.344601,iscal=True,known=True)
  s = Source('PKS1932-46',293.985649,-46.344601,known=True)
  allsrcs.append(s)
  s = Source('CygA',299.8681525,40.7339156,el_limit=20.0,known=True)
  allsrcs.append(s)
  s = Source('A3667',303.1254,-56.8165,freqs=['69','93','117','141','165'],dwelltime_sec=236)
  allsrcs.append(s)
  s = Source('3C409',303.6149855,23.5813630)
  allsrcs.append(s)
  s = Source('3C433',320.935559,25.069967)
  allsrcs.append(s)
  s = Source('3C444',333.607250,-17.026611,el_limit=30,iscal=True,known=True)
  allsrcs.append(s)
  #s = Source('PKS2153-69',329.2749187,-69.6899125,iscal=True,known=True)
  s = Source('PKS2153-69',329.2749187,-69.6899125,known=True)
  allsrcs.append(s)
  s = Source('PKS2356-61',359.768178,-60.916466,iscal=True,known=True)
  allsrcs.append(s)

  # sources from G0011
  s = Source('J0041-09',10.4574,-9.3493,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J0552-2103',88.2536,-21.0809,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J0817-07',124.3463,-7.5151,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J0909-0939',137.234,-9.6725,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J1202-0650',180.7082,-6.849,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J1215-39',183.89,-38.9929,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J1258-01',194.6899,-1.7697,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J1259-04',194.8641,-4.1915,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J1303-24',195.9278,-24.2806,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J1327-31',201.9823,-31.4983,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J1347-32',206.88,-32.8496,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J1358-47',209.6776,-47.8116,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J1407-5100',211.9007,-51.021,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J1638-64',249.5789,-64.3594,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('PSZ1G018.75',255.5725,-0.9984,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J2012-56',303.0754,-56.8261,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J2034-35',308.6712,-35.8201,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J2152-19',328.0913,-19.6330,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J2154-57',328.5834,-57.8575,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J2201-59',329.5707,-60.4365,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J2210-12',332.5536,-12.1642,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J2246-52',341.5637,-52.7166,freqs=[71,121,169])
  allsrcs.append(s)
  s = Source('J2249-64',342.4462,-64.4186,freqs=[71,121,169])
  allsrcs.append(s)

  # for A00001
  s = Source('SSC',202.0,-37.0,freqs=[78,121,156])
  allsrcs.append(s)

  # for D0000
  allsrcs.append(Source('PKS1718',259.0,-65,freqs=[78,121,156]))

  # for D0007
  allsrcs.append(Source('K2',172.5,1.5,freqs=[121,145],dwelltime_sec=296))
  allsrcs.append(Source('K2_F3',335,-14.0,freqs=[121,145],dwelltime_sec=296))

  # For G0016
  allsrcs.append(Source('WASP-18b',24.3542,-45.678,freqs=[69],dwelltime_sec=240,el_limit=30))
  allsrcs.append(Source('HD14004b',89.9568,-48.2397,freqs=[69],dwelltime_sec=240,el_limit=30))

  
  return allsrcs

def CalSources(allsrcs):
  calsrcs = []
  for s in allsrcs:
    if s.iscal: calsrcs.append(s)
  return calsrcs

def InRA_Range(sources, ra_hrs, ha_range_hrs):
  srcs = []
  # calc ranges of RA that are OK for this source.
  ra_min = ra_hrs - ha_range_hrs
#  if ra_min < 0.0: ra_min += 24.0
  ra_max = ra_hrs + ha_range_hrs
#  if ra_max >= 24.0: ra_max -= 24.0
  for s in sources:
    src_ra_hrs = s.ra_degs*24.0/360.0
    if src_ra_hrs > ra_min and src_ra_hrs < ra_max: srcs.append(s)
  return srcs

def InBeam(sources, beam_ra_degs, beam_dec_degs,beam_radius_degs):
    srcs = []
    for s in sources:
        if s.AngularSep(beam_ra_degs, beam_dec_degs)*180/math.pi < beam_radius_degs:
            srcs.append(s)
    return srcs

def FindCal(sources, ra_degs, za_limit_rad=math.pi/4.0, lat_degs=-26.7):
    """Find the subset of sources within a ZA limit given an RA. Return list of tuples (calsrc, za)"""
    outsrcs = []
    for s in sources:
        ang_sep = s.AngularSep(ra_degs,lat_degs)
        if ang_sep < za_limit_rad:
            outsrcs.append((s,ang_sep))
    # sort the result by zenith angle
    newlist=sorted(outsrcs, key=lambda item: item[1])
    if len(newlist) > 0:
        return newlist
    else:
        return []

if __name__ == "__main__":
  allsrcs = LoadSources()
  calsrcs = CalSources(allsrcs)
  print "There are "+str(len(allsrcs))+" sources."
  print "There are "+str(len(calsrcs))+" cal sources."
  #for s in calsrcs: print s.name
  for lst in range(24):
    cals_this_lst = FindCal(calsrcs,lst*15.0)
    print "Cal sources at LST "+str(lst)
    for s,za in cals_this_lst:
        print "%s ZA: %.1f" % ( s,za*180.0/math.pi)

