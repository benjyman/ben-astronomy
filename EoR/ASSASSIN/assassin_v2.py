#!/usr/bin/env python3
#ASSASSIN: (All Sky SignAl Short Spacing INterferometer)
#v2 uses complex beam patterns to simulate model vis in ms format

import matplotlib
matplotlib.use('Agg')
import os,sys
import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.time import TimeDelta
import astropy.units as u
#from PIL import Image
#from PIL import ImageFont
#from PIL import ImageDraw 
import healpy as hp
from pygdsm import GSMObserver
from pygdsm import GlobalSkyModel
from pygdsm import GlobalSkyModel2016
from datetime import datetime, date
import time
