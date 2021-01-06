#!/usr/bin/env python
#check pulse pulse period ratios
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from psrqpy import QueryATNF

#download the catalogue
query = QueryATNF()
table = query.table