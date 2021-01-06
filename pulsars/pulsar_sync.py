#!/usr/bin/env python3
#check pulse pulse period ratios
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from psrqpy import QueryATNF

#download the catalogue
#query = QueryATNF()
#query.save('atnfquery.pkl')

oldquery = QueryATNF(loadquery='atnfquery.pkl')

numstring = 'Version {} of the ATNF catalogue contains {} pulsars'
print(numstring.format(oldquery.get_version, oldquery.num_pulsars))



#pandas
df = oldquery.pandas

num_pulsars_used = 20
pulsars_used = df.iloc[0:num_pulsars_used,:]

print(pulsars_used)

print(df.columns)


#calculate the ratio of the pulse periods and the distance separating each pulsar pair
for pulsar1_index in range(0,num_pulsars_used):
   pulsar1_name = df.loc[pulsar1_index,'PSRJ']
   #pulsar1_period = 
   for pulsar2_index in range(0,num_pulsars_used):
      #only uniques baselines and dont include autos
      if (pulsar2_index > pulsar1_index):
         pulsar2_name = df.loc[pulsar2_index,'PSRJ']
         print("pulsar1 name %s, pulsar2 name %s" % (pulsar1_name,pulsar2_name))   