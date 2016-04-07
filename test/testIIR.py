import numpy as np
import scipy as sp
import pylab as pl
import scipy.signal as sg
import matplotlib.pyplot as plt

# Passband Test
filter_number = 82
delta_stop 	  = 0.15
delta_pass 	  = 0.15
sampling_freq 	  = 100000
h_transistion = 2
l_transistion = 2
m,q,r         = 7,0,7

passband_lower_freq  = 4 + 0.7*q + 2*r
passband_higher_freq = passband_lower_freq + 10

# Analog Filter Frequencies initialization
omega_p1 = (passband_lower_freq)*1000.0
omega_p2 = (passband_higher_freq)*1000.0
omega_s1 = (passband_lower_freq-l_transistion)*1000.0
omega_s2 = (passband_higher_freq+h_transistion)*1000.0

analog_freq = np.array([omega_s1,omega_p1,omega_p2,omega_s2],dtype='f')
digital_freq = (analog_freq/sampling_freq)*2*np.pi

tol = -20*np.log(0.15)
wc1 = (digital_freq[1] + digital_freq[0])/2
wc2 = (digital_freq[2] + digital_freq[3])/2
b,a = sg.iirfilter(4, wc1, wc2,tol,tol,btype='bandpass',analog=True,ftype='cheby1')

print b
print a