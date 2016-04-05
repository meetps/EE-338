import numpy as np
import scipy as sp
import pylab as pl
import scipy.signal as sg
import matplotlib.pyplot as plt

filter_number   = 82
delta_stop 		= 0.15
delta_pass    	= 0.15
sampling_freq 	= 100000
m = 7
q = 0
r = 7
passband_lower_freq  = 2 + 0.7*q + 2*r
passband_higher_freq = passband_lower_freq+10

# analog values for bandpass
omega_p1 = (passband_lower_freq)*1000.0
omega_p2 = (passband_higher_freq)*1000.0
omega_s1 = (passband_lower_freq-1) *1000.0
omega_s2 = (passband_higher_freq+1)*1000.0

#analog f_range array
analog_freq = np.array([omega_s1,omega_p1,omega_p2,omega_s2],dtype='f')

#normalized digital frequencies
digital_freq = (analog_freq/sampling_freq)*2*np.pi

#Kaiser window parameters
del_omega1 = digital_freq[3] - digital_freq[2]
del_omega2 = digital_freq[1] - digital_freq[0]
del_omega  = min(abs(del_omega1),abs(del_omega2))

A = -20*np.log10(delta_stop)

#Order
N_float = (A-8)/(2*2.285*del_omega)
N   = np.ceil(N_float)

if(A<21):
    alpha=0
elif(A<=50):
    alpha=0.5842*(A-21)**0.4+0.07886(A-21)
else:
    alpha=0.1102(A-8.7)

omega_c1 = (digital_freq[1]+digital_freq[0])*0.5
omega_c2 = (digital_freq[3]+digital_freq[2])*0.5

#Ideal BPF h[n]
iterable   = ((np.sin(omega_c2*k)-np.sin(omega_c1*k))/(np.pi*k) for k in range(int(-N),int(N+1)))
h_ideal    = np.fromiter(iterable,float)
h_ideal[N] = ((omega_c2-omega_c1)/np.pi)
beta       = alpha/N

#Generating Kaiser window 
h_kaiser = sg.kaiser(2*N+1,beta)
h_org    = h_ideal*h_kaiser

print "Filter Coefficients:\n",h_org

#print "Hideal",h_ideal
#print m,q_m,r_m,passband_lower_freq,passband_higher_freq
#print "Analog frequencies",analog_freq
#print "Digital frequencies",digital_freq
#print "Equivalent Digital frequencies",f_eqv_analog_a
#print "Equivalent Analog lpf freq",f_eqv_analog_lpf_a
#print np.sqrt(D_1),np.sqrt(D_2),B
#print "order=",N 
#print A_k 
#print "Poles",poles_a
#print "Numerator",numer 
#print "Denominator:",denom


#plotting poles

#plt.figure(1)
#plt.grid(True)
#plt.scatter(poles_a.real,poles_a.imag,marker='x')
#plt.legend()

nyquist_rate = sampling_freq/2

#------------------------------------------------
# Plot the FIR filter coefficients.
#------------------------------------------------

plt.figure(1)
plt.plot(h_org, 'bo-', linewidth=2)
plt.title('Filter Coefficients (%d taps)' % (2*N+1))
plt.grid(True)

#Plot Frequency response
plt.figure(2)
plt.clf()
plt.grid(True)
w,h= sg.freqz(h_org)
plt.plot((w/np.pi)*nyquist_rate, np.absolute(h), linewidth=2)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain')
plt.title('Frequency Response')
plt.ylim(-0.05, 1.2)

# Upper inset plot.
ax1 = plt.axes([0.42, 0.6, .45, .25])
plt.plot((w/np.pi)*nyquist_rate, np.absolute(h), linewidth=2)
plt.xlim(8000.0,13000.0)
plt.ylim(0.9, 1.2)
plt.grid(True)

# Lower inset plot
ax2 = plt.axes([0.42, 0.25, .45, .25])
plt.plot((w/np.pi)*nyquist_rate, np.absolute(h), linewidth=2)
plt.xlim(13000.0, 20000.0)
plt.ylim(0.0, 0.11)
plt.grid(True)

plt.figure(3)
plt.grid(True)
h_Phase = pl.unwrap(np.arctan2(np.imag(h),np.real(h)))
plt.plot(w/max(w),h_Phase)
plt.ylabel('Phase (radians)')
plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
plt.title(r'Phase response')

#Stem Diagram
plt.figure(4)
y = pl.linspace(0,61,61)
plt.stem(y,h_org,linefmt='b-', markerfmt='bo', basefmt='r-')
plt.title('Filter Coefficients (%d taps)' % (2*N+1))
plt.grid(True)
plt.show()