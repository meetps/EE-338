import numpy as np
import scipy as sp
import pylab as pl
import scipy.signal as sg
import matplotlib.pyplot as plt

def printLatexPoly(coeffs):
    '''
    Utility function to print a numpy polynomial into Latex
    taking the coeffs of the polynomial (highest power first)
    as input
    '''    
    latex_poly = ''
    for i in range(len(coeffs)):
        if i == 0:
            latex_poly += str( "%.4f" %coeffs[i]) + ' '
        else: 
            if coeffs[i] > 0.0:
                latex_poly += ' + '
        latex_poly += ' ' + str( "%.4f" %coeffs[i]) + ' Z^{' + str(-i) + '} ' 
    print latex_poly

# Filter Specification Declaration
filter_number   = 82
delta_stop 		= 0.15
delta_pass    	= 0.15
sampling_freq 	= 100000
m,q,r 			= 7,0,7
h_transistion 	= 2
l_transistion 	= 2

passband_lower_freq  = 4 + 0.7*q + 2*r
passband_higher_freq = passband_lower_freq+10

# Analog Filter Frequencies initialization
omega_p1 = (passband_lower_freq)*1000.0
omega_p2 = (passband_higher_freq)*1000.0
omega_s1 = (passband_lower_freq - l_transistion) *1000.0
omega_s2 = (passband_higher_freq + h_transistion)*1000.0
analog_freq = np.array([omega_s1,omega_p1,omega_p2,omega_s2],dtype='f')

# Normalized digital frequencies
digital_freq = (analog_freq/sampling_freq)*2*np.pi

# Kaiser window parameters
del_omega1 = digital_freq[3] - digital_freq[2]
del_omega2 = digital_freq[1] - digital_freq[0]
del_omega  = min(abs(del_omega1),abs(del_omega2))
A 		   = -20*np.log10(delta_stop)

if(A<21):
    alpha=0
elif(A<=50):
    alpha=0.5842*(A-21)**0.4+0.07886(A-21)
else:
    alpha=0.1102(A-8.7)

'''
Order Calculation for the FIR filter
Since a filter designed for the critcal value of 
order N does nto meet the tolerance specifications,
we need to increase the order by an empirical factor 
order_offset which can be concluded by seeing the frequency
response plots and ensuring they meet the tolerance specs. 
'''
order_offset = 5
N_critical = np.ceil((A-8)/(2*2.285*del_omega))
N = N_critical + order_offset

# Cutoff frequency calculation ideal impulse response
omega_c1 = (digital_freq[1]+digital_freq[0])*0.5
omega_c2 = (digital_freq[3]+digital_freq[2])*0.5

# Obtain the ideal bandpass impulse response
iterable   = ((np.sin(omega_c2*k)-np.sin(omega_c1*k))/(np.pi*k) for k in range(int(-N),int(N+1)))
h_ideal    = np.fromiter(iterable,float)
h_ideal[N] = ((omega_c2-omega_c1)/np.pi)
beta       = alpha/N

# Generate Kaiser window 
h_kaiser = sg.kaiser(2*N+1,beta)
h_org    = h_ideal*h_kaiser

print "FIR Filter Coefficients:\n",h_org
print "\nH(Z) = :\n"
printLatexPoly(h_org)
print "\nHideal",h_ideal

# Plot the FIR filter coefficients.
nyquist_rate = sampling_freq/2
plt.figure(1)
plt.plot(h_org, 'bo-', linewidth=2)
plt.title('Filter Coefficients (%d taps)' % (2*N+1))
plt.grid(True)

# Plot Frequency response
plt.figure(2)
plt.clf()
plt.grid(True)
w,h= sg.freqz(h_org)
plt.plot((w/np.pi)*nyquist_rate, np.absolute(h), linewidth=2)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain')
plt.title('Frequency Response')
plt.ylim(-0.05, 1.2)
plt.xlim(0, 100000)

# Zoomed plot 1
ax1 = plt.axes([0.42, 0.6, .45, .25])
plt.plot((w/np.pi)*nyquist_rate, np.absolute(h), linewidth=2)
plt.xlim(16000.0,30000.0)
plt.ylim(0.9, 1.15)
plt.grid(True)

# Zoomed plot 2
ax2 = plt.axes([0.42, 0.25, .45, .25])
plt.plot((w/np.pi)*nyquist_rate, np.absolute(h), linewidth=2)
plt.xlim(30000.0, 34000.0)
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
y = pl.linspace(0,h_org.shape[0],h_org.shape[0])
plt.stem(y,h_org,linefmt='b-', markerfmt='bo', basefmt='r-')
plt.title('Filter Coefficients (%d taps)' % (2*N+1))
plt.grid(True)
plt.show()