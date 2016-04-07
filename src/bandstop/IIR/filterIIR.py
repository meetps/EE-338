import numpy as np
import scipy as sp
import scipy.signal as sg
import matplotlib.pyplot as plt
import pylab as pl

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
filter_number = 82
delta_stop    = 0.15
delta_pass    = 0.15
sampling_freq = 100000
h_transistion = 2
l_transistion = 2
m,q,r         = 7,0,7

stopband_lower_freq  = 4 + 0.7*q + 2*r
stopband_higher_freq = stopband_lower_freq + 10

# Analog Filter Frequencies initialization
omega_p1 = (stopband_lower_freq)*1000.0
omega_p2 = (stopband_higher_freq)*1000.0
omega_s1 = (stopband_lower_freq-l_transistion)*1000.0
omega_s2 = (stopband_higher_freq+h_transistion)*1000.0
analog_freq = np.array([omega_s1,omega_p1,omega_p2,omega_s2],dtype='f')

# Normalized digital frequencies
digital_freq=(analog_freq/sampling_freq)*2*np.pi

# Bilinear Transformation from (-pi,pi) to (-inf,inf)
equiv_analog_freq=np.tan(digital_freq/2)

# Values for frequency transformation
omega_z_s=equiv_analog_freq[0]*equiv_analog_freq[3]
f_B=equiv_analog_freq[3]-equiv_analog_freq[0]

#Frequency Transformation for bandstop
equiv_analog_lowpass_freq=(f_B*equiv_analog_freq)/(omega_z_s-(equiv_analog_freq**2))

# Values for Butterworth low pass filter design
D_1     = (1/(1-delta_stop)**2)-1
D_2     = (1/delta_pass**2)-1
epsilon = np.sqrt(D_1)
mod_equiv_analog_lowpass_freq = abs(equiv_analog_lowpass_freq)

# Order Calculation and getting stringent omega
stringent_omega_s=min(mod_equiv_analog_lowpass_freq[1],mod_equiv_analog_lowpass_freq[2])
N = np.ceil(np.log(np.sqrt(D_2)/np.sqrt(D_1))/np.log(stringent_omega_s/equiv_analog_lowpass_freq[0]))

# Pole calculations
omega_p = equiv_analog_lowpass_freq[0]
omega_c = ((omega_p/(D_1**(1/(2*N))))+(stringent_omega_s/(D_2**(1/(2*N)))))/2
poles   = np.zeros([2*N-1],dtype='complex64')
iterable = ((2*k+1)*np.pi/(2*N) for k in range(2*int(N)))
xp = np.fromiter(iterable,float)
poles = (1.j)*omega_c*np.exp(1.j*xp)

#Lowpass filter transfer function
a=1+0.j
for c in poles:
  if(c.real< 0):
    a=np.poly1d([1,-c],r=0)*a
print "Low pass transfer function denominator\n",a

#Find gain for Butterworth
butter_k=omega_c**N
print "GAIN:",butter_k

'''
For BPF when the transformation is applied to 1/(s-root) we get 
Numerator = s^2 + omega_0^2
Denominator = -c*s^2 + B*s -c*omega_0^2
Using this basic result and finding numerator and denominator to get the bandpass transfer function
'''
analog_numer = butter_k
analog_denom = 1+0.j
for c in poles:
  if(c.real<= 0):
    analog_numer = np.poly1d([1,0,omega_z_s],r=0)*analog_numer 
    analog_denom = np.poly1d([-c,f_B,-c*omega_z_s],r=0)*analog_denom
print "Analog numerator\n",analog_numer
print "\nAnalog denominator\n",analog_denom

# Converting back to digital domain from transfer function
z,p,k=sg.tf2zpk(analog_numer,analog_denom)

plt.figure(5)
plt.grid(True)
plt.scatter(p.real,p.imag,s=50,marker='x')
plt.scatter(z.real,z.imag,s=50,marker='o')
plt.title('Pole Zero plot of Analog Bandstop filter')
plt.ylabel('Imaginary')
plt.xlabel('Real')

'''
For converting the bandpass filter to digital domain using 
s = (1-z^-1)/(1+z^-1)
Numerator = (omega_0^2+1)z^2+2*(omega_0^2-1)z+(omega_0^2+1)
Denominator = (-B-c-c*omega_0^2)z^2+(2*c-2*c*omega_0^2)z+(-c*omega_0^2-c+B)
Using this basic result and finding numerator and denominator to get the bandpass transfer function
'''
digital_numer=butter_k
digital_denom=1+0.j
for c in poles:
  if(c.real<= 0):
    digital_numer=np.poly1d([(omega_z_s+1),2*(omega_z_s-1),(omega_z_s+1)],r=0)*digital_numer 
    digital_denom=np.poly1d([(-f_B-omega_z_s*c-c),(2*c-2*c*omega_z_s),(-c+f_B-c*omega_z_s)],r=0)*digital_denom

print "Analog frequencies :",analog_freq
print "Digital frequencies :",digital_freq
print "Equivalent Digital frequencies : ",equiv_analog_freq
print "Equivalent Analog lpf freq :",equiv_analog_lowpass_freq
print "D1,D2 :",D_1,D_2
print "Poles :",poles
print "Order : ", N
print "stringent_omega :",stringent_omega_s
print "omega_z_s :",omega_z_s
print "B :",f_B

print "Digital Numerator\n",digital_numer 
print "\nDigital Denominator:\n",digital_denom

# Plotting poles of low pass filter
plt.figure(1)
plt.grid(True)
neg_poles=np.zeros([0],dtype='complex64')
for c in poles:
  if(c.real<= 0):
    neg_poles=np.append(neg_poles,c)
plt.scatter(neg_poles.real,neg_poles.imag,s=50,marker='x')
plt.title('Pole Zero plot of Low pass filter')
plt.ylabel('Imaginary')
plt.xlabel('Real')
plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)

nmrz = (digital_numer.c).round(decimals=6)[::-1].real
dmrz = (digital_denom.c).round(decimals=6)[::-1].real

z,p,k = sg.tf2zpk(nmrz,dmrz)

plt.figure(4)
plt.grid(True)
plt.scatter(p.real,p.imag,s=50,marker='x')
plt.scatter(z.real,z.imag,s=50,marker='o')
plt.title('Pole Zero plot of Digital Bandstop filter')
plt.ylabel('Imaginary')
plt.xlabel('Real')

# Printing the latex digital polynomial
print "\nNormalized numerator coefficients array:\n",nmrz/dmrz[0]
printLatexPoly(nmrz/dmrz[0])
print "\nNormalized denominator coefficients array:\n",dmrz/dmrz[0]
printLatexPoly(dmrz/dmrz[0])

#Plot Frequency response
nyq_rate=sampling_freq/2
plt.figure(2)
plt.clf()
plt.grid(True)
w,h= sg.freqz(nmrz,dmrz,worN=512)
plt.plot((w/np.pi)*nyq_rate, np.absolute(h), linewidth=2)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain')
plt.title('Frequency Response')
plt.ylim(-0.05, 1.05)

# Zoomed plot 2
ax1 = plt.axes([0.44, 0.56, .45, .25])
plt.plot((w/np.pi)*nyq_rate, np.absolute(h), linewidth=2)
plt.xlim(16000.0,30000.0)
plt.ylim(-0.01, 0.2)
plt.grid(True)

# Zoomed plot 1
ax1 = plt.axes([0.44, 0.22, .45, .25])
plt.plot((w/np.pi)*nyq_rate, np.absolute(h), linewidth=2)
plt.xlim(28000.0,35000.0)
plt.ylim(0.86, 1.05)
plt.grid(True)

plt.figure(3)
plt.grid(True)
h_Phase = pl.unwrap(np.arctan2(np.imag(h),np.real(h)))
plt.plot(w/max(w),h_Phase)
plt.ylabel('Phase (radians)')
plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
plt.title(r'Phase response')
plt.show()