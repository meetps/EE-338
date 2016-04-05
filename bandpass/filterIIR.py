import numpy as np
import scipy as sp
import pylab as pl
import scipy.signal as sg
import matplotlib.pyplot as plt
import matplotlib.patches as pat

filter_number = 82
delta_stop = 0.15
delta_pass = 0.15
f_sample = 100000
h_transistion = 2
l_transistion = 2

m = 7
q = 0
r = 7
passband_lower_freq  = 4 + 0.7*q + 2*r
passband_higher_freq = passband_lower_freq + 10

# analog values for bandpass
omega_s1=(passband_lower_freq-l_transistion)*1000.0
omega_p1=(passband_lower_freq)*1000.0
omega_p2=(passband_higher_freq)*1000.0
omega_s2=(passband_higher_freq+h_transistion)*1000.0

# analog f_range array
analog_freq=np.array([omega_s1,omega_p1,omega_p2,omega_s2],dtype='f')

# normalized digital frequencies
digital_freq=(analog_freq/f_sample)*2*np.pi

# bilinear transformation
equiv_analog_freq=np.tan(digital_freq/2)

# intermediate values needed for frequency transformation
omega_z_s=equiv_analog_freq[1] * equiv_analog_freq[2]
f_B=equiv_analog_freq[2] - equiv_analog_freq[1]

#Frequency Transformation for bandpass
equiv_analog_lowpass_freq=((equiv_analog_freq**2)-omega_z_s)/(f_B*equiv_analog_freq)

#Values for Chebyschev lpf design
D_1=(1/(1-delta_stop)**2)-1
D_2=(1/delta_pass**2)-1
epsilon=np.sqrt(D_1)


abs_equiv_analog_lowpass_freq=abs(equiv_analog_lowpass_freq)
stringent_omega_s=min(abs_equiv_analog_lowpass_freq[0],abs_equiv_analog_lowpass_freq[3])
N_float = np.arccosh(np.sqrt(D_2)/np.sqrt(D_1)) / np.arccosh(stringent_omega_s/equiv_analog_lowpass_freq[2])
N=np.ceil(N_float)

#print "stringent",stringent_omega_s
#print "omega_z_s",omega_z_s
#print "B",f_B

omega_p  = equiv_analog_lowpass_freq[2]
poles  = np.zeros([2*N-1],dtype='complex64')
iterable = ((2*k+1)*np.pi/(2*N) for k in range(2*int(N)))
A_k = np.fromiter(iterable,float)
B   = np.arcsinh(1/epsilon)/N
poles_real  = omega_p*np.sin(A_k)*np.sinh(B)
poles_imag  = omega_p*np.cos(A_k)*np.cosh(B)
poles       = poles_real + poles_imag*(1.j)

#Lowpass filter transfer function
#Find gain for Chebyshev
a = 1+0.j
chebyshev_k=1

for c in poles:
  if(c.real< 0):
    a=np.poly1d([1,-c],r=0)*a

print "Low pass tf\n",a

for k in range(int(N)):
		chebyshev_k=chebyshev_k*poles[k]

print "GAIN:",chebyshev_k

#Hanalog(s) s_l <--F(s)
#For BPF when the transformation is applied to 1/(s-c)
#we get Bs/(s^2+omega_0^2-B*c*s).Using this basic result and
#finding numerator and denominator separately
#For bandstop numerator and denominator will be interchanged
numer = chebyshev_k.real
denom = 1+0.j
for c in poles:
	if(c.real<= 0):
		numer=np.poly1d([f_B,0],r=0)*numer 
		denom=np.poly1d([1,-f_B*c,omega_z_s],r=0)*denom
#print "analog numerator\n",numer
#print "\nanalog denominator\n",denom
z,p,k=sg.tf2zpk(numer,denom)

plt.figure(5)
plt.grid(True)
plt.scatter(p.real,p.imag,s=50,c='b',marker='x')
plt.scatter(z.real,z.imag,s=50,c='b',marker='o')
plt.title('Pole Zero plot of Analog Bandpass filter')
plt.ylabel('Imaginary')
plt.xlabel('Real')

#Tedious still basic version
#Need to try convolution
#Numerator= B(z^2-1)
#Denominator= (omega_0^2-B+1)z^2+(2*omega_0^2-2)z+(omega_0^2+B+1)

numer=chebyshev_k.real
denom=1+0.j
for c in poles:
  if(c.real<= 0):
    numer=np.poly1d([f_B,0,-f_B],r=0)*numer 
    denom=np.poly1d([(omega_z_s-f_B*c+1),((2*omega_z_s)-2),(omega_z_s+f_B*c+1)],r=0)*denom

z,p,k=sg.tf2zpk(numer,denom)

#print "Check",z,p,k

plt.figure(4)
plt.grid(True)
plt.scatter(p.real,p.imag,s=50,marker='x')
plt.scatter(z.real,z.imag,s=50,marker='o')
plt.title('Pole Zero plot of Digital Bandpass filter')
plt.ylabel('Imaginary')
plt.xlabel('Real')

#print m,q_m,r_m,passband_lower_freq,passband_higher_freq
#print "Analog frequencies",analog_freq
#print "Digital frequencies",digital_freq
#print "Equivalent Digital frequencies",equiv_analog_freq
#print "Equivalent Analog lpf freq",equiv_analog_lowpass_freq
#print "D1,D2",D_1,D_2
#print "order=",N 
#print A_k 
#print "Poles",poles
print "Numerator\n",numer 
print "\nDenominator\n",denom


#plotting poles

# plt.figure(1)
# plt.grid(True)
# neg_poles=np.zeros([0],dtype='complex64')
# for c in poles:
#   if(c.real<= 0):
#     neg_poles=np.append(neg_poles,c)
# plt.scatter(neg_poles.real,neg_poles.imag,s=50,marker='x')
# plt.title('Pole Zero plot of Low pass filter')
# plt.ylabel('Imaginary')
# plt.xlabel('Real')

#nmrz=(numer.c).round(decimals=6)[::-1].real
#dmrz=(denom.c).round(decimals=6)[::-1].real
nmrz=(numer.c).round(decimals=6).real
dmrz=(denom.c).round(decimals=6).real


print "\nNormalized numerator Array:\n",nmrz/dmrz[0]
print "\nNormalized denominator Array:\n",dmrz/dmrz[0]
nyq_rate=f_sample/2


#Lattice Coefficients for IIR
An=dmrz/dmrz[0]
Bn=nmrz/dmrz[0]
N=len(dmrz)
Kn=np.zeros(N)
Cn=np.zeros(N)
for i in np.arange(N-1,-1,-1):
		Kn[i]=An[i]
		Cn[i]=Bn[i]
		An_tilda=An[::-1]
		if(Kn[i] != 1):
				An=(An-(Kn[i]*An_tilda))/(1-(Kn[i]**2))
		Bn=(Bn-Cn[i]*An_tilda)
		An=np.delete(An,len(An)-1)
		Bn=np.delete(Bn,len(Bn)-1)

print "\nLattice Coefficients:Kn\n",Kn[1:]
print "\nLattice Coefficients:Cn\n",Cn

# Plot Frequency response
plt.figure(2)
plt.clf()
plt.grid(True)
w,h= sg.freqz(nmrz,dmrz,worN=512)
plt.plot((w/np.pi)*nyq_rate, np.absolute(h), linewidth=2)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain')
plt.title('Frequency Response')
plt.ylim(-0.05, 1.3)

# Upper inset plot.
ax1 = plt.axes([0.44, 0.3, .45, .25])
plt.plot((w/np.pi)*nyq_rate, np.absolute(h), linewidth=2)
plt.xlim(16000.0,30000.0)
plt.ylim(0.85, 1.3)
plt.grid(True)

plt.figure(3)
plt.grid(True)
h_Phase = pl.unwrap(np.arctan2(np.imag(h),np.real(h)))
plt.plot(w/max(w),h_Phase)
plt.ylabel('Phase (radians)')
plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
plt.title(r'Phase response')
plt.show()