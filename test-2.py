from main import STFRFT
import numpy as np 
import matplotlib.pyplot as plt
from cmath import *
from math import pi, sin

fs = 80
T = 10
N = fs * T
t = np.linspace(0, T, T * fs)
nb_q = 150
delta_t = t[1] - t[0]

# Complex signal
phi1 = 30*np.pi*t + 5*np.pi*t**2 + (2*np.pi/9)*t**3 - (np.pi/30)*t**4
phi2 = 40*np.pi*t + 5*np.pi*t**2 + (2*np.pi/9)*t**3 - (np.pi/30)*t**4
s = 5*np.exp(1j*phi1) + 5*np.exp(1j*phi2)

sigma=0.55


alpha = -0.199
#alpha = pi/2 
alpha = 0.002

a = STFRFT(s, nb_q, alpha, delta_t, sigma)
res = a.V_alpha()

time_axis = t
freq_axis = (np.arange(nb_q)) / nb_q * fs * sin(alpha)

plt.figure(figsize=(12, 8)) 
plt.pcolormesh(time_axis, freq_axis, np.abs(res.T), shading='auto', cmap='jet')
plt.xlabel('time [s]')
plt.ylabel('fractional frequency [Hz]')
plt.title(f'STFrFT (α = {alpha})')
plt.colorbar(label='Amplitude')
plt.tight_layout(pad=3.0)  
plt.show()