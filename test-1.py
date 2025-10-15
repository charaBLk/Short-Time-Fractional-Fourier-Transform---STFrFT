from main import STFRFT
import numpy as np 
import matplotlib.pyplot as plt
from cmath import *
from math import pi, sin, cos, atan


fs = 500
T = 1
t = np.linspace(0, T, T * fs)
s4 = np.exp(1j * 2 * np.pi * (250 * t -  250/2*t**2))  # chirp
nb_q = 300
delta_t = t[1] - t[0]
sigma=0.6

# arccot(250) = 0.004
alpha = pi/2

a = STFRFT(s4, nb_q, alpha, delta_t, sigma)
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