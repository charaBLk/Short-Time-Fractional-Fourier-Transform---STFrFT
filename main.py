import numpy as np 
import matplotlib.pyplot as plt
from cmath import *
from math import pi, sin, cos, atan

class STFRFT:
    
    def __init__(self, s, nb_q, alpha,delta_t, sigma=1):
        
        self.s = s
        self.N = len(s)
        self.nb_q = nb_q
        self.alpha = alpha
        self.delta_t = delta_t
        self.sigma = sigma
        self.delta_u = 2*np.pi * np.abs(np.sin(alpha)) / (self.N * delta_t)
    
    
    def K_alpha(self):
        """
        Computes the discrete kernel of the Fractional Fourier Transform (FrFT) of order alpha.
        """
        
        sin_a = np.sin(self.alpha)
        print(sin_a)
        cos_a = np.cos(self.alpha)
        cot_a = cos_a / sin_a
        delta_u = 2*np.pi * np.abs(sin_a) / (self.N * self.delta_t)

        A = np.sqrt((np.abs(sin_a) - 1j*np.sign(sin_a) * cos_a) / self.nb_q)
        n = np.arange(self.N)
        q = np.arange(self.nb_q)
        N_, Q_ = np.meshgrid(n, q, indexing='ij')
        
        res = np.exp(
            1j*2*np.pi*(
                ((N_ * self.delta_t )**2) *cot_a/2
                - np.sign(sin_a)*N_ * Q_ / self.nb_q
            )
        )
        return res
    
    def V_alpha(self):
        """
        Computes the Short-Time Fractional Fourier Transform (STFrFT) 
        in the form of a matrix V(n, q) with n ≤ N and q ≤ nb_q.

        Parameters
        ----------
        s : ndarray
            Input signal (complex or real), length N.
        sigma : float
            Standard deviation of the Gaussian window.
        alpha : float
            Fractional Fourier Transform order (in radians).
        nb_q : int
            Number of fractional frequency bins.
        delta_t : float
            Time step (1 / sampling frequency).

        Returns
        -------
        res : ndarray, shape (N, nb_q)
            STFrFT matrix.
        """

    

        
        t_win = np.linspace(-0.5, 0.5, 64)
        g = np.exp(-np.pi * (t_win / self.sigma)**2)

        L = len(g)
        k = (L - 1) // 2
        # Initialisation of the result matrix
        res = np.zeros((self.N, self.nb_q), dtype=np.complex128)
        
        # compute the kernel coefficients
        coef = self.K_alpha()  # shape (N, nb_q)
        
        #loop over time indices n
        for n in range(self.N):
            # boundries to avoid window overflow
            m_start = max(-k, -n)
            m_end   = min(k, self.N - n - 1)
            
            # indexes for signal and window
            signal_idx = np.arange(n + m_start, n + m_end + 1)
            window_idx = np.arange(m_start + k, m_end + k + 1)
            
            # loop over fractional frequency bins q
            for q in range(self.nb_q):
                res[n, q] = np.sum(self.s[signal_idx] * g[window_idx] * coef[signal_idx, q])
        
        # Normalization
        res /= np.sum(np.abs(g))
        
        return res
    
    
    
    