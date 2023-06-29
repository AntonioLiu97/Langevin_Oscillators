import numpy as np
from scipy.integrate import odeint
from scipy.fftpack import fft, fftfreq
from scipy import signal

class Harmonic_Oscillator:

    def __init__(self, F, gamma, m = 1, k = 1, sigma_x = 0, sigma_v = 0, c_drive = 1):
        '''
        Defines a driven, damped harmonic oscillator:
            F(t) - kx - gamma dx/dt = d^2x/dt^2 * m
        
        To add noise, we write the second order ODE as two first order ODE,
        and add Brownian noise to x and v separately
        dv = F(t)dt - kx(t)/m dt - gamma v(t)/m dt + dW_v
        dx = v(t)dt + dW_x
            
        k: spring constant
        F: driving function
        gamma: damping constant
        m: oscillator mass
        sigma_x: additive noise to x (Brownian bombardment)
        sigma_v: additive noise to v (noise in force)
        c_drive: drive coupling strenght (linearly scales F)
        '''        
        self.F = F
        self.gamma = gamma
        self.m = m
        self.k = k
        self.sigma_x = sigma_x
        self.sigma_v = sigma_v
        self.c_drive = c_drive
        
    def _equations(self, y, t):
        """
        Express the second order ODE as two first order ODEs
        
        Parameters:
        y: List [x, v] where x is the displacement and v is the velocity.
        t: Time.
        
        Returns:
        List of the first derivatives [dx/dt, dv/dt].
        """
        x, v = y
        dxdt = v 
        dvdt = self.c_drive * self.F(t) / self.m - self.k * x / self.m - self.gamma * v / self.m
        return [dxdt, dvdt]        
    
    def run_without_noise(self, initial_condition, T, dt):
        """
        inputs
            initial_condition: tuple (initial position, velocity)
            T: total time to integrate over
            dt: time step
        
        output: 
            x_array: a numpy array containing the displacement of the oscillator at each time step
        """
        t = np.arange(0, T, dt)
        solution = odeint(self._equations, initial_condition, t)
        x_array = solution[:, 0]
        return x_array, t    
    
    def run_with_noise(self, initial_condition, T, dt):
        """
        Use first order method to intergrate with noise
        """
        t = np.arange(0, T, dt)
        x_array = np.zeros_like(t)
        v_array = np.zeros_like(t)
        x_array[0], v_array[0] = initial_condition

        for i in range(1, t.shape[0]):
            x_array[i] = x_array[i-1] + dt * v_array[i-1] \
                        + self.sigma_x * np.sqrt(dt) * np.random.normal()
            v_array[i] = v_array[i-1] + dt * (self.c_drive * self.F(t[i-1]) / self.m - self.k * x_array[i-1] / self.m - self.gamma * v_array[i-1] / self.m) \
                        + self.sigma_v * np.sqrt(dt) * np.random.normal()

        return x_array, t
    
    
class Radial_Oscillator:

    def __init__(self, F, alpha = 0.1, 𝜔 = 1, 𝜇 = 1, k = 0.3, sigma_r = 0, sigma_𝜙 = 0):
        '''
        Defines an oscillator in polar coodinate driven in phase,
        with noise directly added in polar-coordinates
            dr = alpha (𝜇 r dt - r^3 dt + dw_r)
            d𝜙 = 𝜔 dt + k sin(F(t) - 𝜙) dt + dw_𝜙
            
        k: coupling to drive
        𝜔: intrinsic frequency
        𝜇: Hopf control parameter (criticality at 0)
        sigma(s): noise strength        
        '''        
        self.F = F
        self.k = k
        self.𝜇 = 𝜇
        self.𝜔 = 𝜔
        self.sigma_r = sigma_r
        self.sigma_𝜙 = sigma_𝜙
        self.alpha = alpha
        
    def _equations(self, Y, t):
        r, 𝜙 = Y
        drdt = self.alpha * (self.𝜇 * r - r**3)
        d𝜙dt = self.𝜔 + self.k * np.sin(self.F(t) - 𝜙)
        return [drdt, d𝜙dt]    
    
    def run_without_noise(self, initial_condition, T, dt):
        t = np.arange(0, T, dt)
        Y = odeint(self._equations, initial_condition, t)
        r_array, 𝜙_array = Y[:, 0], Y[:, 1]
        x_array = r_array * np.sin(𝜙_array)
        return x_array, t    
        
    def run_with_noise(self, initial_condition, T, dt):
        """
        Use first order method to intergrate with noise
        """
        t = np.arange(0, T, dt)
        r_array = np.zeros_like(t)
        𝜙_array = np.zeros_like(t)
        r_array[0], 𝜙_array[0] = initial_condition
        for i in range(1, t.shape[0]):
            r_array[i] = r_array[i-1] + self.alpha * (self.𝜇 * r_array[i-1] * dt - r_array[i-1]**3 * dt) \
                        + self.sigma_r * np.sqrt(dt) * np.random.normal()
            𝜙_array[i] = 𝜙_array[i-1] + self.𝜔 * dt + dt * self.k * np.sin(self.F(t[i-1]) - 𝜙_array[i-1])  \
                        + self.sigma_𝜙 * np.sqrt(dt) * np.random.normal()

        x_array = r_array * np.sin(𝜙_array)
        return x_array, t
    
    
