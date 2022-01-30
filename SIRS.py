# Forth order runge kutta

class SIRS:

	def __init__(self,S_0,I_0,R_0):

		"""
		Parameters
        ----------
        S0: Initial number of Suceptible people
        I0:               ... Infected
        R0:               ... Recovered people
        a:  Rate of Transmission                    [1/Time]
        b:      ... Recovery                        [1/Time]
        c:      ... Immunity loss                   [1/Time]
        """

		#Inital population
		self.S_0 = S_0
		self.I_0 = I_0
		self.R_0 = R_0
		self.N = self.S_0 + self.I_0 + self.R_0
		
		#Rate of; Transmission, Recovery and Immunity
		self.a = 4.0
		self.b = 1.0
		self.c = 0.5

		self.Nt = 10
		self.N_len = 10000
		self.h = float(self.Nt)/self.N_len
		self.t = np.linspace(0, self.Nt, self.N_len)

	def eulersmethod(self):

		#Array allocation
		S = np.zeros(self.N_len,'float'); I=np.zeros(self.N_len,'float'); R=np.zeros(self.N_len,'float')
		S[0] = self.S_0
		I[0] = self.I_0
		R[0] = self.R_0

		#Eulers method
		for i in range(0,len(S)-1,1):

			S[i+1] = S[i] + self.h*(self.c*R[i] - (self.a*S[i]*I[i])/self.N)
			I[i+1] = I[i] + self.h*((self.a*S[i]*I[i])/self.N -self.b*I[i])
			R[i+1] = R[i] + self.h*(self.b*I[i] -self.c*R[i])

		import matplotlib.pylab as plt

		plt.figure(1)
		plt.plot(self.t,S)
		plt.plot(self.t,I)
		plt.plot(self.t,R)
		plt.xlabel('Time')
		plt.ylabel('Number of individuals in respective population')
		plt.legend(['Suceptible','Infected','Recovered'])
		plt.show()

	def analytical(self):
		s_star = b/a
		i_star = (1-(b/a))/(1+(b/c))
		r_star = (b/c)*(1-(b/a))/(1+(b/c))

		return s_star,i_star,r_star
	

	def MonteCarlo(self):

		S = self.S_0;	I = self.I_0;	R = self.R_0
		S_array = np.zeros(self.N_len,'float'); I_array =np.zeros(self.N_len,'float'); R_array=np.zeros(self.N_len,'float')
		S_array[0] = self.S_0; I_array[0] = self.I_0; R_array[0] = self.R_0;
	
		#Choose time step.
		import random
		delta_t = np.min([4/(self.a*self.N),1/(self.b*self.N),1/(self.c*self.N)])
		for i in range(0,self.N_len-1,1):
			tmp = random.random()
		
			if tmp < (self.a*S*I*delta_t)/self.N:
				S = S-1
				I = I+1

			if tmp < (self.b*I*delta_t):
				I = I-1
				R = R+1

			if tmp < (self.c*R*delta_t):
				R = R-1
				S = S+1
	
			S_array[i+1] = S
			I_array[i+1] = I
			R_array[i+1] = R
		
		import matplotlib.pylab as plt

		plt.figure(1)
		plt.plot(self.t,S_array)
		plt.plot(self.t,I_array)
		plt.plot(self.t,R_array)
		plt.ylabel('Number of individuals in respective population')
		plt.legend(['Suceptible','Infected','Recovered'])
		plt.show()


if __name__ == '__main__':
	import numpy as np
	s = SIRS(S_0=300,I_0=100,R_0=0)
	s.eulersmethod()
	#s.MonteCarlo()
