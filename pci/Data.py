import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from scipy.optimize import curve_fit
import random
import sys

class AsymmetricMeasurment:
	"""
	A class for dealing with asymmetric probability distribution of parameters given by their confidence intervals.

	e.g. Value = 10 ( -1, +2 ) (68% confidence interval)
	"""
	def __init__(self, me=10.0, err_n=1.0, err_p=1.0, N=10000, confidence = 68, precision = 500, creation_method='by_parameters', data=[]):
		"""
		Instantiates the class AsymmetricMeasurment.
	   
		Parameters
		----------
		me (measurment) : float
			The measured value, the mode of the distribution.
		err_n : float
			The negative error, the minus error (insert only the value, without the (-) sign).
		err_p : float
			The positive error, the plus error (insert only the value, without the (+) sign).
		N : int
			The number of points to constructed the distribution.
		confidence : int
			The statistical significance of the errors given by the percentile (%) confience intervals.
			The implemented options are 68, 90 and 95:
				Confidence(%)	Δχ2	ΔlogL
				68%	1.00	-0.50
				90%	2.71	-1.36
				95%	4.00	-2.00
				97%	4.60	-2.30

		Returns
		-------
		Nothing.
		"""		


	   
		if (err_n*err_p) < 0:
			raise RuntimeError("Both 'err_n' and 'err_p' have to be positive")

		try:
			self._get_const(confidence) 
		except KeyError:
			print('confidence interval value needs to be chose among 68, 90, 95. Assuming 68!.')
			confidence = 68

		if N < 1000:
			raise RuntimeError("N needs to be greater than 1000!")




		self.me = me
		self.err_n = err_n
		self.err_p = err_p
		self.conf_inter = confidence
		self.N = int(N)
		self.creation_method = (creation_method)
		self.bins = int(precision)

		if not any(data):
			self.data = np.asarray([])
		else:
			self.data = np.asarray(data)

		if(str(self.creation_method) == 'by_parameters'):
			delta = 8*np.max([self.err_n,self.err_p])/self._get_const(self.conf_inter)
			self.x_lim =  [self.me - delta, self.me + delta]
			self.x_values = np.linspace(self.x_lim[0], self.x_lim[1], self.N)
			

			self.norm = 1.0
			self.norm = self._calculate_norm()

			
			self.pdf_values = np.asarray(self.pdf(self.x_values))
			self.cdf_values = np.asarray(self._calculate_cdf_values())
			self.log_likelihood_values = np.asarray(self.log_likelihood(self.x_values))
			
			self._generate_data()

		elif str(self.creation_method) == 'by_function':
			self.N = self.data.size

			self._fit()
			#self.sigma_n, self.sigma_p = self.estimate()

			delta = 8*np.max([self.err_n,self.err_p])/self._get_const(self.conf_inter)
			self.x_lim =  [self.me - delta, self.me + delta]
			self.x_values = np.linspace(self.x_lim[0], self.x_lim[1], self.N)

			self.norm = 1.0
			self.norm = self._calculate_norm()

			
			self.pdf_values = np.asarray(self.pdf(self.x_values))
			self.cdf_values = np.asarray(self._calculate_cdf_values())
			self.log_likelihood_values = np.asarray(self.log_likelihood(self.x_values))
			

	def __str__(self):
		output = "Value = {:.2e} ( - {:.2e} , + {:.2e} )\n({:.0f}% confidence level)"
		return output.format(self.me, self.err_n, self.err_p, self.conf_inter)



	
	def _get_const(self,conf_inter ):
		c_dic = {
				'68': 1,
				'90': 1.64,
				'95': 2,
				'97': 2.14}
		const = c_dic[str(conf_inter)]
		return const

	def _integrate(self):
		delta_x = self.x_lim[1] - self.x_lim[0]
		c = delta_x / (self.N - 1)
		# x_values = np.linspace(self.x_limits[0], self.x_limits[1], self.N, dtype=float)
		area = np.sum(c * self.pdf(self.x_values))
		return area

	def _calculate_norm(self):
		area = self._integrate()
		norm = 1/area
		return norm

	def pdf(self, x):
		"""
		Measures the Probability Density Function (CDF) for a given (x) value.
        
        Returns
        --------
        pdf_x : float
            Value of PDF at x.
        """
		c = self._get_const(self.conf_inter)
		
		par_1 = (2.0 * self.err_p * self.err_n) / ((self.err_p + self.err_n)*c)
		par_2 = (self.err_p  - self.err_n) / ((self.err_p  + self.err_n)*c)
		par_3= (-1.0/2.0) * ((self.me - x)/(par_1 + par_2*(x - self.me)))**2.0
		
		par_4 = self.norm / (2.0 * np.pi)**0.5
		value = par_4 * np.exp(par_3)
		#print(v1, v2)
		return value
	
	def log_likelihood(self, x):
		"""
		Measures the logarithm of the likelihood (Variable Width Gaussian, from R. Barlow's 2004 paper "Asymmetric Statistical Errors") for a given (x) value.
        
        Returns
        --------
        log_L_x : float
            Value of log_likelihood at x.
        """

		c = self._get_const(self.conf_inter)
		par_1 = (2.0 * self.err_p * self.err_n) / ((self.err_p + self.err_n)*c)
		par_2 = (self.err_p  - self.err_n) / ((self.err_p  + self.err_n)*c)
		value = (-1.0/2.0) * ((self.me - x)/(par_1 + par_2*(x - self.me)))**2.0
		return value

	def _calculate_cdf_values(self):
		
		cdf_values = np.asarray([])
		for i in range(self.N):
			area = np.trapz(self.pdf_values[:i], x = self.x_values[:i])
			cdf_values = np.append(cdf_values, area)
		return cdf_values

	def cdf(self, x):
		"""
		Measures the Cumulative Density Function (CDF) for a given (x) value.
        
        Returns
        --------
        cdf_x : float
            Value of CDF at x.
        """
		cdf = interpolate.interp1d(self.x_values, self.cdf_values, kind='nearest')
		return cdf(x)

	def _calculate_inverse_cdf(self, x):
		inverse_cdf = interpolate.interp1d(self.cdf_values, self.x_values, kind='nearest')
		return inverse_cdf(x)


	def generate_random(self):
		"""
		Generates a random number using the probability density function of the object (see .plot_pdf() function).
        
        Returns
        --------
        n : float
            Random number.
        """
		rnd_prob = np.random.uniform(self.cdf_values[0], self.cdf_values[-1])
		n_rand = self._calculate_inverse_cdf(rnd_prob)
		return n_rand

	def _generate_data(self):
		rnd_prob = np.random.uniform(0, 1, self.N)
		self.data = self._calculate_inverse_cdf(rnd_prob)

	def get_confidence(self, conf_inter=68):
		"""
		Gets the confidence intervals. Implemented options are 68, 90 and 95% (e.g conf_inter = 90)
        
        Returns
        --------
        conf_inter : list
            2-size list with the x values at the selected confidence interval.
        """


		d_ll = self._get_const(conf_inter)**2 / -2
		
		flag_up = self.x_values > self.me
		
		conf_up =(self.x_values[flag_up])[np.where(abs(self.log_likelihood_values[flag_up] - d_ll) == np.min(abs(self.log_likelihood_values[flag_up]  - d_ll)))]
		flag_low = self.x_values < self.me
		conf_low =(self.x_values[flag_low])[np.where(abs(self.log_likelihood_values[flag_low] - d_ll) == np.min(abs(self.log_likelihood_values[flag_low]  - d_ll)))]
		return [conf_low, conf_up]

	@staticmethod
	def fit_func(x, norm, me, err_n, err_p, c):
		
		par_1 = (2.0 * err_p * err_n) / ((err_p + err_n)*c)
		par_2 = (err_p - err_n) / ((err_p + err_n)*c)
		par_3 = (-1.0 / 2.0) * ((me - x) / (par_1 + par_2 * (x - me))) ** 2.0
		par_4 = norm / (2.0 * np.pi) ** 0.5
		value = par_4 * np.exp(par_3)
		
		return value

	def _fit(self, expected_values=None):
		y, x, _ = plt.hist(self.data, bins=int(self.bins))
		plt.clf()
		plt.close('all')
		x = (x[1:] + x[:-1]) / 2  # for len(x)==len(y)
		
		mod = None
		max_y = max(y)
		for i in range(len(y)):
			if y[i] == max_y:
				mod = x[i]

		
		min_data = min(self.data)
		max_data = max(self.data)
		norm = 1000.0
		c = self._get_const(self.conf_inter)
		if not expected_values:
			expected_values = norm, mod, (mod - min_data) * 0.1, (max_data - mod) * 0.1, c

		expected = (expected_values[0], expected_values[1], expected_values[2], expected_values[3], expected_values[4])
		
		bounds = ([0, -np.inf, 0, 0, c], [np.inf, np.inf, np.inf, np.inf, c+1e-10])
		params, cov = curve_fit(self.fit_func, x, y, expected,  bounds = bounds, method='trf')
		self.norm = params[0]
		self.me = params[1]
		#print("params", params)
		if params[2] > 0.0:
			self.err_n = (params[2])
			self.err_p = (params[3])
		else:
			self.err_n = (params[3])
			self.err_p = (params[2])

	def plot_pdf(self, show=True, save=False):

		plt.clf()
		plt.plot(self.x_values, self.pdf_values, color="blue", label = r'$%f^{+%f}_{-%f}$  (%d%% confidence level)' % (self.me, self.err_p,self.err_n, self.conf_inter))
		ymin,ymax = plt.ylim()
		plt.ylim(ymin, 1.2*ymax)
		plt.vlines(self.me, ymin =0,ymax=ymax,color = 'black')
		plt.vlines(self.me + self.err_p, ymin =0,ymax=ymax, color = 'black', linestyle = '--')
		plt.vlines(self.me - self.err_n, ymin =0,ymax=ymax, color = 'black',  linestyle = '--')
		plt.legend(fontsize = 12)
		plt.xlabel("x", fontsize=12)
		plt.ylabel("Normalized Probability", fontsize=12)

		if save:
			plt.savefig("plot_pdf.png", dpi=300)

		if show:
			plt.show()

	def plot_log_likelihood(self, show=True, save=False):
		plt.clf()
		plt.plot(self.x_values, self.log_likelihood(self.x_values))
		plt.ylim([-5, 0.5])
		

		plt.xlabel("x", fontsize=14)
		plt.ylabel(r'ln $\Delta$ L', fontsize=14)

		plt.axhline(y=self._get_const(self.conf_inter)**2 / -2, color="black", ls="--", label = r"ln  $\Delta$L = "+str(((self._get_const(self.conf_inter))**2)/-2))
		plt.axhline(y=0, color="black", label = r"ln  $\Delta$L = 0")
		plt.legend()
		if save:
			plt.savefig("plot_log_likelihood.png", dpi=300)

		if show:
			plt.show()

	def plot_cdf(self, show=True, save=False):
		plt.plot(self.x_values, self.cdf(self.x_values))

		if save:
			plt.savefig("plot_cdf.png", dpi=300)

		if show:
			plt.show()

	def plot_data(self, bins=None, show=True, save=False):
		if not bins:
			bins = self.bins

		plt.clf()
		plt.hist(self.data, bins=bins, density=True, color="green", alpha=0.7)

		if save:
			plt.savefig("plot_data.png", dpi=300)

		if show:
			plt.show()

	def plot_data_and_pdf(self, bins=None, show=True, save=False):
		if not bins:
			bins = self.bins

		plt.clf()
		plt.hist(self.data, bins=bins, density=True, color="green", alpha=0.6)
		plt.plot(self.x_values, self.pdf_values, color="blue")

		if save:
			plt.savefig("plot_data_and_pdf.png", dpi=300)

		if show:
			plt.show()

	def __add__(self, other):
		if isinstance(other, self.__class__):
			if(other.data.size != self.data.size):
				print('The size (N) of the data of both AsymmetricMeasurment have to be equal')
				sys.exit()
			add = self.data + other.data
			mod = self.me + other.me
			#print(len(add))
		elif isinstance(other, (int, float)):
			add = self.data + float(other)
			mod = self.me + float(others)
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		temp_obj = AsymmetricMeasurment(creation_method='by_function', data=add, confidence =self.conf_inter)
		return temp_obj

	def __radd__(self, other):
		if isinstance(other, self.__class__):
			if(other.data.size != self.data.size):
				print('The size (N) of the data of both AsymmetricMeasurment have to be equal')
				sys.exit()
			add = other.data + self.data
			mod = other.me + self.me
		elif isinstance(other, (int, float)):
			add = float(other) + self.data
			mod = float(othr) + self.me
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		temp_obj = AsymmetricMeasurment(creation_method='by_function', data=add, confidence =self.conf_inter)
		return temp_obj

	def __sub__(self, other):
		if isinstance(other, self.__class__):
			if(other.data.size != self.data.size):
				print('The size (N) of the data of both AsymmetricMeasurment have to be equal')
				sys.exit()
			add = self.data - other.data
			mod - self.me - other.me
		elif isinstance(other, (int, float)):
			add = self.data - float(other)
			mod = self.me - float(other)
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		temp_obj = AsymmetricMeasurment(creation_method='by_function', data=add, N=self.N, confidence =self.conf_inter)
		return temp_obj

	def __rsub__(self, other):
		if isinstance(other, self.__class__):
			if(other.data.size != self.data.size):
				print('The size (N) of the data of both AsymmetricMeasurment have to be equal')
				sys.exit()
			add = other.data - self.data
			mod = other.data - self.me
		elif isinstance(other, (int, float)):
			add = float(other) - self.data
			mod = float(other) - self.me
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		temp_obj = AsymmetricMeasurment(creation_method='by_function', data=add, confidence =self.conf_inter)
		return temp_obj

	def __mul__(self, other):
		if isinstance(other, self.__class__):
			if(other.data.size != self.data.size):
				print('The size (N) of the data of both AsymmetricMeasurment have to be equal')
				sys.exit()
			add = self.data * other.data
			mod = self.me * other.me
		elif isinstance(other, (int, float)):
			add = self.data * float(other)
			mod = self.me * float(other)
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		temp_obj = AsymmetricMeasurment(creation_method='by_function', data=add, confidence =self.conf_inter)
		return temp_obj

	def __rmul__(self, other):
		if isinstance(other, self.__class__):
			if(other.data.size != self.data.size):
				print('The size (N) of the data of both AsymmetricMeasurment have to be equal')
				sys.exit()
			add = other.data * self.data
			mod = other.me * self.me
		elif isinstance(other, (int, float)):
			add = float(other) * self.data
			mod = float(other) * self.me
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		temp_obj = AsymmetricMeasurment(creation_method='by_function', data=add, confidence =self.conf_inter)
		return temp_obj

	def __truediv__(self, other):
		if isinstance(other, self.__class__):
			if(other.data.size != self.data.size):
				print('The size (N) of the data of both AsymmetricMeasurment have to be equal')
				sys.exit()
			add = self.data / other.data
			mod = self.me / other.me
		elif isinstance(other, (int, float)):
			add = self.data / float(other)
			mod = self.me / float(other)
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()

		
		temp_obj = AsymmetricMeasurment(creation_method='by_function', data=add, confidence =self.conf_inter)
		return temp_obj

	def __rtruediv__(self, other):
		if isinstance(other, self.__class__):
			if(other.data.size != self.data.size):
				print('The size (N) of the data of both AsymmetricMeasurment have to be equal')
				sys.exit()
			add = other.data / self.data
			mod = other.me / self.me
		elif isinstance(other, (int, float)):
			add = float(other) / self.data
			mod = float(other) / self.me
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		
		temp_obj = AsymmetricMeasurment(creation_method='by_function', data=add, confidence =self.conf_inter)
		return temp_obj

	def __pow__(self, other):
		if isinstance(other, self.__class__):
			if(other.data.size != self.data.size):
				print('The size (N) of the data of both AsymmetricMeasurment have to be equal')
				sys.exit()
			add = self.data ** other.data
			mod = self.me ** other.me
		elif isinstance(other, (int, float)):
			add = self.data ** float(other)
			mod = self.me ** float(other)
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()

		temp_obj = AsymmetricMeasurment(creation_method='by_function', data=add, confidence =self.conf_inter)
		return temp_obj

	def __rpow__(self, other):
		if isinstance(other, self.__class__):
			if(other.data.size != self.data.size):
				print('The size (N) of the data of both AsymmetricMeasurment have to be equal')
				sys.exit()
			add = other.data ** self.data
			mod = other.me ** self.me
		elif isinstance(other, (int, float)):
			add = float(other) ** self.data
			mod = float(other) ** self.me
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()

		
		temp_obj = AsymmetricMeasurment(creation_method='by_function', data=add, confidence =self.conf_inter)
		return temp_obj

	@staticmethod
	def Propagator(func, asymmetric_measurments= [], pars = [], N = 10000, confidence = 68, precision = 300):
		'''
		A function to propagate errors of asymmetric measurment in user-defined functions,
		it uses Monte Carlo simulations and the probability distribution functions of the AsymmetricMeasurment objects.

		Parameters
		----------
		func (function) : "<class 'function'>"
			The function in which the errors will be propagated needs to have a very specific format, 
			the AsymmetricMeasurment parameters need to be given as the first parameter as a (one) list,
			the other (float-like) parameters have to be the second parameter
			and be given as another (one) list. The function needs to return only 
			one parameter. See an example:

			def func(asymmetric_measurments: list, parameters: list):
				as1, as2, as3 = asymmetric_measurments
				p1,p2,p3 = parameters
				value = as1 * as2 * as3 -p1 * p2 * p3
				return value 

		asymmetric_measurments: list
			The list containing the asymmetric measurments (instances of the AsymmetricMeasurment class),

		pars: list
			The list containing the float(or int)-like parameters.

		N: int
			Number of Monte Carlo iterations.

		confidence : int
			The statistical significance of the errors given by the percentile (%) confience intervals.
			The implemented options are 68, 90 and 95.
			

		precision: float
			Number of bins for the histogram fit.
			Standart value is 300, for complex functions it may be intersting to increase (e.g 500, 1000, 5000?),
			however it will take longer to fit the pdf.



		Returns
		-------
		Measurment: __main__.AsymmetricMeasurment
			An instance of the AsymmetricMeasurment class resulted from the propagation of errors.
		'''		


		if str(type(func)) != "<class 'function'>":
			raise RuntimeError("The 'func' parameter needs to be a function ({})".format(str(type(func))))

		if len(asymmetric_measurments) == 0:
			raise RuntimeError("The 'asymmetric_measurments' list cannot be empty")

		if str(type(asymmetric_measurments)) != "<class 'list'>":
			raise RuntimeError("The 'asymmetric_measurments' parameters needs to be a list ({})".format(str(type(asymmetric_measurments))))

		for i in range(len(asymmetric_measurments)):
			if str(type(asymmetric_measurments[i])) != "<class 'pci.Data.AsymmetricMeasurment'>":
				raise RuntimeError("The elements of the 'asymmetric_measurments' need to be  AsymmetricMeasurments class objects ({})".format(str(type(asymmetric_measurments[i]))))
			n = (asymmetric_measurments[i]).N
			if (n != (asymmetric_measurments[i-1]).N):
				raise RuntimeError('The size (N) of the data in all AsymmetricMeasurments objects have to be equal')

		if (asymmetric_measurments[0].N != N):
			N = asymmetric_measurments[0].N

		as_ = [j.data for j in asymmetric_measurments]
		new_data = func(as_, pars)

		temp_obj = AsymmetricMeasurment(creation_method='by_function', data=new_data, confidence =confidence, precision=precision)
		return temp_obj




if __name__ == "__main__":
	a=AsymmetricMeasurment(10,1,1, confidence=68)