"""
Analyse and plot disc fragment properties.
"""
import numpy as np
import matplotlib.pyplot as plt
from libanalysis import PhantomAnalysis as pa

class Analyser():

	def __init__(self, idisc, wd):
		if idisc < 10:
			disc = 'disc00%i' %idisc
		elif idisc < 100:
			disc = 'disc0%i' %idisc
		else:
			disc = 'disc%i' %idisc
		self.file = disc 								#disc filename
		self.idisc = idisc								#ID of disc
		self.wd = wd									#directory where data is stored
		self.disc = pa('%s/%s' %(self.wd,self.file))
		self.readClumpcat()								#read in clump info from clumpcat file
		self.annulus = 3								#annulus in AU to analyse around clump centre

		#define constants
		self.udens = self.disc.units['udens'] #gcm-3
		self.umass = self.disc.units['umass'] #g
		self.udist = self.disc.units['udist'] #cm
		self.utime = self.disc.units['utime']
		self.uenerg = self.umass*self.udist**2/self.utime**2


	def readClumpcat(self):
		'''
		Read in fragment properties from clumpcat file
		'''
		f = open('%s/raw_clumpcat_%s' %(self.wd,self.file))
		data = []
		for line in f.readlines()[1:]:
			data.append(map(float, line.split()))
		self.clumpcat = data

	def azAverage(self,rads,vals,nbins=50):
		'''
		Azimuthally average values, centred on clump center
		'''
		try:
			avVals = []
			bins = np.linspace(0,self.annulus,nbins)
			for i, bin in enumerate(bins[:-1]):
				av = np.max(vals[(rads>bins[i]) & (rads<=bins[i+1])])
				avVals.append(av)
		except:
			#if bin size is too small, and some bins have no particles, make bins bigger
			nbins=25			
			avVals = []
			bins = np.linspace(0,self.annulus,nbins)
			for i, bin in enumerate(bins[:-1]):
				try:
					av = np.max(vals[(rads>bins[i]) & (rads<=bins[i+1])])
				except:
					av = 0
				avVals.append(av)
		return bins[:-1], avVals

	def fitPolytrope(self,rads,vals):
		'''
		fit polytropic line to plotted data
		'''
		coeffs = np.polyfit(rads, vals, deg=2)
		#polyvals = coeffs[0]*rads**3 + coeffs[1]*rads**2 + coeffs[2]*rads + coeffs[3]
		polyvals = coeffs[0]*rads**2 + coeffs[1]*rads + coeffs[2]
		return rads, polyvals

	@staticmethod
	def polytropicProfile(n=[1.5]):
		dx = 0.001
		xs = []
		thetas = []
		for j, index in enumerate(n):
			theta = [1]
			dtheta_dx = [0]
			x = [1E-6]
			i=0
			while theta[i] > 0 and i<5000:
				dtheta_dx.append(dtheta_dx[i] - (2*dtheta_dx[i]/x[i] + theta[i]**index)*dx)
				theta.append(theta[i] + dtheta_dx[i]*dx)
				x.append(x[i] + dx)
				i+=1
			xs.append(np.array(x))
			thetas.append(np.array(theta))
		return xs, thetas

	def plotClump(self,ID):
		'''
		Plot radial profiles of clump properties - single clump (dens, temp, dustfrac).
		'''
		clumpxyz = self.clumpcat[ID][:3]
		r2 = (self.disc.xyzh[0]-clumpxyz[0])**2 + (self.disc.xyzh[1]-clumpxyz[1])**2
		members = np.sqrt(r2) < self.annulus 			#members are all particles within radial annulus

		gas = self.disc.itype == 1
		dust = self.disc.itype == 2
		dustfrac = self.disc.dustfrac*1e8  #as I originally set dust-to-gas=1e-10

		#Calculate temperatures from thermal energies
		k = 1.38064852e-16    #ergs
		mH = 1.6735575e-24    #grams
		gmw = 2.381           #mean mass taken from Phantom
		N = sum(gas*self.disc.massofgas)*self.umass/mH/gmw #number of atoms
		temp = 2.*self.disc.utherm*self.uenerg/3./k/N

		radplot, densav = self.azAverage(np.sqrt(r2[members&gas]), self.disc.density[members&gas]*self.udens)
		_, tempav = self.azAverage(np.sqrt(r2[members&gas]), temp[members&gas])
		raddust, fdustav = self.azAverage(np.sqrt(r2[members&gas]), dustfrac[members&gas],nbins=20)


		_, denspoly = self.fitPolytrope(radplot,densav)
		_, temppoly = self.fitPolytrope(radplot,tempav)
		raddust, fdustpoly = self.fitPolytrope(raddust,fdustav)

		polyx, polytheta = Analyser.polytropicProfile(n=1.5)

		#Plot density vs. radius
		plt.figure(1)
		plt.scatter(np.sqrt(r2[members & gas]), self.disc.density[members & gas]*self.udens, s=0.1)
		plt.plot(radplot, densav)
		plt.plot(radplot, denspoly)
		plt.ylim([0, max(self.disc.density[members&gas]*self.udens*1.1)])
		plt.ylabel(r'$\rho$ (g cm$^{-3}$)')
		plt.xlabel('Dist from clump centre (au)')
		plt.title('%s : Gas density vs. radius' %self.file)
		plt.legend(['%.1f' %self.disc.time])

		#Plot temperature vs. radius
		plt.figure(2)
		plt.scatter(np.sqrt(r2[members & gas]), temp[members & gas], s=0.1)
		plt.plot(radplot, tempav)
		plt.plot(radplot, temppoly)
		plt.ylim([0, max(temp[members&gas]*1.1)])
		plt.ylabel('Temp (K)')
		plt.xlabel('Dist from clump centre (au)')
		plt.title('%s : Gas temp vs. radius' %self.file)
		plt.legend(['%.1f' %self.disc.time])

		#Plot dust fraction vs. radius
		plt.figure(3)
		plt.scatter(np.sqrt(r2[members & gas]), dustfrac[members & gas], s=0.1)
		plt.plot(raddust, fdustav)
		plt.plot(raddust, fdustpoly)
		plt.ylim([0, max(dustfrac[members&gas]*1.1)])
		plt.ylabel('Dust-to-gas ratio')
		plt.xlabel('Dist from clump centre (au)')
		plt.title('%s : Dust fraction vs. radius' %self.file)
		plt.legend(['%.1f' %self.disc.time])

		plt.show()

	@staticmethod
	def plotMultiClumps(discs,IDs):
		'''
		Plot radial profiles of clump properties (dens, temp, dustfrac).
		'''

		#Initially plot polytropic indices n=1.5 for all clumps
		polyx, polytheta = Analyser.polytropicProfile(n=np.ones(len(IDs))*1.5)

		keepplotting = True
		while keepplotting:
			timelegend = []
			for i, [d, ID] in enumerate(zip(discs,IDs)):
				clumpxyz = d.clumpcat[ID][:3]
				r2 = (d.disc.xyzh[0]-clumpxyz[0])**2 + (d.disc.xyzh[1]-clumpxyz[1])**2
				members = np.sqrt(r2) < d.annulus 			#members are all particles within radial annulus

				gas = d.disc.itype == 1
				dust = d.disc.itype == 2
				dustfrac = d.disc.dustfrac*1e8  #as I originally set dust-to-gas=1e-10

				#Calculate temperatures from thermal energies
				k = 1.38064852e-16    #ergs
				mH = 1.6735575e-24    #grams
				gmw = 2.381           #mean mass taken from Phantom
				N = sum(gas*d.disc.massofgas)*d.umass/mH/gmw #number of atoms
				temp = 2.*d.disc.utherm*d.uenerg/3./k/N

				radplot, densav = d.azAverage(np.sqrt(r2[members&gas]), d.disc.density[members&gas]*d.udens)
				_, tempav = d.azAverage(np.sqrt(r2[members&gas]), temp[members&gas])
				raddust, fdustav = d.azAverage(np.sqrt(r2[members&gas]), dustfrac[members&gas],nbins=20)

				#_, denspoly = d.fitPolytrope(radplot,densav)
				#_, temppoly = d.fitPolytrope(radplot,tempav)
				#raddust, fdustpoly = d.fitPolytrope(raddust,fdustav)

				timelegend.append('%.2f yrs' %d.disc.time)

				#Plot density vs. radius
				plt.figure(1)
				plt.scatter(np.sqrt(r2[members & gas]), d.disc.density[members & gas]*d.udens, s=0.1)
				plt.plot(polyx[i], max(d.disc.density[members&gas]*d.udens)*polytheta[i]**1.5) #plot polytropic profile
				plt.ylim([1e-12, max(d.disc.density[members&gas]*d.udens*2)])
				plt.xlim([0,3.])
				plt.ylabel(r'$\rho$ (g cm$^{-3}$)')
				plt.xlabel('Dist from clump centre (au)')
				plt.title('Gas density vs. radius')
				plt.yscale('log')
				plt.legend(timelegend)

				#Plot temperature vs. radius
				plt.figure(2)
				plt.scatter(np.sqrt(r2[members & gas]), temp[members & gas], s=0.1)
				plt.plot(polyx[i], max(temp[members&gas])*polytheta[i]) #plot polytropic profile
				plt.ylim([10, max(temp[members&gas]*1.1)])
				plt.xlim([0,3.])
				plt.ylabel('Temp (K)')
				plt.xlabel('Dist from clump centre (au)')
				plt.yscale('log')
				plt.title('Gas temp vs. radius')
				plt.legend(timelegend)

				'''
				#Plot dust fraction vs. radius
				plt.figure(3)
				plt.scatter(np.sqrt(r2[members & gas]), dustfrac[members & gas], s=0.1)
				plt.plot(raddust, fdustav)
				#plt.plot(raddust, fdustpoly)
				plt.ylim([0, max(dustfrac[members&gas]*1.1)])
				plt.xlim([0,3.])
				plt.ylabel('Dust-to-gas ratio')
				plt.xlabel('Dist from clump centre (au)')
				#plt.yscale('log')
				plt.title('Dust fraction vs. radius')
				plt.legend(timelegend)
				'''
			plt.show()
			plt.clf()

			keepplotting = bool(input('Change polytropic indices and replot? (y=1,n=0)?'))
			if keepplotting:
				polyi = input('Select polytropic indices for dumps (1 - %i): ' %len(IDs))
				try:
					polyi = [float(polyi)]
				except:
					polyi = map(float, polyi)
				polyx, polytheta = Analyser.polytropicProfile(n=polyi)

