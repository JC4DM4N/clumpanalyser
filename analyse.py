"""
Analyse and plot disc fragment properties.
"""
import numpy as np
from matplotlib import cm
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
		#self.readMembers()
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

        def readMembers(self):
                """
                Read in binary clumpmembers file from Clumpfind output.
                """
                f = open('%s/raw_clumpmembers_%s' %(self.wd,self.file), 'rb')
                #Skip first and last entries from this array.
                data = np.fromfile(f, dtype='i')[1:-1]
                self.nclumps = max(data)
                members = {}
                # I think we don't want ID==0 as this refers to no clump (CHECK...)
                for clump in range(self.nclumps):
                        membershold = np.argwhere(data==clump+1).flatten()
                        members[clump] = membershold
                self.members = members

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

	@staticmethod
	def writeClumpstoDumps(discs,IDs):
		'''
		Write text file of particle properties within each clump
		'''
                timelegend = []
                for i, [d, ID] in enumerate(zip(discs,IDs)):
                        clumpxyz = d.clumpcat[ID][:3]
                        r2 = (d.disc.xyzh[0]-clumpxyz[0])**2 + (d.disc.xyzh[1]-clumpxyz[1])**2
                        members = np.sqrt(r2) < d.annulus                       #members are all particles within radial annulus

                        gas = d.disc.itype == 1
                        dust = d.disc.itype == 2

                        dustfrac = d.disc.dustfrac*1e8  #as I originally set dust-to-gas=1e-10

                        #Calculate temperatures from thermal energies
                        k = 1.38064852e-16    #ergs
                        mH = 1.6735575e-24    #grams
                        gmw = 2.381           #mean mass taken from Phantom
                        N = sum(gas*d.disc.massofgas)*d.umass/mH/gmw #number of atoms
                        temp = 2.*d.disc.utherm*d.uenerg/3./k/N

			utime = d.utime/(60*60*24*365.25)

			#create arrays of particle masses
			mass = np.zeros(len(d.disc.xyzh[0,:]))
			mass[d.disc.itype == 1] = d.disc.massofgas
			mass[d.disc.itype == 2] = d.disc.massofdust

			clumpdata = zip(d.disc.xyzh[0,members], d.disc.xyzh[1,members], d.disc.xyzh[2,members], d.disc.xyzh[3,members], 
					d.disc.density[members], mass[members], temp[members], dustfrac[members], d.disc.itype[members])
			clumpdata = np.asarray(clumpdata)
			header = ("time: %s utime (yrs^-1): %s \n x, y, z, h, density, mass, temp, " %(str(d.disc.time), str(utime)) +
				  "dustfrac, itype \n %s, %s, %s, %s, %s, %s, 0.0, 0.0, 0.0 \n \n" %(str(d.udist), str(d.udist),
				 								     str(d.udist), str(d.udist),
												     str(d.udens), str(d.umass)))
			np.savetxt('clumpfiles/clumpdata_%.0f.txt' %d.disc.time, clumpdata, header=header)

	@staticmethod
	def plotDustFrac(discs,IDs):
        	timelegend = []
                for i, [d, ID] in enumerate(zip(discs,IDs)):
                	clumpxyz = d.clumpcat[ID][:3]
                        r2 = (d.disc.xyzh[0]-clumpxyz[0])**2 + (d.disc.xyzh[1]-clumpxyz[1])**2
                        members = np.sqrt(r2) < d.annulus                       #members are all particles within radial annulus
                        inner = np.sqrt(r2) < 0.25

                        gas = d.disc.itype == 1
                        dust = d.disc.itype == 2
                        dustfrac = d.disc.dustfrac*1e8  #as I originally set dust-to-gas=1e-10

                        #Calculate temperatures from thermal energies
                        k = 1.38064852e-16    #ergs
                        mH = 1.6735575e-24    #grams
                        gmw = 2.381           #mean mass taken from Phantom
                        N = sum(gas*d.disc.massofgas)*d.umass/mH/gmw #number of atoms
                        temp = 2.*d.disc.utherm*d.uenerg/3./k/N

                        timelegend.append('%.2f yrs' %d.disc.time)
			#Plot dust fraction vs. radius
                        plt.figure(3)
                        plt.scatter(np.sqrt(r2[members & gas]), dustfrac[members & gas], s=1.0,
                                    c=temp[members & gas], cmap=cm.magma)
                        plt.yscale('log')
                        plt.ylim([0.001,1.0])
                        plt.xlim([0,3.])
                        plt.ylabel('Dust-to-gas ratio')
                        plt.xlabel('Dist from clump centre (au)')
                        plt.title('Dust fraction vs. radius')
                        cbar = plt.colorbar()
                        cbar.set_label('Temp (K)')
			plt.clim([0,400])
                        plt.legend(['%.2f yrs' %d.disc.time])
			plt.savefig('/disk2/cadman/phantom_runs/dustyfragment/fragplots/5cmhires_dustfrac/%i.png' %(len(IDs)-i))
			#plt.show()
			plt.clf()

	@staticmethod
	def plotMultiClumps(discs,IDs):
		'''
		Plot radial profiles of clump properties (dens, temp, dustfrac).
		'''

		#Initially plot polytropic indices n=1.5 for all clumps
		polyx, polytheta = Analyser.polytropicProfile(n=np.ones(len(IDs))*1.5)

		#plot colors
		colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

		keepplotting = True
		while keepplotting:
			timelegend = []
			for i, [d, ID] in enumerate(zip(discs,IDs)):
				clumpxyz = d.clumpcat[ID][:3]
				r2 = (d.disc.xyzh[0]-clumpxyz[0])**2 + (d.disc.xyzh[1]-clumpxyz[1])**2
				members = np.sqrt(r2) < d.annulus 			#members are all particles within radial annulus
				inner = np.sqrt(r2) < 0.25

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
				#plt.scatter(np.sqrt(r2[members & gas]), d.disc.density[members & gas]*d.udens, s=0.1)
				plt.plot(radplot, densav,'-',c=colors[i], label=timelegend[i])
				plt.plot(polyx[i], max(d.disc.density[members&gas&inner]*d.udens)*polytheta[i]**1.5,'--',c=colors[i]) #polytropic profile
				plt.ylim([1e-13, max(d.disc.density[members&gas]*d.udens*2)])
				plt.xlim([0,3.])
				plt.ylabel(r'$\rho$ (g cm$^{-3}$)')
				plt.xlabel('Dist from clump centre (au)')
				plt.title('Gas density vs. radius')
				plt.yscale('log')
				plt.legend()

				#Plot temperature vs. radius
				plt.figure(2)
				#plt.scatter(np.sqrt(r2[members & gas]), temp[members & gas], s=0.1)
				plt.plot(radplot, tempav, '-',c=colors[i], label=timelegend[i])
				plt.plot(polyx[i], max(temp[members&gas&inner])*polytheta[i], '--', c=colors[i]) #polytropic profile
				plt.ylim([10, max(temp[members&gas]*1.1)])
				plt.xlim([0,3.])
				plt.ylabel('Temp (K)')
				plt.xlabel('Dist from clump centre (au)')
				plt.yscale('log')
				plt.title('Gas temp vs. radius')
				plt.legend()
				
				#Plot dust fraction vs. radius
				plt.figure(3)
				plt.scatter(np.sqrt(r2[members & gas]), dustfrac[members & gas], s=1.0, 
					    c=temp[members & gas], cmap=cm.magma)
				plt.plot(raddust, fdustav, '-', c=colors[i], label=timelegend[i])
				plt.ylim([0, max(dustfrac[members&gas]*1.1)])
				plt.xlim([0,3.])
				plt.ylabel('Dust-to-gas ratio')
				plt.xlabel('Dist from clump centre (au)')
				#plt.yscale('log')
				plt.title('Dust fraction vs. radius')
				cbar = plt.colorbar()
				cbar.set_label('Temp (K)')
				plt.legend()
				
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

