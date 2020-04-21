import os
from track import Tracker
from analyse import Analyser

def readParams():
	f = open('tracker.params','r')
	params = f.read().splitlines()
	datadir = params[0].split('#')[0].split()[0]				#working directory with all the disc files
	initial = int(params[1].split('#')[0])					#final dump ID
	final = int(params[2].split('#')[0])					#final dump ID
	ctrack = map(int, params[3].split('#',1)[0].split())			#list of the clump IDs in final disc I wish to track
	dt = float(params[4].split('#',1)[0])					#timestep between disc dumps
	return [datadir, initial, final, ctrack, dt]

params = readParams()

#ask user which of input clumps we actually wish to track on this run
ctrack = input('Select which of clumps %s to track: ' %params[3])

try:
	#if only one clump chosen
	ctrack = [int(ctrack)]
except:
	#if more than one clump chosen
	ctrack = map(int, ctrack)

params[3] = ctrack

#Run clump tracker rountine

tracks = Tracker.runTracker(params)

wd = params[0]

#Iterate through dumps, create analyser instances, and plot as desired...

discs=[]
iclumps=tracks[ctrack[0]][:,1]
iclumps=map(int, iclumps)

#Check necessary directories for saving plots exist
if not os.path.exists('%s/clumpfiles' %wd):
    os.mkdir('%s/clumpfiles' %wd)
if not os.path.exists('%s/fragplots' %wd):
    os.mkdir('%s/fragplots' %wd)

#Ask what user wishes to plot/analyse
writeClump = False
plotDust = True

if writeClump or plotDust:
    for i, idump in enumerate(tracks[ctrack[0]]):
        disc = Analyser(idump[0],wd)

        if writeClump:
            #Write clump particle data to text files...
	    print('Writing clump in %i to .txt file...' %idump[0])
            disc.writeClumptoDump(iclumps[i])
        if plotDust:
            #Plot dust fraction coloured by temperature...
	    print('Plotting clump dust fraction from disc %i...' %idump[0])
            disc.plotDustFrac(iclumps[i],len(iclumps)-i)


#Plot clump data for desired dumps files (temperature, density, dustfrac)

keepplotting=True
while keepplotting:
	for frag in tracks:
		discs = []
		iclumps = []
		#ask user which disc dumps they now wish to analyse for each clump
		dumps = input('Select which disc files to plot polytropic profiles for clump %i: ' %frag)
		try:
			dumps = [int(dumps)]
		except:
			dumps = map(int, dumps)
		for idump in dumps:
			#create analyser instanes for each desired dump
			discs.append(Analyser(idump,wd))
			iclumps.append(int(tracks[frag][tracks[frag][:,0]==idump,1]))

		print 'Plotting polytropic profile for Clump %i...' %int(frag)
		Analyser.plotMultiClumps(discs,iclumps)

	keepplotting = bool(input('Would you like to plot polytropic profiles of other dumps? (y=1,n=0)?'))

