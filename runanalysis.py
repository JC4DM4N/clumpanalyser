from track import Tracker
from analyse import Analyser

def readParams():
	f = open('tracker.params','r')
	params = f.read().splitlines()
	import pdb; pdb.set_trace()
	datadir = params[0].split('#')[0].split()[0]				#working directory with all the disc files
	final = int(params[1].split('#')[0])					#final dump ID
	ctrack = map(int, params[2].split('#',1)[0].split())			#list of the clump IDs in final disc I wish to track
	dt = float(params[3].split('#',1)[0])					#timestep between disc dumps
	return [datadir, final, ctrack, dt]

params = readParams()

#ask user which of input clumps we actually wish to track on this run
ctrack = input('Select which of clumps %s to track: ' %params[2])

try:
	#if only one clump chosen
	ctrack = [int(ctrack)]
except:
	#if more than one clump chosen
	ctrack = map(int, ctrack)

params[2] = ctrack

###Run clump tracker rountine

tracks = Tracker.runTracker(params)

wd = params[0]


###Create Analyser instances for each disc file

discs=[]
iclumps=tracks[ctrack[0]][:,1]
iclumps=map(int, iclumps)
for i, idump in enumerate(tracks[ctrack[0]]):
    discs.append(Analyser(idump[0],wd))


###Write clump particle data to text files...

Analyser.writeClumpstoDumps(discs,iclumps)


###Plot dust fraction coloured by temperature...

Analyser.plotDustFrac(discs,iclumps)


###Plot clump data for desired dumps files (temperature, density, dustfrac)

keepplotting=True
while keepplotting:
	for frag in tracks:
		discs = []
		iclumps = []
		#ask user which disc dumps they now wish to analyse for each clump
		dumps = input('Select which disc files to plot for clump %i: ' %frag)
		try:
			dumps = [int(dumps)]
		except:
			dumps = map(int, dumps)
		for idump in dumps:
			discs.append(Analyser(idump,wd))
			iclumps.append(int(tracks[frag][tracks[frag][:,0]==idump,1]))

		print 'Plotting for Clump %i...' %int(frag)
		Analyser.plotMultiClumps(discs,iclumps)

	keepplotting = bool(input('Would you like to plot other dumps? (y=1,n=0)?'))

