from track import Tracker
from analyse import Analyser

def readParams():
	f = open('tracker.params','r')
	params = f.read().splitlines()
	datadir = params[0]						#working directory with all the disc files
	final = int(params[1])					#final dump ID
	ctrack = map(int, params[2].split())	#list of the clump IDs in final disc I wish to track
	dt = float(params[3])					#timestep between disc dumps
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

tracks = Tracker.runTracker(params)

wd = params[0]

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
'''
wd = '/disk2/cadman/phantom_runs/dustyfragment/500cm_slow_2'
disc = Analyser(200,wd)
Analyser.plotClump(disc,2)
import pdb; pdb.set_trace()
'''