import numpy as np

"""
Takes the clump information from the final disc dump and iterates backward until this
    clump can no longer be located (i.e. when the clump has first formed.)

We will compare 3 things for each clump
i- 50% of the members of the previous clumps are present in the new clump
ii- The new clump has a COM in the location we expect it to be
(iii- The most bound particle from the previous clump is still present)
If the first two conditions are satisfied, we have located the clump in the previous
    dump and can keep moving backward.
If a clump cannot be located in the previous (10) dumps to where it was last found, 
    then this clump no longer exists and we have found its origin.
"""

class Tracker(object):

	def __init__(self,idisc,wd):
		if idisc < 10:
			disc = 'disc00%i' %idisc
		elif idisc < 100:
			disc = 'disc0%i' %idisc
		else:
			disc = 'disc%i' %idisc
		self.disc = disc
		self.wd = wd
		self.readMembers()
		self.readProperties()

	def readMembers(self):
		"""
		Read in binary clumpmembers file from Clumpfind output.
		"""
		f = open('%s/raw_clumpmembers_%s' %(self.wd,self.disc), 'rb')
		#Skip first and last entries from this array.
		data = np.fromfile(f, dtype='i')[1:-1]
		self.nclumps = max(data)
		members = {}
		# I think we don't want ID==0 as this refers to no clump (CHECK...)
		for clump in range(self.nclumps):
			membershold = np.argwhere(data==clump+1).flatten()
			members[clump] = membershold
		self.members = members

	def readProperties(self):
		"""
		Read in clumpcat file with locations and info on clumps
		First 3 values in clumpcat rows are x,y,z positions
		Next 3 values are vx,vy,vz
		"""
		f = open('%s/raw_clumpcat_%s' %(self.wd,self.disc),'r')
		self.locs = []
		self.vels = []
		self.mbp = []	
		for line in f.readlines()[1:]:
			vals = np.array(line.split(' '))
			vals = vals[vals!='']
			self.locs.append([float(vals[0]), float(vals[1]), float(vals[2])])
			self.vels.append([float(vals[3]), float(vals[4]), float(vals[5])])
			self.mbp.append(float(vals[9]))
		self.locs = np.array(self.locs)
		self.vels = np.array(self.vels)
		self.mbp = np.array(self.mbp)

	@staticmethod
	def compareMembers(clump1,clump2,step=1.):
		"""
		Check if at least 50% of the members in clump1 are also in clump2.
		If they are then they are the same clump, return True.
		Otherwise return False.
		"""
		same = False
		if step<=5:
			if (len(np.intersect1d(clump1,clump2)) > (0.51-step/100.)*len(clump1) or
				len(np.intersect1d(clump1,clump2)) > (0.51-step/100.)*len(clump2)):
				same = True
		else:
			#if we are getting desperate, only require 25% of the members are the same
			if (len(np.intersect1d(clump1,clump2)) > (0.25)*len(clump1) or
				len(np.intersect1d(clump1,clump2)) > (0.25)*len(clump2)):
				same = True
		return same

	@staticmethod
	def comparePositions(loc1,vel1,loc2,dt,step=1.):
		"""
		Check if clump2 is roughly in the location we would expect it to be
		"""
		same = False
		if step<=5:
			#check if x and y coordinates are in expected positions (don't require z pos)
			if np.allclose(loc1[:2] - vel1[:2]*dt*step, loc2[:2], rtol=0.2+step/100.):
				same = True
		else:
			#if we are getting desparate, only require x OR y co-ordincates match
			if sum(np.isclose(loc1[:2] - vel1[:2]*dt*step, loc2[:2], rtol=0.2+step/100.)) >= 1:
				same = True
		return same

	@staticmethod
	def compareMBP(clump1,mbp1,clump2,mbp2):
		#Might not use this as all clumps have mbp in 
		same = False
		if mbp1 in clump2 or mbp2 in clump2:
			same = True
		return same

	@staticmethod
	def findClumpsinDump(disc_curr,disc_prev,ID,dt,step=1.):
		"""
		routine to find clump (ID) from disc_curr in disc_prev
		Iterate through all the clumps in disc_prev, and compare members, positions and MBP to clump ID
			in disc_curr.members[ID].
		If the clump is found, append to prevID and make disc_curr=disc_prev and move on.
		If the clump is not found in dump i-1, move onto i-2 and make step+=1. Do the same method to try and find the
			clump again.
		If step > 5, make requirements less strict and start comparing MBP.
		If two potential clumps are identified in the previous dump then check which one is the most likely candidate.
		"""
		same = 0
		found = False
		prevID = []
		# iterate through clumps in previous disc, identifying cadidates for clump progenitor
		for clump2 in disc_prev.members:
			same = 0
			#compare if the two clumps have the same members
			if Tracker.compareMembers(disc_curr.members[ID],disc_prev.members[clump2],step):
				same+=1
			#compare if they have the same positions
			if Tracker.comparePositions(disc_curr.locs[ID],disc_curr.vels[ID],disc_prev.locs[clump2],dt,step):
				same+=1
			if step >= 5:
				#if we are getting desparate, allow same mbp to +1 same
				if Tracker.compareMBP(disc_curr.members[ID],disc_curr.mbp[ID],
							  disc_prev.members[clump2],disc_prev.mbp[clump2]):
					same+=1
			if same==2:
				#candidate has been found in previous dump, save it and move onto next clump
				found=True
				prevID.append(clump2)
		if not found:
			#clump has not been found
			prevID.append(np.nan)

		# make sure the central point mass hasn't been identified as a candidate
		if found and 0 in prevID:
			prevID.remove(0)
			# if that was the only cadidate then set found = False again
			if len(prevID) == 0:
				found = False
				prevID = [np.nan]

		#check if multiple candidates have been found, if they have then do further checks.
		if len(prevID) > 1:
			chosen = False
			#do either candidate have the same mbp as the original clump?
			for clump2 in prevID:
				if Tracker.compareMBP(disc_curr.members[ID],disc_curr.mbp[ID],
									  disc_prev.members[clump2],disc_prev.mbp[clump2]):
					prevID = clump2
					chosen = True
			#if we still can't decide, i.e. neither have mbp in original clump, just choose the larger candidate.
			#i.e. which of the candidates have more members?
			if chosen==False:
				bigger = max([disc_prev.members[clump2].size for clump2 in prevID])
				prevID = prevID[prevID==bigger]
		else:
			# only one candidate has been identified
			prevID = prevID[0]

		return prevID, found

	@staticmethod
	def runTracker(params):
		"""
		Iterate through dump files, running findClumpsinDump method to identify desired clump in previous dump
		"""
		wd, initial, final, ctrack, dt = params
		#curr_disc = fdisc
		track = {}
		for i, ID in enumerate(ctrack):
			idump = final
			ogID = ID
			track[ogID] = [idump,ID]
			#go backward to dump 0, trying to find these members and positions in each dump
			while idump >= initial:
				curr_disc = Tracker(idump,wd)
				#if idump-1==198 and ogID==5:
				#	#this clump is causing problems, skip it
				#	idump=198
				prev_disc = Tracker(idump-1,wd)
				found = False
				prevID, found = Tracker.findClumpsinDump(curr_disc,prev_disc,ID,dt)
				if found:
					print ('Clump %i found in disc%i with ID %i at pos %.4f %.4f' %(ctrack[i], idump-1, prevID,
																		 			prev_disc.locs[prevID][0], 
																		 			prev_disc.locs[prevID][1]))
					ID = prevID
					idump -= 1
					track[ogID] += [idump, ID]
					continue
				else:
					print ('Clump %i not found in disc%i' %(ctrack[i], idump-1))
					#look for it in the previous 10 discs
					for step in [2.,3.,4.,5.,6.,7.,8.,9.,10.]:
						prev_disc = Tracker(idump-step,wd)
						prevID, found = Tracker.findClumpsinDump(curr_disc,prev_disc,ID,dt,step)
						if found:
							print ('Clump %i found in disc%i with ID %i at pos %.4f %.4f' %(ctrack[i], idump-step, prevID,
																		 					prev_disc.locs[prevID][0], 
																		 					prev_disc.locs[prevID][1]))
							ID = prevID
							idump -= step
							track[ogID] += [idump, ID]
							break
						else:
							print ('Clump %i not found in disc%i' %(ctrack[i], idump-step))
					if found:
						continue
					if not found:
						print ('Clump %i can no longer be found, moving on...' %(ctrack[i]))
						break
			track[ogID] = np.reshape(track[ogID], [len(track[ogID])/2,2])
		return track
