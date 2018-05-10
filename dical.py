import math
import random
import os
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt


random.seed (4711)


def writeIsolationMigrationDemographyFile (filename, migrationAfterSplit):
	# this is a model where ancestral population of certain size splits into two (individual sizes), with gene flow possible after split (or not)
	ofs = open (filename, 'w')
	ofs.write (" # boundary points of the epochs [0,t_1,t_2,infinity) [intervals of constant demography]\n")
	ofs.write ("[ ?0 ]\n")
	ofs.write ("# EPOCH 1\n")
	ofs.write ("# population structure\n")
	ofs.write ("{{0},{1}}\n")
	ofs.write ("# population sizes\n")
	ofs.write ("?1	?2\n")
	ofs.write ("# instaneous migration rates at beginnig of epoch\n")
	ofs.write ("null\n")
	ofs.write ("# migration rates during epoch\n")
	if (migrationAfterSplit):
		ofs.write ("0	?3\n")
		ofs.write ("?3	0\n")
	else:
		ofs.write ("0	?3\n")
		ofs.write ("?3	0\n")
	ofs.write ("# EPOCH 2\n")
	ofs.write ("# population structure\n")
	ofs.write ("{{0,1}}\n")
	ofs.write ("# population sizes\n")
	ofs.write ("?4\n")
	ofs.write ("# instaneous migration rates at beginnig of epoch\n")
	ofs.write ("null\n")
	ofs.write ("# migration rates during epoch\n")
	ofs.write ("0\n")
	# be nice
	ofs.close ()


def writePieceWiseConstantDemographyFile (filename, numPieces, endChangeTimes, globalParamOffset):

	for p in range(len(endChangeTimes)-1):
		if (endChangeTimes[p] > endChangeTimes[p+1]):
			raise Exception ("In constant demography: Change times not ordered correctly.")

	# this is a model where ancestral population of certain size splits into two (individual sizes), with gene flow possible after split (or not)
	ofs = open (filename, 'w')

	ofs.write (" # boundary points of the epochs [0,t_1,t_2,infinity) [intervals of constant demography]\n")
	numFixedTimes = len(endChangeTimes)
	# only as many as not given
	estimateFirstTimes = ["?%d" % (x+globalParamOffset) for x in range(numPieces-numFixedTimes-1)]
	# the given ones
	fixedEndTimes = [ "%.6e" % x for x in endChangeTimes]
	# and put together
	timeString = ", ".join (estimateFirstTimes + fixedEndTimes)
	ofs.write ("[ " + timeString + " ]\n")

	sizeParamOffset = numPieces-1-numFixedTimes
	for p in range(numPieces):
		ofs.write ("# EPOCH %d\n" % (p+1))
		ofs.write ("# population structure\n")
		ofs.write ("{{0}}\n")
		ofs.write ("# population sizes\n")
		ofs.write ("?%d\n" % (globalParamOffset + sizeParamOffset + p))
		ofs.write ("# instaneous migration rates at beginnig of epoch\n")
		ofs.write ("null\n")
		ofs.write ("# migration rates during epoch\n")
		ofs.write ("0\n")

	# be nice
	ofs.close ()


def writeOnePopExpRateRecentFile (filename, numPieces):
	if (numPieces < 1):
		raise Exception ("In exp rates: Need at least one epoch.")
	ofs = open (filename, 'w')
	ofs.write ("# exponential growth rate EPOCH 1\n")
	ofs.write ("?0\n")
	for p in range(1,numPieces):
		ofs.write ("# exponential growth rate EPOCH %d\n" % (p+1))
		ofs.write ("0\n")

	# be nice
	ofs.close ()


def writeMutRecoParameterFile (filename, mutRate, recoRate):
	ofs = open (filename, 'w')
	ofs.write ("# mutation rate\n")
	ofs.write ("%s\n" % str(mutRate))
	ofs.write ("# recombination rate\n")
	ofs.write ("%s\n" % str(recoRate))
	ofs.write ("# mutation matrix (2 alleles)\n")
	ofs.write ("0	1\n")
	ofs.write ("1	0\n")
	# be nice
	ofs.close ()


def writeConfigFile (filename, numLoci, numAlleles, numDemes, numHapsInDeme):
	if (numDemes != len(numHapsInDeme)):
		raise Exception ("In config file: Num demes does not match given haplotype numbers.")
	# open file
	ofs = open (filename, 'w')
	# config header (vcf ignores numLoci anyways)
	ofs.write ("# numLoci\tnumAlleles\tnumDemes\t (numLoci is ignored for vcf)\n")
	ofs.write ("%d\t%d\t%d\n" % (numLoci, numAlleles, numDemes))
	ofs.write ("# one line for each haplotype to how many in which deme (each diploid is considered as two seperate haplotypes)\n")
	for d in range(len(numHapsInDeme)):
		for h in range(numHapsInDeme[d]):
			nVec = [0]*numDemes
			nVec[d] = 1
			ofs.write ("\t".join([str(x) for x in nVec]) + "\n")
	# be nice
	ofs.close ()


def makeReferenceFile (filename, numLoci):
	# ancestral is all A
	refSeq = ["A"]*numLoci
	ofs = open (filename, 'w')
	ofs.write ("".join (refSeq) + "\n")
	ofs.close ()


def theta (Nref, perGenMutProb):
	return Nref * 4 * perGenMutProb

def rho (Nref, perGenRecoProb):
	return Nref * 4 * perGenRecoProb

def DiploidSizeToCoalSize (Nref, size):
	return size/float(Nref)

def CoalSizeToDiploidSize (Nref, size):
	return size * Nref

def GenTimeToCoalTime (Nref, gen):
	return gen / float(2 * Nref)

def CoalTimeToGenTime (Nref, coalTime):
	return coalTime * 2 * Nref

def PerGenPercentToExpRate (Nref, perGenRate):
	return perGenRate / float(100) * 2 * Nref

def ExpRateToPerGenPercent (Nref, expRate):
	return expRate / float(2*Nref) * 100

def PerGenMigToMigRate (Nref, perGenMigProb):
	return perGenMigProb * 4 * Nref

def MigRateToPerGenMig (Nref, migRate):
	return migRate / float(4*Nref)


def logUniformStartingPointsFile (filename, numPoints, bounds):
	# open file
	ofs = open (filename, 'w')
	# get the logBounds
	logBounds = [[math.log(b[0]), math.log(b[1])] for b in bounds]
	# and write the points
	for i in range(numPoints):
		# random points in log scale
		logPoint = [random.uniform(b[0], b[1]) for b in logBounds]
		point = [math.exp(x) for x in logPoint]
		ofs.write ("\t".join (["%.10f" % x for x in point]) + "\n")
	# be nice
	ofs.close ()


def logGrid (firstTime, lastTime, numPoints):
	return [math.exp(math.log(firstTime) + x/float(numPoints-1) * (math.log(lastTime) - math.log(firstTime))) for x in range(numPoints)]


def batchify (realLines):

	# get the preTraces
	preTraces = {}
	for line in realLines:
		# print ("======================")
		# print (preTraces)
		metaString = line.split()[-1]
		[gen, emStep, ind] = [int(x) for x in metaString.split("[")[1].split("]")[0].split("_")]
		if (gen not in preTraces.keys()):
			preTraces[gen] = {}
		if (ind not in preTraces[gen].keys()):
			preTraces[gen][ind] = {}
		# finally put it in
		if (emStep in preTraces[gen][ind]):
			raise Exception ("In batchify: Output inconsistent")
		preTraces[gen][ind][emStep] = line

	# convert it into real traces
	traces = {}
	for (gen,v1) in preTraces.items():

		if (gen not in traces.keys()):
			traces[gen] = {}

		for (ind,v2) in v1.items():

			if (ind in traces[gen].keys()):
				raise Exception ("In batchify: Output inconsistent")

			# make sure that steps for this individual are contiguous
			allSteps = v2.keys()
			minStep = min(allSteps)
			maxStep = max(allSteps)
			if ((minStep != 0) or (len(set(allSteps)) != (maxStep+1)) or (len(allSteps) != (maxStep+1))):
				raise Exception ("In batchify: Output inconsistent")

			# fill it up
			traces[gen][ind] = [0]*(maxStep+1)
			for (step, values) in v2.items():

				# make a pair with (likeliood, point)
				daValues = values.split("\t")
				# likelihood is first
				likelihood = float(daValues[0])
				# everything from 3rd until second to last should belong to point
				point = [float(x) for x in daValues[2:-1]]

				traces[gen][ind][step] = (likelihood, point)

			# see whether all points of equal length
			thisPoints = traces[gen][ind]
			for i in range(1,len(thisPoints)):
				# print(thisPoints[i])
				if (len(thisPoints[i-1]) != len(thisPoints[i])):
					raise Exception ("In batchify: Output inconsistent")

	return traces


def returnMLE (outputFilename):
	# now input
	ifs = open (outputFilename)
	# get only lines without "#"
	realLines = []
	for line in ifs:
		if ((line == "\n") or (line.startswith ("#"))):
			continue
		realLines.append (line)


	# batchify
	traces = batchify (realLines)


	# show me what you got
	mle = []
	maxLike = float("-inf")
	# maxLike = float("+inf")
	for (gen,v1) in traces.items():
		for (ind, v2) in v1.items():
			# print ("++++++++++ (%d, %d): " % (gen,ind))
			for i in range(len(v2)):
				# look at likelihoods
				# print (v2[i])
				thisOne = v2[i]
				if (thisOne[0] > maxLike):
				# if (thisOne[0] < maxLike):
					maxLike = thisOne[0]
					mle = thisOne[1]

	# should have a good mle
	return (mle, maxLike)


# class ExpGrowthEstimateOnsetTimeAnalysis:

# 	metaNumStartPoints = 10
# 	metaNumPoints = 5
# 	metaKeepBest = 2
# 	metaNumIterations = 3
# 	numEMiterations = 3
# 	numIterationsMstep = 3

# 	cl = "lol"

# 	numLociPerHmmStep = 500

# 	Nref = 10000
# 	mu = 1.25e-8
# 	r = 1.25e-8

# 	numEpochs = 3
# 	# firstFixedTimeInGen = 1000
# 	# lastFixedTimeInGen = 3000
# 	fixedTimeInGen = 2000

# 	# bounds related stuff
# 	minPerGenGrowth = 0.0001
# 	maxPerGenGrowth = 2
# 	minBoundTimeInGen = 100
# 	maxBoundTimeInGen = 1000
# 	maxSize = 100000
# 	minSize = 100

# 	jarFile = "cleanDiCal2.jar"


# 	def __init__ (self, uniqueBasename, numCores, vcfFiles, refFiles, sampleSize, randomSeed):
# 		self.vcfFileList = vcfFiles
# 		self.refFileList = refFiles
# 		self.seeed = randomSeed
# 		self.numCores = numCores

# 		self.diCalOutputFileName = "%s.dical_out" % uniqueBasename
# 		self.paramFile = "%s.param" % uniqueBasename
# 		self.demoFile = "%s.demo" % uniqueBasename
# 		self.rateFile = "%s.rates" % uniqueBasename
# 		self.configFile = "%s.config" % uniqueBasename

# 		# put together the bounds
# 		# one parameter per epoch
# 		self.bounds = []
# 		self.bounds.extend ([[PerGenPercentToExpRate (self.Nref, self.minPerGenGrowth), PerGenPercentToExpRate (self.Nref, self.maxPerGenGrowth)]])
# 		self.bounds.extend ([[GenTimeToCoalTime (self.Nref, self.minBoundTimeInGen), GenTimeToCoalTime (self.Nref, self.maxBoundTimeInGen)]])
# 		# self.bounds.extend ([[GenTimeToCoalTime (self.Nref, 10*self.minBoundTimeInGen), GenTimeToCoalTime (self.Nref, 10*self.maxBoundTimeInGen)]])
# 		self.bounds.extend ([[DiploidSizeToCoalSize (self.Nref, self.minSize), DiploidSizeToCoalSize (self.Nref, self.maxSize)]]*self.numEpochs)

# 		# random starting points
# 		self.metaStartFile = "%s.rand" % uniqueBasename
# 		samplingRegionFactor = 0.25
# 		samplingRegion = [[math.exp (math.log(bs[0]) + samplingRegionFactor * (math.log(bs[1]) - math.log(bs[0]))), math.exp (math.log(bs[0]) + (1-samplingRegionFactor) * (math.log(bs[1]) - math.log(bs[0])))] for bs in self.bounds]
# 		logUniformStartingPointsFile (self.metaStartFile, self.metaNumStartPoints, samplingRegion)

# 		# param file
# 		writeMutRecoParameterFile (self.paramFile, theta (self.Nref, self.mu), rho (self.Nref, self.r))

# 		# config file
# 		# numLoci doesn't even matter, also, two alleles and one deme
# 		writeConfigFile (self.configFile, 123456, 2, 1, [sampleSize])

# 		# demography file
# 		# firstTime = GenTimeToCoalTime (self.Nref, self.firsTimeInGen)
# 		# lastTime = GenTimeToCoalTime (self.Nref, self.lastTimeInGen)
# 		# themTimes = [math.exp(math.log(firstTime) + x/float(self.numEpochs-2) * (math.log(lastTime) - math.log(firstTime))) for x in range(self.numEpochs-1)]
# 		# at the moment, there is only one
# 		givenTimes = [GenTimeToCoalTime (self.Nref, self.fixedTimeInGen)]
# 		writePieceWiseConstantDemographyFile (self.demoFile, self.numEpochs, givenTimes, 1)
# 		writeOnePopExpRateRecentFile (self.rateFile, self.numEpochs)


# 	def run (self):

# 		if (self.numCores > 1):
# 			parallelEmSteps = math.ceil (self.numCores / 2)
# 			parallelString = " --parallel %d --metaParallelEmSteps %d" % (self.numCores, parallelEmSteps)
# 		else:
# 			parallelString = ""

# 		dicalCmd = ("java -Xmx30G -jar %s" +
# 			" --paramFile %s" +
# 			" --demoFile %s" +
# 			" --ratesFile %s" +
# 			" --configFile %s" +
# 			" --vcfFile %s" +
# 			" --vcfFilterPassString 'PASS'" +
# 			" --vcfReferenceFile %s" +
# 			" --seed %d" +
# 			" --lociPerHmmStep %d" +
# 			" --compositeLikelihood %s" +
# 			" --metaStartFile %s" +
# 			" --metaNumIterations %d" +
# 			" --metaKeepBest %d" +
# 			" --metaNumPoints %d" +
# 			" --numberIterationsEM %d" +
# 			" --numberIterationsMstep %d" +
# 			" --printEmPath" +
# 			# " --intervalType simple" +
# 			" --intervalType logUniform" +
# 			" --intervalParams '11,0.01,4'" +
# 			" --bounds '%s'" +
# 			"%s" +
# 			" >%s") % (
# 			self.jarFile,
# 			self.paramFile,
# 			self.demoFile,
# 			self.rateFile,
# 			self.configFile, 
# 			self.vcfFileList,
# 			self.refFileList,
# 			self.seeed,
# 			self.numLociPerHmmStep,
# 			self.cl,
# 			self.metaStartFile,
# 			self.metaNumIterations,
# 			self.metaKeepBest,
# 			self.metaNumPoints,
# 			self.numEMiterations,
# 			self.numIterationsMstep,
# 			";".join([",".join([str(x) for x in y]) for y in self.bounds]),
# 			parallelString,
# 			self.diCalOutputFileName)

# 		return dicalCmd


# 	def returnMLE (self):
# 		# get the MLE
# 		(mlePoint, maxLike) = returnMLE (self.diCalOutputFileName)

# 		# convert it to real scales
# 		realMlePoint = ([ExpRateToPerGenPercent (self.Nref, mlePoint[0])] +
# 			[CoalTimeToGenTime (self.Nref, mlePoint[1])] +
# 			[CoalSizeToDiploidSize (self.Nref, x) for x in mlePoint[2:]])

# 		# and return it
# 		return (realMlePoint, maxLike)


# class ExpGrowthFixedOnsetTimeAnalysis:

# 	metaNumStartPoints = 10
# 	metaNumPoints = 5
# 	metaKeepBest = 2
# 	metaNumIterations = 3
# 	numEMiterations = 3
# 	numIterationsMstep = 3

# 	cl = "lol"

# 	numLociPerHmmStep = 500

# 	Nref = 10000
# 	mu = 1.25e-8
# 	r = 1.25e-8

# 	numEpochs = 3
# 	firstFixedTimeInGen = 1000
# 	lastFixedTimeInGen = 3000
# 	# fixedTimeInGen = 2000

# 	# bounds related stuff
# 	minPerGenGrowth = 0.0001
# 	maxPerGenGrowth = 2
# 	minBoundTimeInGen = 100
# 	maxBoundTimeInGen = 1000
# 	maxSize = 100000
# 	minSize = 100

# 	jarFile = "cleanDiCal2.jar"


# 	def __init__ (self, uniqueBasename, numCores, vcfFiles, refFiles, sampleSize, randomSeed):
# 		self.vcfFileList = vcfFiles
# 		self.refFileList = refFiles
# 		self.seeed = randomSeed
# 		self.numCores = numCores

# 		self.diCalOutputFileName = "%s.out" % uniqueBasename
# 		self.paramFile = "%s.param" % uniqueBasename
# 		self.demoFile = "%s.demo" % uniqueBasename
# 		self.rateFile = "%s.rates" % uniqueBasename
# 		self.configFile = "%s.config" % uniqueBasename

# 		# put together the bounds
# 		# one parameter per epoch
# 		self.bounds = []
# 		self.bounds.extend ([[PerGenPercentToExpRate (self.Nref, self.minPerGenGrowth), PerGenPercentToExpRate (self.Nref, self.maxPerGenGrowth)]])
# 		# self.bounds.extend ([[GenTimeToCoalTime (self.Nref, self.minBoundTimeInGen), GenTimeToCoalTime (self.Nref, self.maxBoundTimeInGen)]])
# 		# self.bounds.extend ([[GenTimeToCoalTime (self.Nref, 10*self.minBoundTimeInGen), GenTimeToCoalTime (self.Nref, 10*self.maxBoundTimeInGen)]])
# 		self.bounds.extend ([[DiploidSizeToCoalSize (self.Nref, self.minSize), DiploidSizeToCoalSize (self.Nref, self.maxSize)]]*self.numEpochs)

# 		# random starting points
# 		self.metaStartFile = "%s.rand" % uniqueBasename
# 		samplingRegionFactor = 0.25
# 		samplingRegion = [[math.exp (math.log(bs[0]) + samplingRegionFactor * (math.log(bs[1]) - math.log(bs[0]))), math.exp (math.log(bs[0]) + (1-samplingRegionFactor) * (math.log(bs[1]) - math.log(bs[0])))] for bs in self.bounds]
# 		logUniformStartingPointsFile (self.metaStartFile, self.metaNumStartPoints, samplingRegion)

# 		# param file
# 		writeMutRecoParameterFile (self.paramFile, theta (self.Nref, self.mu), rho (self.Nref, self.r))

# 		# config file
# 		# numLoci doesn't even matter, also, two alleles and one deme
# 		writeConfigFile (self.configFile, 123456, 2, 1, [sampleSize])

# 		# demography file
# 		# firstTime = GenTimeToCoalTime (self.Nref, self.firsTimeInGen)
# 		# lastTime = GenTimeToCoalTime (self.Nref, self.lastTimeInGen)
# 		# themTimes = [math.exp(math.log(firstTime) + x/float(self.numEpochs-2) * (math.log(lastTime) - math.log(firstTime))) for x in range(self.numEpochs-1)]
# 		self.epochTimes = logGrid (GenTimeToCoalTime (self.Nref, self.firstFixedTimeInGen), GenTimeToCoalTime (self.Nref, self.lastFixedTimeInGen), self.numEpochs-1)
# 		writePieceWiseConstantDemographyFile (self.demoFile, self.numEpochs, self.epochTimes, 1)
# 		writeOnePopExpRateRecentFile (self.rateFile, self.numEpochs)


# 	def run (self):

# 		if (self.numCores > 1):
# 			parallelEmSteps = math.ceil (self.numCores / 2)
# 			parallelString = " --parallel %d --metaParallelEmSteps %d" % (self.numCores, parallelEmSteps)
# 		else:
# 			parallelString = ""

# 		dicalCmd = ("java -Xmx70G -jar %s" +
# 			" --paramFile %s" +
# 			" --demoFile %s" +
# 			" --ratesFile %s" +
# 			" --configFile %s" +
# 			" --vcfFile %s" +
# 			" --vcfFilterPassString 'PASS'" +
# 			" --vcfReferenceFile %s" +
# 			" --seed %d" +
# 			" --lociPerHmmStep %d" +
# 			" --compositeLikelihood %s" +
# 			" --metaStartFile %s" +
# 			" --metaNumIterations %d" +
# 			" --metaKeepBest %d" +
# 			" --metaNumPoints %d" +
# 			" --numberIterationsEM %d" +
# 			" --numberIterationsMstep %d" +
# 			" --printEmPath" +
# 			# " --intervalType simple" +
# 			" --intervalType logUniform" +
# 			" --intervalParams '11,0.01,4'" +
# 			" --bounds '%s'" +
# 			"%s" +
# 			" >%s") % (
# 			self.jarFile,
# 			self.paramFile,
# 			self.demoFile,
# 			self.rateFile,
# 			self.configFile, 
# 			self.vcfFileList,
# 			self.refFileList,
# 			self.seeed,
# 			self.numLociPerHmmStep,
# 			self.cl,
# 			self.metaStartFile,
# 			self.metaNumIterations,
# 			self.metaKeepBest,
# 			self.metaNumPoints,
# 			self.numEMiterations,
# 			self.numIterationsMstep,
# 			";".join([",".join([str(x) for x in y]) for y in self.bounds]),
# 			parallelString,
# 			self.diCalOutputFileName)

# 		return dicalCmd


# 	def returnMLE (self):
# 		# get the MLE
# 		(mlePoint, maxLike) = returnMLE (self.diCalOutputFileName)

# 		# convert it to real scales
# 		realMlePoint = ([ExpRateToPerGenPercent (self.Nref, mlePoint[0])] +
# 			[CoalSizeToDiploidSize (self.Nref, x) for x in mlePoint[1:]])

# 		# and return it
# 		return (realMlePoint, maxLike)


class PieceWiseConstantAnalysis:

	# metaNumStartPoints = 40
	# metaNumPoints = 12
	# metaKeepBest = 4
	# metaNumIterations = 3
	numEMiterations = 10
	numIterationsMstep = 3

	cl = "lol"

	numLociPerHmmStep = 1000

	Nref = 10000
	mu = 1.25e-8
	r = 1.25e-8

	numEpochs = 12
	firsTimeInGen = 200
	lastTimeInGen = 20000
	maxSize = 200000
	minSize = 100

	yearsPerGen = 25

	jarFile = "cleanDiCal2.jar"


	def __init__ (self, uniqueBasename, numCores, vcfFiles, refFiles, sampleSize, randomSeed):
		self.vcfFileList = vcfFiles
		self.refFileList = refFiles
		self.seeed = randomSeed
		self.numCores = numCores

		self.diCalOutputFileName = "%s.dical_out" % uniqueBasename
		self.paramFile = "%s.param" % uniqueBasename
		self.demoFile = "%s.demo" % uniqueBasename
		self.configFile = "%s.config" % uniqueBasename

		# put together the bounds
		# one parameter per epoch
		self.bounds = [[DiploidSizeToCoalSize (self.Nref, self.minSize), DiploidSizeToCoalSize (self.Nref, self.maxSize)]]*self.numEpochs

		# random starting point(s)
		# self.metaStartFile = "%s.rand" % uniqueBasename
		samplingRegionFactor = 0.25
		samplingRegion = [[math.exp (math.log(bs[0]) + samplingRegionFactor * (math.log(bs[1]) - math.log(bs[0]))), math.exp (math.log(bs[0]) + (1-samplingRegionFactor) * (math.log(bs[1]) - math.log(bs[0])))] for bs in self.bounds]
		# logUniformStartingPointsFile (self.metaStartFile, self.metaNumStartPoints, samplingRegion)
		logBounds = [[math.log(b[0]), math.log(b[1])] for b in samplingRegion]
		# and get a point
		logPoint = [random.uniform(b[0], b[1]) for b in logBounds]
		point = [math.exp(x) for x in logPoint]
		point = [1]*self.numEpochs
		self.startPointString = "'" + ",".join([("%.8f" % x) for x in point]) + "'"


		# param file
		writeMutRecoParameterFile (self.paramFile, theta (self.Nref, self.mu), rho (self.Nref, self.r))

		# config file
		# numLoci doesn't even matter, also, two alleles and one deme
		writeConfigFile (self.configFile, 123456, 2, 1, [sampleSize])

		# demography file
		firstTime = GenTimeToCoalTime (self.Nref, self.firsTimeInGen)
		lastTime = GenTimeToCoalTime (self.Nref, self.lastTimeInGen)
		self.epochTimes = logGrid (firstTime, lastTime, self.numEpochs-1)
		writePieceWiseConstantDemographyFile (self.demoFile, self.numEpochs, self.epochTimes, 0)


	def run (self):

		if (self.numCores > 1):
			parallelEmSteps = min (math.ceil (self.numCores / 2), 2)
			# parallelString = " --parallel %d --metaParallelEmSteps %d" % (self.numCores, parallelEmSteps)
			parallelString = " --parallel %d" % (self.numCores)
		else:
			parallelString = ""

		dicalCmd = ("java -Xmx15G -jar %s" +
			" --paramFile %s" +
			" --demoFile %s" +
			" --configFile %s" +
			" --vcfFile %s" +
			" --vcfFilterPassString 'PASS'" +
			" --vcfReferenceFile %s" +
			" --seed %d" +
			" --lociPerHmmStep %d" +
			" --compositeLikelihood %s" +
			" --startPoint %s" +
			# " --metaStartFile %s" +
			# " --metaNumIterations %d" +
			# " --metaKeepBest %d" +
			# " --metaNumPoints %d" +
			" --numberIterationsEM %d" +
			# " --disableCoordinateWiseMStep" +
			" --numberIterationsMstep %d" +
			" --printEmPath" +
			" --intervalType simple" +
			# " --intervalType logUniform" +
			# " --intervalParams '11,0.01,4'" +
			" --bounds '%s'" +
			"%s" +
			" >%s") % (
			self.jarFile,
			self.paramFile,
			self.demoFile,
			self.configFile, 
			self.vcfFileList,
			self.refFileList,
			self.seeed,
			self.numLociPerHmmStep,
			self.cl,
			self.startPointString,
			# self.metaStartFile,
			# self.metaNumIterations,
			# self.metaKeepBest,
			# self.metaNumPoints,
			self.numEMiterations,
			self.numIterationsMstep,
			";".join([",".join([str(x) for x in y]) for y in self.bounds]),
			parallelString,
			self.diCalOutputFileName)

		return dicalCmd


	def returnMLE (self):
		# get the MLE
		(mlePoint, maxLike) = returnMLE (self.diCalOutputFileName)

		# convert it to real scales
		realMlePoint = [CoalSizeToDiploidSize (self.Nref, x) for x in mlePoint]

		# and return it
		return (realMlePoint, maxLike)


	def writeResultsCSV (self, filename):
		# get the MLE
		(mlePoint, maxLike) = self.returnMLE ()

		# also get them times rescaled to years
		epochTimesInGen = [0] + [CoalTimeToGenTime (self.Nref, x) * self.yearsPerGen for x in self.epochTimes]

		assert (len(mlePoint) == len(epochTimesInGen))
		# now write it to file in hopefully right format
		ofs = open (filename, "w")
		for i in range(len(mlePoint)):
			ofs.write ("%.3f\t%.3f\n" % (epochTimesInGen[i], mlePoint[i]))
		ofs.close()
		# done


def realAnalysis ():

	# -------------------- piece-wise constant --------------------
	# set the directory
	laptopBase = "/Users/steinrue/labsharecri/"
	realBase = "/home/steinrue/labshare"
	dirExtension = "projects/coalHMMopionPiece/analysis/test/pieceWiseTestThree"
	daDir = os.path.join (laptopBase, dirExtension)
	# daDir = "/Users/steinrue/labsharecri/projects/coalHMMopionPiece/preLuigiStuff/pieceTest"
	# os.makedirs (daDir)
	os.chdir (daDir)

	uniqueBasename = "catfishFive"

	dataDir = "/home/steinrue/labshare/projects/coalHMMopionPiece/data/test/increase"

	numContigs = 10

	# # make reference
	# numLoci = 50000
	# refFilename = "%s.ref" % uniqueBasename
	# makeReferenceFile (refFilename, numLoci)

	# we have a reference
	refFilename = os.path.join (dataDir, "output.ref")
	refList = [refFilename]*numContigs
	refFiles = "'" + ",".join(refList) + "'"
	

	vcfList = []
	for d in range(numContigs):
		vcfList.append (os.path.join (dataDir, "contig.%d.vcf" % d))
	vcfFiles = "'" + ",".join(vcfList) + "'"

	numCores = 24

	# try to put some analysis together
	# arguments: <fileBasename> <numCorse> <list_of_vcf> <list_of_refFiles> <sample_size> <seed>
	analysis = PieceWiseConstantAnalysis (uniqueBasename, numCores, vcfFiles, refFiles, 10, 4712)
	print ("---------- COMMAND ----------")
	cmd = analysis.run()
	print (cmd)
	ofs = open ("%s.sh" % uniqueBasename, "w")
	ofs.write("#!/bin/bash\n")
	ofs.write("#PBS -N %s\n" % uniqueBasename)
	ofs.write("#PBS -S /bin/bash\n")
	ofs.write("#PBS -l walltime=10:00:00\n")
	ofs.write("#PBS -l nodes=1:ppn=%d\n" % numCores)
	ofs.write("#PBS -l mem=16gb\n")
	ofs.write("#PBS -o %s\n" % os.path.join (realBase, dirExtension, "%s.out" % uniqueBasename))
	ofs.write("#PBS -e %s\n" % os.path.join (realBase, dirExtension, "%s.err" % uniqueBasename))
	# try to change the workdir
	ofs.write("cd $PBS_O_WORKDIR\n")
	ofs.write (cmd + "\n")
	ofs.close()
	# print (analysis.run())
	print ("---------- MLE ----------")
	# print (analysis.returnMLE())
	analysis.writeResultsCSV ("%s.csv" % uniqueBasename)

	# # -------------------- exponential growth, estimate onset --------------------
	# # set the directory
	# daDir = "/Users/steinrue/labsharecri/projects/coalHMMopionPiece/preLuigiStuff/expTestEstimateTime"
	# # os.makedirs (daDir)
	# os.chdir (daDir)

	# uniqueBasename = "tortoise"

	# # make reference
	# numLoci = 50000
	# refFilename = "%s.ref" % uniqueBasename
	# makeReferenceFile (refFilename, numLoci)

	# vcfFilename = "%s.vcf" % uniqueBasename

	# vcfFiles = "'%s,%s,%s'" % (vcfFilename, vcfFilename, vcfFilename)
	# refFiles = "'%s,%s,%s'" % (refFilename, refFilename, refFilename)

	# # try to put some analysis together
	# # arguments: <fileBasename> <numCorse> <list_of_vcf> <list_of_refFiles> <sample_size> <seed>
	# analysis = ExpGrowthEstimateOnsetTimeAnalysis (uniqueBasename, 3, vcfFiles, refFiles, 12, 4711)
	# print ("---------- COMMAND ----------")
	# print (analysis.run())
	# print ("---------- MLE ----------")
	# # print (analysis.returnMLE())

	# # -------------------- exponential growth, fixed onset --------------------
	# # set the directory
	# daDir = "/Users/steinrue/labsharecri/projects/coalHMMopionPiece/preLuigiStuff/expTestFixed"
	# # os.makedirs (daDir)
	# os.chdir (daDir)

	# uniqueBasename = "barnacle"

	# # make reference
	# numLoci = 50000
	# refFilename = "%s.ref" % uniqueBasename
	# makeReferenceFile (refFilename, numLoci)

	# vcfFilename = "%s.vcf" % uniqueBasename

	# vcfFiles = "'%s,%s,%s'" % (vcfFilename, vcfFilename, vcfFilename)
	# refFiles = "'%s,%s,%s'" % (refFilename, refFilename, refFilename)

	# # try to put some analysis together
	# # arguments: <fileBasename> <numCorse> <list_of_vcf> <list_of_refFiles> <sample_size> <seed>
	# analysis = ExpGrowthFixedOnsetTimeAnalysis (uniqueBasename, 3, vcfFiles, refFiles, 12, 4711)
	# print ("---------- COMMAND ----------")
	# print (analysis.run())
	# print ("---------- MLE ----------")
	# # print (analysis.returnMLE())


def test():

	# writePieceWiseConstantDemographyFile ("constFixed.demo", 4, [0.01,0.1,1], 1)
	# writePieceWiseConstantDemographyFile ("constPartial.demo", 4, [0.1,1], 1)
	# writePieceWiseConstantDemographyFile ("constEstimate.demo", 4, [], 1)

	# print (logGrid (1,100,4))
	# print (logGrid (1,100,1))
	pass

def plotStuff ():
	Nref = 10000
	(mlePoint, maxLike) = returnMLE ("catfishFive.dical_out")

	# print (mlePoint)

	# convert it to real scales
	realMlePoint = [CoalSizeToDiploidSize (Nref, x) for x in mlePoint]

	# plot some stuff
	plt.clf()
	# first the estimate
	firsTimeInGen = 200
	lastTimeInGen = 20000
	numEpochs = 12
	# needs some random modifications to works, cuase with 'post', the step-function would jsut omit the last value
	xs = [100] + logGrid (firsTimeInGen, lastTimeInGen, numEpochs-1) + [40000]
	# print (xs)
	# print (realMlePoint)
	realMlePoint =  realMlePoint + [realMlePoint[-1]]
	plt.step (xs, realMlePoint, where='post')
	ax = plt.gca()
	ax.set_xscale('log')
	# also plot the truth
		# pop sizes
	N_1 = 10000
	N_2 = 20000
	N_A = 10000
	# times
	# generation_time = 25
	tOne = 1000
	tTwo = 2000
	xs = [100] + [1000,2000] + [40000]
	ys = [10000,20000,10000] + [10000]
	plt.step (xs, ys, "r",where='post')

	plt.savefig ("popSize.pdf")


def main():

	# realAnalysis ()

	# test()

	plotStuff()







if __name__ == "__main__":
	main()


