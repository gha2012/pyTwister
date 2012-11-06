from pymol import cmd
from pymol import stored
import numpy
import sys
import locale
import random

locale.setlocale(locale.LC_NUMERIC, "")
out = sys.stdout
def __init__(self):
	"""pyTwister: A PyMOL version of Twister"""

def format_num(num):
	"""Format a number according to given places.
	Adds commas, etc."""

	try:
		return locale.format("%.2f", (num), True)

	except (ValueError, TypeError):
		return str(num)

def get_max_width(table, index):
	"""Get the maximum width of the given column index"""
	return max([len(format_num(row[index])) for row in table])

def pprint_table(out, table):
	"""Prints out a table of data, padded for alignment
	@param out: Output stream (file-like object)
	@param table: The table to print. A list of lists.
	Each row must have the same number of columns. """
	col_paddings = []

	for i in range(len(table[0])):
		col_paddings.append(get_max_width(table, i))

	for row in table:
		# left col
		print >> out, row[0].ljust(col_paddings[0] + 1),
		# rest of the cols
		for i in range(1, len(row)):
			col = format_num(row[i]).rjust(col_paddings[i] + 2)
			print >> out, col,
		print >> out		

def getDihedral(a,b,c,d):
	v1 = getNormedVector(a, b)
	v2 = getNormedVector(b, c)
	v3 = getNormedVector(c, d)
	v1v2 = numpy.cross(v1,v2)
	v2v3 = numpy.cross(v2,v3)
	return getAngle(v1v2,v2v3)

def getNormedVector(a,b):
	return (b-a)/numpy.linalg.norm(b-a)

def getAngle(a,b):
	return numpy.rad2deg(numpy.arccos(numpy.dot(a/numpy.linalg.norm(a),b/numpy.linalg.norm(b))))

def createPseudoatom (coordinates, objectName):
	x=float(coordinates[0])
	y=float(coordinates[1])
	z=float(coordinates[2])
	#print x,y,z,objectName
	posString="[%3.2f,%3.2f,%3.2f]" % (x,y,z)
	cmd.pseudoatom(pos=posString, object=objectName)

def getCaPositions(selection):	 
	#print selection
	stored.xyz = []
	cmd.iterate_state(1,selection,"stored.xyz.append([x,y,z])")
	caPositions=numpy.array(stored.xyz)
	return caPositions

def getAx(caPositions):
	#Sergej Strelkovs way of calculating the axis positions
	axisPositions = numpy.zeros(shape=(len(caPositions),3))
	for i in range (1, len(caPositions)-2):
		m=(caPositions[i-1]+caPositions[i+1])/2-caPositions[i]
		m1=(caPositions[i]+caPositions[i+2])/2-caPositions[i+1]
		d=caPositions[i+1]-caPositions[i]	
		mu=numpy.dot(d,m)/(numpy.dot(m,m)-numpy.dot(m,m1))
		nu=numpy.dot(d,m1)/(numpy.dot(m,m1)-numpy.dot(m1,m1))
		
		currentAxisPosition=caPositions[i]+m*mu
		
		if i == 1:
			axisPositions[1]=currentAxisPosition
		else:
			axisPositions[i]=(axisPositions[i]+currentAxisPosition)/2
		#insert for next cycle, will be avgd with currentAxisPosition
		axisPositions[i+1]=caPositions[i+1]+m1*nu
	
	#extrapolate for first and last residue
	axisPositions[0]=2*axisPositions[1]-axisPositions[2]
	axisPositions[len(axisPositions)-1]=2*axisPositions[len(axisPositions)-2]-axisPositions[len(axisPositions)-3]
	#print axisPositions
	return axisPositions

def calculateRadii(caPositions, axisPoints):
	radii = numpy.zeros(shape=(len(caPositions),1))
	for i in range (0,len(axisPoints)):
		radius=numpy.linalg.norm(caPositions[i]-axisPoints[i])
		radii[i]=radius
	return radii

def calculateRisesPerResidue(axisPoints):
	risesPerResidue = numpy.zeros(shape=(len(axisPoints),1))
	for i in range (1,len(axisPoints)-1):
		rise=(numpy.linalg.norm(axisPoints[i-1]-axisPoints[i])+numpy.linalg.norm(axisPoints[i]-axisPoints[i+1]))/2
		risesPerResidue[i] = rise
	return risesPerResidue

def calculatePhaseYieldsPerResidue(twistingAxis, fixedAxis):
	phaseYieldsPerResidue = numpy.zeros(shape=(len(twistingAxis),1))
	for i in range (1,len(twistingAxis)-1):
		a=getDihedral(twistingAxis[i-1], fixedAxis[i-1], fixedAxis[i], twistingAxis[i])
		b=getDihedral(twistingAxis[i], fixedAxis[i], fixedAxis[i+1], twistingAxis[i+1])
		phaseYieldsPerResidue[i] = ((a+b)/2)
	return phaseYieldsPerResidue

def calculatePitchesPerResidue(axisPoints,rise,phase):
	pitchesPerResidue = numpy.zeros(shape=(len(axisPoints),1))
	for i in range (1,len(axisPoints)-1):
		pitch=rise[i]*360/phase[i]
		pitchesPerResidue[i] = pitch
	return pitchesPerResidue
	
def calculateResiduesPerTurn(axisPoints,rise,phase):
	residuesPerTurn = numpy.zeros(shape=(len(axisPoints),1))
	for i in range (1,len(axisPoints)-1):
		residuesPerTurn[i] = 360/phase[i]
	return residuesPerTurn

def calculateCrickAngle(An, On, Cn, nextOn):
	OC=getNormedVector(On, Cn)
	OA=getNormedVector(On, An)
	OnextOn=getNormedVector(On, nextOn)
	cross=numpy.cross(OC,OnextOn)
	mixed=numpy.dot(OA, cross)
	angle=getAngle(OC,OA)
	if (mixed<0):
		angle*=1
	else: #(mixed>=0):
		angle*=-1
	return angle

def calculateavgCrickAngles(helicesInSelection, avgCoiledCoilParameters, chains):
	#calculate Crick angles avgd over all helices and assign positions
	crickAnglesForAllHelices=[]
	for i in range (0,len(chains)):
		#print i
		crickAngles=[]
		crickAngles.append(0)
		for j in range (1,len(helicesInSelection[i]["avgAlphaHelixAxisPoints"])-1):
			angle=calculateCrickAngle(helicesInSelection[i]['caPositions'][j], helicesInSelection[i]["avgAlphaHelixAxisPoints"][j], avgCoiledCoilParameters['coiledCoilAxisPoints'][j], helicesInSelection[i]["avgAlphaHelixAxisPoints"][j+1])
			#print angle
			crickAngles.append(angle)
		crickAnglesForAllHelices.append(crickAngles)
	
	avgCrickAngles=[]
	for i in range (0, len(crickAngles)):
		avg=0
		x=0
		y=0
		for j in range(0, len(chains)):
			#print crickAnglesForAllHelices[j][i]
			x+=numpy.cos(crickAnglesForAllHelices[j][i]*numpy.pi/180)
			y+=numpy.sin(crickAnglesForAllHelices[j][i]*numpy.pi/180)
		
		avg=numpy.arctan2(y,x)*180/numpy.pi
		#print avg
		avgCrickAngles.append(avg)
	avgCrickAngles.append(0)
	#print "avgCrickAngles:", avgCrickAngles
	return avgCrickAngles

def assignHeptadPositions(helicesInSelection, avgCrickAngles, residueNumbersInSelection, selection):
	pos=[]
	pos.append("?")
	for i in range (1,len(helicesInSelection[0]["avgAlphaHelixAxisPoints"])):
		#print avgCrickAngles[i]
		if (i > 1 and avgCrickAngles[i-1] < 0 and avgCrickAngles[i] > 0 and avgCrickAngles[i-1] + avgCrickAngles[i] < 0):
			#print avgCrickAngles[i]
			pos.append("a")
			cmd.color("red", selection+" & resi "+str(residueNumbersInSelection[i]))
		elif (avgCrickAngles[i] < 0 and avgCrickAngles[i+1] > 0 and (avgCrickAngles[i] + avgCrickAngles[i+1]) > 0):
			pos.append("d")
			cmd.color("blue", selection+" & resi "+str(residueNumbersInSelection[i]))
		else:
			pos.append("?")
	pos.append("?")
	#print pos
	for i in range (1, len(pos)-1):
		if pos[i]=="?":
			if pos[i-1] == "a":
				pos[i]="b"
			elif pos[i-1] == "b":
				pos[i]="c"
			elif pos[i-1] == "d":
				pos[i]="e"
			elif pos[i-1] == "e":
				pos[i]="f"
			elif pos[i-1] == "f":
				pos[i]="g"

	#print pos
	#tidy up beginning
	for i in range (4, 0, -1):
		if pos[i]=="?":
			if pos[i+1] == "a":
				pos[i]="g"
			if pos[i+1] == "b":
				pos[i]="a"
			if pos[i+1] == "c":
				pos[i]="b"
			if pos[i+1] == "d":
				pos[i]="c"
			if pos[i+1] == "e":
				pos[i]="d"
			if pos[i+1] == "f":
				pos[i]="e"
			if pos[i+1] == "g":
				pos[i]="f"
	return pos

def makeResultsTable(helicesInSelection, residuesInSelection, residueNumbersInSelection, avgCoiledCoilParameters, avgAlphaHelixParameters, heptadPositions, avgCrickAngles):
	table = []
	table = [["Res.", "cc-rad", "cc-rise", "cc-pitch", "cc-phaseYield", "pos","Crick-angle", "a-radius", "a-rise", "a-res/turn", "a-phaseYield"]]
	for i in range (1,len(helicesInSelection[0]["avgAlphaHelixAxisPoints"])-1):
		#print i
		list=[]
		list.append(stored.residuesInSelection[i]+ " " +stored.residueNumbersInSelection[i])
		list.append(avgCoiledCoilParameters['coiledCoilRadiusPerResidue'][i])
		list.append(avgCoiledCoilParameters['coiledCoilRisesPerResidue'][i])
		list.append(avgCoiledCoilParameters['coiledCoilPitch'][i])
		list.append(avgCoiledCoilParameters['coiledCoilPhaseYieldsPerResidue'][i])
		list.append(heptadPositions[i])
		list.append(avgCrickAngles[i])
		list.append(avgAlphaHelixParameters['alphaHelixRadiusPerResidue'][i])
		list.append(avgAlphaHelixParameters['alphaHelixRisesPerResidue'][i])
		list.append(avgAlphaHelixParameters['alphaHelixResiduesPerTurn'][i])
		list.append(avgAlphaHelixParameters['alphaHelixPhaseYieldsPerResidue'][i])
		table.append(list)
	
	avgs = []
	stdevs = []
	avgCoiledCoilRadius = numpy.mean(avgCoiledCoilParameters['coiledCoilRadiusPerResidue'])
	stdevCoiledCoilRadius = numpy.std(avgCoiledCoilParameters['coiledCoilRadiusPerResidue'])
	
	avgCoiledCoilRise = numpy.mean(avgCoiledCoilParameters['coiledCoilRisesPerResidue'])
	stdevCoiledCoilRise = numpy.std(avgCoiledCoilParameters['coiledCoilRisesPerResidue'])
	
	avgCoiledCoilPitch = numpy.mean(avgCoiledCoilParameters['coiledCoilPitch'])
	stdevCoiledCoilPitch = numpy.std(avgCoiledCoilParameters['coiledCoilPitch'])
	
	avgCoiledCoilPhaseYield = numpy.mean(avgCoiledCoilParameters['coiledCoilPhaseYieldsPerResidue'])
	stdevCoiledCoilPhaseYield = numpy.std(avgCoiledCoilParameters['coiledCoilPhaseYieldsPerResidue'])
	
	avgAlphaHelixRadius = numpy.mean(avgAlphaHelixParameters['alphaHelixRadiusPerResidue'])
	stdevAlphaHelixRadius = numpy.std(avgAlphaHelixParameters['alphaHelixRadiusPerResidue'])
	
	avgAlphaHelixRise = numpy.mean(avgAlphaHelixParameters['alphaHelixRisesPerResidue'])
	stdevAlphaHelixRise = numpy.std(avgAlphaHelixParameters['alphaHelixRisesPerResidue'])
	
	avgResPerTurn = numpy.mean(avgAlphaHelixParameters['alphaHelixResiduesPerTurn'])
	stdevResPerTurn = numpy.std(avgAlphaHelixParameters['alphaHelixResiduesPerTurn'])
	
	avgAlphaHelixPhaseYield = numpy.mean(avgAlphaHelixParameters['alphaHelixPhaseYieldsPerResidue'])
	stdevAlphaHelixPhaseYield = numpy.std(avgAlphaHelixParameters['alphaHelixPhaseYieldsPerResidue'])
	
	avgs = ["Avg", avgCoiledCoilRadius, avgCoiledCoilRise, avgCoiledCoilPitch, avgCoiledCoilPhaseYield, "-","-", avgAlphaHelixRadius, avgAlphaHelixRise, avgResPerTurn, avgAlphaHelixPhaseYield]
	stdevs = ["Std", stdevCoiledCoilRadius, stdevCoiledCoilRise, stdevCoiledCoilPitch, stdevCoiledCoilPhaseYield, "-","-", stdevAlphaHelixRadius, stdevAlphaHelixRise, stdevResPerTurn, stdevAlphaHelixPhaseYield]
	empty = ["","","","","","","","","","",""]
	
	table.append(empty)
	table.append(avgs)
	table.append(stdevs)
	return table

def getOLC(string):
	if string == "GLY": return "G"
	elif string == "ALA": return "A"
	elif string == "VAL": return "V"
	elif string == "LEU": return "L"
	elif string == "ILE": return "I"
	elif string == "PRO": return "P"
	elif string == "ASN": return "N"
	elif string == "ASP": return "D"
	elif string == "GLN": return "Q"
	elif string == "GLU": return "E"
	elif string == "HIS": return "H"
	elif string == "SER": return "S"
	elif string == "THR": return "T"
	elif string == "ARG": return "R"
	elif string == "LYS": return "K"
	elif string == "TRP": return "W"
	elif string == "TYR": return "Y"
	elif string == "PHE": return "F"
	elif string == "CYS": return "C"
	elif string == "MET": return "M"
	else: return "?"

def printSequence(residueNumbersInSelection, residuesInSelection, heptadPositions):
	charactersPerLine = 60
	lines=[]
	numberOfResidues = int(residueNumbersInSelection[len(residueNumbersInSelection)-1])-int(residueNumbersInSelection[0])+1
	firstResidue = int(residueNumbersInSelection[0])
	#print residueNumbersInSelection[len(residueNumbersInSelection)-1]
	#print residueNumbersInSelection[0]
	#print len(heptadPositions), len(residueNumbersInSelection), numberOfResidues
	h = 0
	while h < numberOfResidues:
		line1=""
		line2=""
		line3=""
		for i in range (0,3):
			j = 0
			if numberOfResidues - h >= charactersPerLine:
				charactersPerLine = charactersPerLine
			else:
				charactersPerLine = numberOfResidues - h
			while j < charactersPerLine:
				if i == 0 and (firstResidue+h+j) % 10 == 0:
					#print firstResidue+h+j
					line1+=str(firstResidue+h+j)
					j+=len(str(firstResidue+h+j))-1
				elif i == 0 and (firstResidue+h+j) % 10 > len(str(firstResidue+h+j))-1:
					line1+=" "
				if i == 1:
					line2+=getOLC(residuesInSelection[h+j])
				if i == 2:
					line3+=heptadPositions[h+j]
				j+=1
		lines.append(line1)
		lines.append(line2)
		lines.append(line3)
		lines.append(" ")
		h+=charactersPerLine
	print ""
	print "Sequence with heptad positions:"
	print ""
	for line in lines:
		print line

def showAxesInPymol(helicesInSelection, avgCoiledCoilParameters):
	cmd.set("dash_gap","0")
	cmd.set("dash_radius","0.2")
	for i in range (0, len(helicesInSelection)):
		#print len(helicesInSelection),i
		for j in range (0,len(helicesInSelection[i]["avgAlphaHelixAxisPoints"])):
			createPseudoatom(helicesInSelection[i]["avgAlphaHelixAxisPoints"][j],"helixAxisPoint_"+str(i)+"_"+str(j))
			createPseudoatom(avgCoiledCoilParameters["coiledCoilAxisPoints"][j],"coiledCoilAxisPoint_"+str(i)+"_"+str(j))
			if j > 0:
				cmd.distance("AxisLine_"+str(i)+"_"+str(j),"helixAxisPoint_"+str(i)+"_"+str(j-1),"helixAxisPoint_"+str(i)+"_"+str(j))
				cmd.distance("CoiledCoilAxisLine_"+str(i)+"_"+str(j),"coiledCoilAxisPoint_"+str(i)+"_"+str(j-1),"coiledCoilAxisPoint_"+str(i)+"_"+str(j))
		cmd.group("AlphaAxis_"+str(i), "AxisLine_"+str(i)+"*")
		cmd.group("AlphaAxisPoints_"+str(i), "helixAxisPoint_"+str(i)+"*")
		cmd.group("coiledCoil", "coiledCoil*")
		cmd.set("dash_color","yellow", "Alpha*")
		cmd.set("dash_color","red", "Coiled*")
		#print i
	cmd.hide("labels")
	cmd.hide("nonbonded")
 
def pyTwister(selection, chains):
	#make list of residueNames
	stored.residueNumbersInSelection=[]
	cmd.iterate(selection + " & name ca", "stored.residueNumbersInSelection.append(resi)")
	stored.residuesInSelection=[]
	cmd.iterate(selection + " & name ca", "stored.residuesInSelection.append(resn)")
	stored.residueOneLetterCodes=[]
	
	#list of helices with Parameters
	helicesInSelection=[]
	
	#calculate alphaHelixParameters
	for i in range (0, len(chains)):
		dict={}
		#get ca positions
		dict['caPositions']=getCaPositions(selection + " & chain " + chains[i] + " & name ca")
		dict['avgAlphaHelixAxisPoints']=getAx(dict['caPositions'])
		dict['alphaHelixRadiusPerResidue']=calculateRadii(dict['caPositions'], dict['avgAlphaHelixAxisPoints'])
		dict['alphaHelixRisesPerResidue']=calculateRisesPerResidue(dict['avgAlphaHelixAxisPoints'])
		dict['alphaHelixPhaseYieldsPerResidue']=calculatePhaseYieldsPerResidue(dict['caPositions'], dict['avgAlphaHelixAxisPoints'])
		dict['alphaHelixResiduesPerTurn']=calculateResiduesPerTurn(dict['avgAlphaHelixAxisPoints'],dict['alphaHelixRisesPerResidue'],dict['alphaHelixPhaseYieldsPerResidue'])
		helicesInSelection.append(dict)
   
	#calculate avg of alpha Helix parameters - avg of AxisPoints is coiledCoilAxis
	avgAlphaHelixParameters={}
	for key in ['avgAlphaHelixAxisPoints','alphaHelixRadiusPerResidue', 'alphaHelixRisesPerResidue', 'alphaHelixPhaseYieldsPerResidue', 'alphaHelixResiduesPerTurn']:
		avg=0
		i=0
		while (i < len(helicesInSelection)):
			avg+=numpy.array(helicesInSelection[i][key])
			i+=1
		avg/=len(helicesInSelection)
		avgAlphaHelixParameters[key]=avg
	
	#calculate coiled coil parameters
	coiledCoils=[]
	for i in range (0, len(chains)):
		dict={}
		dict['coiledCoilAxisPoints']=avgAlphaHelixParameters['avgAlphaHelixAxisPoints']
		dict['coiledCoilRadiusPerResidue']=calculateRadii(helicesInSelection[i]['avgAlphaHelixAxisPoints'], dict['coiledCoilAxisPoints'])
		dict['coiledCoilRisesPerResidue']=calculateRisesPerResidue(dict['coiledCoilAxisPoints'])
		dict['coiledCoilPhaseYieldsPerResidue']=-calculatePhaseYieldsPerResidue(helicesInSelection[i]['avgAlphaHelixAxisPoints'], dict['coiledCoilAxisPoints'])
		dict['coiledCoilPitch']=-calculatePitchesPerResidue(dict['coiledCoilAxisPoints'],dict['coiledCoilRisesPerResidue'],dict['coiledCoilPhaseYieldsPerResidue'])
		coiledCoils.append(dict)
	
	#calculate avg of coiled coil parameters
	avgCoiledCoilParameters={}
	for key in ['coiledCoilAxisPoints','coiledCoilRadiusPerResidue', 'coiledCoilRisesPerResidue', 'coiledCoilPhaseYieldsPerResidue', 'coiledCoilPitch']:
		avg=0
		i=0
		while (i < len(coiledCoils)):
			avg+=numpy.array(coiledCoils[i][key])
			i+=1
		avg/=len(coiledCoils)
		avgCoiledCoilParameters[key]=avg
	
	#calculate Crick angles for alpha-Helices
	avgCrickAngles=[]
	avgCrickAngles=calculateavgCrickAngles(helicesInSelection, avgCoiledCoilParameters, chains)
	
	#assign heptad positions
	heptadPositions=assignHeptadPositions(helicesInSelection, avgCrickAngles, stored.residueNumbersInSelection, selection)
	
	#print sequence with heptad positions
	printSequence(stored.residueNumbersInSelection, stored.residuesInSelection, heptadPositions)
	
	#print a nice table of the results
	table=makeResultsTable(helicesInSelection, stored.residuesInSelection, stored.residueNumbersInSelection, avgCoiledCoilParameters, avgAlphaHelixParameters, heptadPositions, avgCrickAngles)
	pprint_table(out, table)	
	
	#show Axes in pymol
	showAxesInPymol(helicesInSelection, avgCoiledCoilParameters)
	
cmd.extend( "pyTwister", pyTwister );