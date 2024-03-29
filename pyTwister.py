from pymol import cmd
from pymol import stored
import numpy
import math
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

def createPseudoatom (coordinates, objectName):
    x=float(coordinates[0])
    y=float(coordinates[1])
    z=float(coordinates[2])
    #print x,y,z,objectName
    posString="[%3.2f,%3.2f,%3.2f]" % (x,y,z)
    cmd.pseudoatom(pos=posString, object=objectName)

def calculateShortestDistance (P1,P2,P3,P4):
    """Calculate the line segment PaPb that is the shortest route between
        two lines P1P2 and P3P4. Calculate also the values of mua and mub where
        Pa = P1 + mua (P2 - P1)
        Pb = P3 + mub (P4 - P3)
        Return FALSE if no solution exists."""
    p21=P2-P1
    p43=P4-P3
    p13=P1-P3
    d1343=numpy.dot(p13,p43)
    d4321=numpy.dot(p43,p21)
    d1321=numpy.dot(p13,p21)
    d4343=numpy.dot(p43,p43)
    d2121=numpy.dot(p21,p21)
    denom = d2121 * d4343 - d4321 * d4321
    numer = d1343 * d4321 - d1321 * d4343
    mua = numer / denom
    mub = (d1343 + d4321 * mua) / d4343
    return [mua,mub]

def getCaPositions(selection):   
    #print selection
    stored.xyz = []
    cmd.iterate_state(1,selection,"stored.xyz.append([x,y,z])")
    caPositions=numpy.array(stored.xyz)
    return caPositions

def calculateBisections(caPositions):
    bisect=[]
    bisect.append([0,0,0])
    for i in range (1,len(caPositions)-1):
        ca0ca1=0
        normCa0Ca1=0
        ca1ca2=0
        normCa1Ca2=0
        ca0ca1=caPositions[i]-caPositions[i-1]
        normCa0Ca1=ca0ca1/numpy.linalg.norm(ca0ca1)
        #cmd.pseudoatom(pos=[normCa0Ca1[0],normCa0Ca1[1],normCa0Ca1[2]], object="normCa0Ca1")
        ca1ca2=caPositions[i+1]-caPositions[i]
        normCa1Ca2=ca1ca2/numpy.linalg.norm(ca1ca2)
        #cmd.pseudoatom(pos=[normCa1Ca2[0],normCa1Ca2[1],normCa1Ca2[2]], object="normCa1Ca2")
        bisect.append((normCa1Ca2-normCa0Ca1)+caPositions[i])
        #cmd.pseudoatom(pos=[biSectCa[0],biSectCa1[1],biSectCa1[2]], object="bisect")
    bisect=numpy.array(bisect)
    return bisect

def calculateHelixAxisPoints(caPositions, bisections):
    helixAxisPoints=[]
    #will be changed to extrapolated value later
    helixAxisPoints.append([0,0,0])
    for i in range (1,len(bisections)):
        #print i,len(bisections)
        if i < len(bisections)-1:
            #print i, len(caPositions), len(bisections)
            mua,mub=calculateShortestDistance (caPositions[i],bisections[i],caPositions[i+1],bisections[i+1])
            #print mua,mub
        elif i == len(bisections):
            #print i, len(caPositions), len(bisections)
            mua,mub=calculateShortestDistance (caPositions[i],bisections[i],caPositions[i-1],bisections[i-1])
            print mua, mub
        helixAxisPoints.append(caPositions[i]-mua*(caPositions[i]-bisections[i]))
    helixAxisPoints[0]=2*helixAxisPoints[1]-helixAxisPoints[2]
    #print helixAxisPoints
    helixAxisPoints.append(2*helixAxisPoints[len(helixAxisPoints)-1]-helixAxisPoints[len(helixAxisPoints)-2])
    return helixAxisPoints

def calculateAlphaHelicalRadius(caPositions, helixAxisPoints):
    alphaHelicalRadii=[]
    #alphaHelicalRadii.append(0)
    #alphaHelicalRadii.append(0)
    for i in range (0,len(helixAxisPoints)):
        radius=numpy.linalg.norm(caPositions[i]-helixAxisPoints[i])
        alphaHelicalRadii.append(radius)
    #print alphaHelicalRadii
    return alphaHelicalRadii

def calculateAlphaHelicalRisePerResidue(helixAxisPoints):
    alphaHelicalRisePerResidue=[]
    #alphaHelicalRisePerResidue.append(0)
    alphaHelicalRisePerResidue.append(0)
    for i in range (1,len(helixAxisPoints)-1):
        rise=(numpy.linalg.norm(helixAxisPoints[i-1]-helixAxisPoints[i])+numpy.linalg.norm(helixAxisPoints[i]-helixAxisPoints[i+1]))/2
        alphaHelicalRisePerResidue.append(rise)
    return alphaHelicalRisePerResidue

def calculateAlphaHelicalPhaseYieldPerResidue(caPositions, helixAxisPoints,correct):
    alphaHelicalPhaseYieldPerResidue=[]
    alphaHelicalPhaseYieldPerResidue.append(0)
    #alphaHelicalPhaseYieldPerResidue.append(0)
    for i in range (1,len(helixAxisPoints)-1):
        #generate pseudoatoms for dihedral calculation
        id=random.random()*1000
        cmd.pseudoatom(pos=str([caPositions[i-1].item(0),caPositions[i-1].item(1),caPositions[i-1].item(2)]), object="ca_i-1"+str(id))
        cmd.pseudoatom(pos=str([caPositions[i].item(0),caPositions[i].item(1),caPositions[i].item(2)]), object="ca_i"+str(id))
        cmd.pseudoatom(pos=str([caPositions[i+1].item(0),caPositions[i+1].item(1),caPositions[i+1].item(2)]), object="ca_i+1"+str(id))
        cmd.pseudoatom(pos=str([helixAxisPoints[i-1].item(0),helixAxisPoints[i-1].item(1),helixAxisPoints[i-1].item(2)]), object="axis_i-1"+str(id))
        cmd.pseudoatom(pos=str([helixAxisPoints[i].item(0),helixAxisPoints[i].item(1),helixAxisPoints[i].item(2)]), object="axis_i"+str(id))
        cmd.pseudoatom(pos=str([helixAxisPoints[i+1].item(0),helixAxisPoints[i+1].item(1),helixAxisPoints[i+1].item(2)]), object="axis_i+1"+str(id))
        a=cmd.get_dihedral("ca_i-1"+str(id),"axis_i-1"+str(id),"axis_i"+str(id),"ca_i"+str(id),state=0)
        b=cmd.get_dihedral("ca_i"+str(id),"axis_i"+str(id),"axis_i+1"+str(id),"ca_i+1"+str(id),state=0)
        #print i,a,b, (a+b)/2
        cmd.delete("ca_i-1"+str(id))
        cmd.delete("ca_i"+str(id))
        cmd.delete("ca_i+1"+str(id))
        cmd.delete("axis_i-1"+str(id))
        cmd.delete("axis_i"+str(id))
        cmd.delete("axis_i+1"+str(id))
        alphaHelicalPhaseYieldPerResidue.append(correct+(a+b)/2)
    return alphaHelicalPhaseYieldPerResidue

def calculateAlphaHelicalPitchPerResidue(helixAxisPoints,rise,phase):
    alphaHelicalPitchPerResidue=[]
    alphaHelicalPitchPerResidue.append(0)
    #alphaHelicalPitchPerResidue.append(0)
    alphaResiduesPerTurn=[]
    alphaResiduesPerTurn.append(0)
    #alphaResiduesPerTurn.append(0)
    for i in range (1,len(helixAxisPoints)-1):
        #print rise[i], phase[i]
        pitch=rise[i]*360/phase[i]
        residuesPerTurn=360/phase[i]
        #print residuesPerTurn
        alphaHelicalPitchPerResidue.append(pitch)
        alphaResiduesPerTurn.append(residuesPerTurn)
    return alphaResiduesPerTurn

def calculateCrickAngle(An, On, Cn):
    #print An, On, Cn
    OC=Cn-On
    normOC=OC/numpy.linalg.norm(OC)
    #print normOC
    OA=An-On
    normOA=OA/numpy.linalg.norm(OA)
    #print normOA
    dot=numpy.dot(normOC,normOA)
    cross=numpy.cross(normOA,normOC)
    #print cross
    angle=numpy.arccos(dot)*180/numpy.pi
    
    if (cross[2]<0):
        angle*=1
    #else:
    if (cross[2]>0):
        angle*=-1
    #print angle
    return angle

def pyTwister(selection, chains):
    '''
    DESCRIPTION
    Brief description what this function does goes here
    '''
    #some preliminary stuff
    #cmd.reinitialize()
    #cmd.load("Helix-coot-0.pdb")
    #cmd.delete("Alpha*")
    #cmd.delete("Axis*")
    #cmd.delete("coiledCoil*")
    cmd.set("dash_gap","0")
    cmd.set("dash_radius","0.2")
    #make list of residueNames
    stored.residueNumbersInSelection=[]
    cmd.iterate(selection + " & name ca", "stored.residueNumbersInSelection.append(resi)")
    stored.residuesInSelection=[]
    cmd.iterate(selection + " & name ca", "stored.residuesInSelection.append(resn)")
    stored.residueOneLetterCodes=[]
    
    helicesInSelection=[]
    
    #calculate alphaHelixParameters
    for i in range (0, len(chains)):
        dict={}
        #get ca positions
        dict['caPositions']=getCaPositions(selection + " & chain " + chains[i] + " & name ca")
        
        #reverse ca dictionary for the reverse calculations C-term to N-term will be averaged later
        dict['reverseCaPositions']=dict['caPositions'][::-1]
        
        #calculate bisections
        dict['bisections']=calculateBisections(dict['caPositions'])
        dict['reverseBisections']=calculateBisections(dict['reverseCaPositions'])
        
        #calculate helix axis points by adjusting the length of the bisection vectors
        dict['helixAxisPoints']=calculateHelixAxisPoints(dict['caPositions'], dict['bisections'])
        dict['reverseHelixAxisPoints']=calculateHelixAxisPoints(dict['reverseCaPositions'], dict['reverseBisections'])
        
        #some artistics to get the reversedReversed helixAxisPoints
        #dict['reverseHelixAxisPoints'].append([0,0,0]) #appends trailing 0,0,0 before turning
        #dict['reverseHelixAxisPoints'].remove([0,0,0]) #removes leading 0,0,0 before turning
        dict['reversedReverseHelixAxisPoints']=dict['reverseHelixAxisPoints'][::-1]
        
        #average helixAxisPoints from forward and reverse calculation
        dict['avgHelixAxisPoints']=(numpy.array(dict['helixAxisPoints'])+numpy.array(dict['reversedReverseHelixAxisPoints']))/2
        
        #calculate alphaHelix parameters with averaged helix axis
        dict['alphaHelicalRadiusPerResidue']=calculateAlphaHelicalRadius(dict['caPositions'], dict['avgHelixAxisPoints'])
        dict['alphaHelicalRisesPerResidue']=calculateAlphaHelicalRisePerResidue(dict['avgHelixAxisPoints'])
        dict['alphaHelicalPhaseYieldsPerResidue']=calculateAlphaHelicalPhaseYieldPerResidue(dict['caPositions'], dict['avgHelixAxisPoints'],0)
        dict['alphaHelicalResiduesPerTurn']=calculateAlphaHelicalPitchPerResidue(dict['avgHelixAxisPoints'],dict['alphaHelicalRisesPerResidue'],dict['alphaHelicalPhaseYieldsPerResidue'])
        helicesInSelection.append(dict)
   
    #calculate average of alpha helical parameters - average of AxisPoints is coiledCoilAxis
    averageAlphaHelicalParameters={}
    for key in ['avgHelixAxisPoints','alphaHelicalRadiusPerResidue', 'alphaHelicalRisesPerResidue', 'alphaHelicalPhaseYieldsPerResidue', 'alphaHelicalResiduesPerTurn']:
        average=0
        i=0
        while (i < len(helicesInSelection)):
            #print key
            average+=numpy.array(helicesInSelection[i][key])
            #print helicesInSelection[i][key], key
            i+=1
        average/=len(helicesInSelection)
        #print average
        averageAlphaHelicalParameters[key]=average
    
    #calculate coiled coil parameters
    coiledCoils=[]
    for i in range (0, len(chains)-1):
        coiledCoilParameters={}
        coiledCoilParameters['coiledCoilAxisPoints']=averageAlphaHelicalParameters['avgHelixAxisPoints']
        coiledCoilParameters['coiledCoilRadiusPerResidue']=calculateAlphaHelicalRadius(helicesInSelection[i]['avgHelixAxisPoints'], coiledCoilParameters['coiledCoilAxisPoints'])
        coiledCoilParameters['coiledCoilRisesPerResidue']=calculateAlphaHelicalRisePerResidue(coiledCoilParameters['coiledCoilAxisPoints'])
        coiledCoilParameters['coiledCoilPhaseYieldsPerResidue']=calculateAlphaHelicalPhaseYieldPerResidue(helicesInSelection[i]['avgHelixAxisPoints'], coiledCoilParameters['coiledCoilAxisPoints'],0)
        coiledCoilParameters['coiledCoilPitch']=numpy.array(calculateAlphaHelicalPitchPerResidue(coiledCoilParameters['coiledCoilAxisPoints'],coiledCoilParameters['coiledCoilRisesPerResidue'],numpy.absolute(numpy.array(coiledCoilParameters['coiledCoilPhaseYieldsPerResidue']))))*numpy.array(coiledCoilParameters['coiledCoilRisesPerResidue'])
        coiledCoils.append(coiledCoilParameters)
    
    #calculate average of coiled coil parameters
    averageCoiledCoilParameters={}
    for key in ['coiledCoilAxisPoints','coiledCoilRadiusPerResidue', 'coiledCoilRisesPerResidue', 'coiledCoilPhaseYieldsPerResidue', 'coiledCoilPitch']:
        average=0
        i=0
        while (i < len(coiledCoils)):
            average+=numpy.array(coiledCoils[i][key])
            i+=1
        average/=len(coiledCoils)
        averageCoiledCoilParameters[key]=average
    
    #calculate Crick angle and assign positions
    crickAnglesForAllHelices=[]
    i=0
    while (i < len(chains)):
        #print i
        crickAngles=[]
        crickAngles.append(0)
        #crickAngles.append(0)
        for j in range (1,len(helicesInSelection[i]["avgHelixAxisPoints"])-1):
            angle=calculateCrickAngle(helicesInSelection[i]['caPositions'][j], helicesInSelection[i]["avgHelixAxisPoints"][j], averageCoiledCoilParameters['coiledCoilAxisPoints'][j])
            crickAngles.append(angle)
        crickAnglesForAllHelices.append(crickAngles)
        i+=1
    avgCrickAngles=[]
    for i in range (0, len(crickAngles)):
        avg=0
        x=0
        y=0
        for j in range(0, len(chains)):
            x+=math.cos(crickAnglesForAllHelices[j][i]*numpy.pi/180)
            y+=math.sin(crickAnglesForAllHelices[j][i]*numpy.pi/180)
        avg=math.atan2(y,x)*180/numpy.pi
        avgCrickAngles.append(avg)
    avgCrickAngles.append(0)#print avgCrickAngles
    pos=[]
    pos.append("x")
    #pos.append("x")
    i = 1
    while (i < len(helicesInSelection[0]["avgHelixAxisPoints"])):
        #print avgCrickAngles[i]
        addUp=0
        if (avgCrickAngles[i-1] < 0 and avgCrickAngles[i] > 0 and math.fabs(avgCrickAngles[i-1]) > math.fabs(avgCrickAngles[i]) and math.fabs(avgCrickAngles[i]) < 50):
            #print avgCrickAngles[i]
            pos.append("a")
            cmd.color("red", selection+" & resi "+str(stored.residueNumbersInSelection[i]))
            if i+1<len(helicesInSelection[0]["avgHelixAxisPoints"])-1:
                pos.append("b")
                addUp+=1
            if i+2<len(helicesInSelection[0]["avgHelixAxisPoints"])-1:
                pos.append("c")
                addUp+=1
        elif (avgCrickAngles[i] < 0 and avgCrickAngles[i+1] > 0 and math.fabs(avgCrickAngles[i]) < math.fabs(avgCrickAngles[i+1]) and math.fabs(avgCrickAngles[i]) < 50):
            pos.append("d")
            cmd.color("blue", selection+" & resi "+str(stored.residueNumbersInSelection[i]))
            if i+1<len(helicesInSelection[0]["avgHelixAxisPoints"])-1:
                pos.append("e")
                addUp+=1
            if i+2<len(helicesInSelection[0]["avgHelixAxisPoints"])-1:
                pos.append("f")
                addUp+=1
            if i+3<len(helicesInSelection[0]["avgHelixAxisPoints"])-1:
                pos.append("g")
                addUp+=1
        else:
            pos.append("?")
        #print addUp
        i+=addUp
        i+=1
    
    #make a nice table of results
    table = [["Res.", "cc-rad", "cc-rise", "cc-pitch", "cc-phaseYield", "pos","Crick-angle", "a-radius", "a-rise", "a-res/turn", "a-phaseYield"]]
    for i in range (1,len(helicesInSelection[0]["avgHelixAxisPoints"])-1):
        #print i
        list=[]
        list.append(stored.residuesInSelection[i]+ " " +stored.residueNumbersInSelection[i])
        list.append(averageCoiledCoilParameters['coiledCoilRadiusPerResidue'][i])
        list.append(averageCoiledCoilParameters['coiledCoilRisesPerResidue'][i])
        list.append(averageCoiledCoilParameters['coiledCoilPitch'][i])
        list.append(averageCoiledCoilParameters['coiledCoilPhaseYieldsPerResidue'][i])
        list.append(pos[i])
        list.append(avgCrickAngles[i])
        list.append(averageAlphaHelicalParameters['alphaHelicalRadiusPerResidue'][i])
        list.append(averageAlphaHelicalParameters['alphaHelicalRisesPerResidue'][i])
        list.append(averageAlphaHelicalParameters['alphaHelicalResiduesPerTurn'][i])
        list.append(averageAlphaHelicalParameters['alphaHelicalPhaseYieldsPerResidue'][i])
        table.append(list)
    averages = []
    stdevs = []
    avgCoiledCoilRadius = numpy.mean(averageCoiledCoilParameters['coiledCoilRadiusPerResidue'])
    stdevCoiledCoilRadius = numpy.std(averageCoiledCoilParameters['coiledCoilRadiusPerResidue'])
    
    avgCoiledCoilRise = numpy.mean(averageCoiledCoilParameters['coiledCoilRisesPerResidue'])
    stdevCoiledCoilRise = numpy.std(averageCoiledCoilParameters['coiledCoilRisesPerResidue'])
    
    avgCoiledCoilPitch = numpy.mean(averageCoiledCoilParameters['coiledCoilPitch'])
    stdevCoiledCoilPitch = numpy.std(averageCoiledCoilParameters['coiledCoilPitch'])
    
    avgCoiledCoilPhaseYield = numpy.mean(averageCoiledCoilParameters['coiledCoilPhaseYieldsPerResidue'])
    stdevCoiledCoilPhaseYield = numpy.std(averageCoiledCoilParameters['coiledCoilPhaseYieldsPerResidue'])
    
    avgAlphaHelicalRadius = numpy.mean(averageAlphaHelicalParameters['alphaHelicalRadiusPerResidue'])
    stdevAlphaHelicalRadius = numpy.std(averageAlphaHelicalParameters['alphaHelicalRadiusPerResidue'])
    
    avgAlphaHelicalRise = numpy.mean(averageAlphaHelicalParameters['alphaHelicalRisesPerResidue'])
    stdevAlphaHelicalRise = numpy.std(averageAlphaHelicalParameters['alphaHelicalRisesPerResidue'])
    
    avgResPerTurn = numpy.mean(averageAlphaHelicalParameters['alphaHelicalResiduesPerTurn'])
    stdevResPerTurn = numpy.std(averageAlphaHelicalParameters['alphaHelicalResiduesPerTurn'])
    
    avgAlphaHelicalPhaseYield = numpy.mean(averageAlphaHelicalParameters['alphaHelicalPhaseYieldsPerResidue'])
    stdevAlphaHelicalPhaseYield = numpy.std(averageAlphaHelicalParameters['alphaHelicalPhaseYieldsPerResidue'])
    
    averages = ["Avg", avgCoiledCoilRadius, avgCoiledCoilRise, avgCoiledCoilPitch, avgCoiledCoilPhaseYield, "-","-", avgAlphaHelicalRadius, avgAlphaHelicalRise, avgResPerTurn, avgAlphaHelicalPhaseYield]
    stdevs = ["Std", stdevCoiledCoilRadius, stdevCoiledCoilRise, stdevCoiledCoilPitch, stdevCoiledCoilPhaseYield, "-","-", stdevAlphaHelicalRadius, stdevAlphaHelicalRise, stdevResPerTurn, stdevAlphaHelicalPhaseYield]
    empty = ["","","","","","","","","","",""]
    
    table.append(empty)
    table.append(averages)
    table.append(stdevs)
    #print table
    pprint_table(out, table)
    
    #show axes in pymol
    for i in range (0, len(helicesInSelection)):
        #print len(helicesInSelection),i
        
        for j in range (0,len(helicesInSelection[i]["avgHelixAxisPoints"])):
            createPseudoatom(coiledCoilParameters["coiledCoilAxisPoints"][j],"coiledCoilAxisPoint_"+str(i)+"_"+str(j))
            createPseudoatom(helicesInSelection[i]["avgHelixAxisPoints"][j],"AxisPoint_"+str(i)+"_"+str(j))
            if j > 0:
                cmd.distance("AxisLine_"+str(i)+"_"+str(j),"AxisPoint_"+str(i)+"_"+str(j-1),"AxisPoint_"+str(i)+"_"+str(j))
                cmd.distance("CoiledCoilAxisLine_"+str(i)+"_"+str(j),"coiledCoilAxisPoint_"+str(i)+"_"+str(j-1),"coiledCoilAxisPoint_"+str(i)+"_"+str(j))
        cmd.group("AlphaAxis_"+str(i), "AxisLine_"+str(i)+"*")
        cmd.group("AlphaAxisPoints_"+str(i), "AxisPoint_"+str(i)+"*")
        cmd.group("coiledCoil", "coiledCoil*")
        cmd.set("dash_color","yellow", "Alpha*")
        cmd.set("dash_color","red", "Coiled*")
        #print i
    cmd.hide("labels")
    cmd.hide("nonbonded")
    
cmd.extend( "pyTwister", pyTwister );



