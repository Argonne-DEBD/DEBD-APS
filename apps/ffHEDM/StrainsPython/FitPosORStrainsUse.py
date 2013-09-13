# Code to do fitting of position, orientation and strains in a sample.
#
# Author: Hemant Sharma
# Version: 1
# Dated: 2013/05/28
# - Initial test version of the code to see if it is viable to implement in python.
#
# Version: 2
# Dated: 2013/06/13
# - Everything is working and fitting using minuit is now implemented.
#
# Version: 3
# Dated: 2013/06/21
# - Code ready for deployment, multiprocessing is not done yet.
#
# Version: 4
# - Multiprocessing is implemented using SWIFT now.

from __future__ import division
from __future__ import with_statement
import sys, getopt, time, numpy, math
numpy.set_printoptions(threshold=numpy.nan)
from math import pi
import math
import os
import csv
from numpy import genfromtxt
from numpy import zeros
from numpy import ones
from numpy import reshape
from numpy import shape
from numpy import array
from numpy import arccos
from numpy import nonzero
from numpy import arctan2 as atan2
from numpy import transpose
from numpy import append
from numpy import absolute as abs
from numpy import linalg as LA
from GenHKLs2 import GenerateHKLs2
from CalcAngleErrors import CalcAngleErrors
from OrientMat2Euler import OrientMat2Euler
from Euler2OrientMat import Euler2OrientMat
from FitErrorsPosT import FitErrorsPosT
from CheckForGoodSpots import CheckForGoodSpots
from FitErrorsOrientStrains import FitErrorsOrientStrains
from FitErrorsStrains import FitErrorsStrains
from CalcStrainTensor import CalcStrainTensor
from DoFitMinuit import DoFitMinuit

def degrees(x):
    return x*180.0/pi

def radians(x):
    return x*pi/180.0

def BestPosFileRead(BestPosFileName):
    statinfo = os.stat(BestPosFileName)
    filesize = int(statinfo.st_size)
    if filesize > 0:
        array1 = genfromtxt(BestPosFileName,filling_values=0,usecols=range(0,14),skip_header=4,delimiter=',',dtype=float)
        return array1
    else:
        print 'The file did not exist or was empty, skipping.'
        array1 = array([None])
        return array1

def FindTopInfoBestPosFile(FileName):
    try:
        with open(FileName,'rb') as f:
            array1 = []
            reader = csv.reader(f)
            try:
                for row in reader:
                    array1.append(row)
                array1 = array1[0:3]
                return array1
            except csv.Error, e:
                array1 = 1
    except IOError:
        array1 = 1

def AllInfoRead(AllInfoFile):
    AllSpots = genfromtxt(AllInfoFile, skip_header=1, dtype=float)
    AllSpots = AllSpots[AllSpots[:,2].argsort(),]
    return AllSpots

def SpotIDsRead(filename):
    info = genfromtxt(filename,dtype=int)
    return info

def ParamsTestFileRead(ParamsTestFileName):
    file = open( ParamsTestFileName, "r" )
    InputParaMeters = []
    for line in file:
        InputParaMeters.append( line )
    file.close()
    return InputParaMeters


def DoFitPosORStrains(argv):
    starttime = time.time()
    ParamsTestFileName = ''
    SpotsfileName = ''
    InputFile = ''
    try:
        opts, args = getopt.getopt(argv,"hp:s:i:",["paramsfile=","spotsfile=","inputfile="])
    except getopt.GetoptError:
        print '\n\nUsage: FitPosORStrainsPython.py -p paramstest.txt -s SpotID -i InputAllExtraInfoFittingAll.csv\n\n'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '\n\nUsage: FitPosORStrainsPython.py -p paramstest.txt -s SpotID -i InputAllExtraInfoFittingAll.csv\n\n'
            sys.exit()
        elif opt in ("-p", "--paramsfile"):
            ParamsTestFileName = arg
        elif opt in ("-s", "--spotsfile"):
            SpId = int(arg)
        elif opt in ("-i", "--inputfile"):
            InputFile = arg

    # Read paramstest.txt file and generate initialization parameters
    InputParaMeters = ParamsTestFileRead(ParamsTestFileName)
    NumberParameters = len(InputParaMeters)
    RingSizes = []
    nRings = []
    OmegaRange = []
    BoxSize = []
    for eachline in InputParaMeters:
        splitresult1 = eachline.split(';')
        splitresult2 = splitresult1[0]
        splitresult = splitresult2.split( )
        if 'LatticeConstant' in splitresult:
            LatticeConstant = float(splitresult[1])
        if 'CellStruct' in splitresult:
            Str = float(splitresult[1])
        if 'Wavelength' in splitresult:
            Wavelength = float(splitresult[1])
        if 'Distance' in splitresult:
            Lsd = float(splitresult[1])
        if 'MinMatchesToAcceptFrac' in splitresult:
            MinFracToAccept = float(splitresult[1]) - 0.15
        if 'MarginOme' in splitresult:
            tol = float(splitresult[1])
        if 'MarginRadial' in splitresult:
            marginPos = float(splitresult[1])
        if 'Rsample' in splitresult:
            Rsample = float(splitresult[1])
        if 'Hbeam' in splitresult:
            Hbeam = float(splitresult[1])
        if 'ExcludePoleAngle' in splitresult:
            MinEta = float(splitresult[1])
        if 'Wedge' in splitresult:
            wedge = float(splitresult[1])
        if 'MarginRadial' in splitresult:
            marginPos = float(splitresult[1])
        if 'RingRadii' in splitresult:
            RingSizes.append(float(splitresult[1]))
        if 'RingNumbers' in splitresult:
            nRings.append(float(splitresult[1]))
        if 'OmegaRange' in splitresult:
            OmegaRange.append(float(splitresult[1]))
            OmegaRange.append(float(splitresult[2]))
        if 'BoxSize' in splitresult:
            BoxSize.append(float(splitresult[1]))
            BoxSize.append(float(splitresult[2]))
            BoxSize.append(float(splitresult[3]))
            BoxSize.append(float(splitresult[4]))
        if 'OutputFolder' in splitresult:
            OutFolder = splitresult[1]
    nrrows = (len(BoxSize))/4
    BoxSize = numpy.reshape(BoxSize, (nrrows,-1))
    OmegaRange = numpy.reshape(OmegaRange, (nrrows,-1))
    a = LatticeConstant
    b = LatticeConstant
    c = LatticeConstant
    alph = 90
    bet = 90
    gamm = 90

    OutputFolder =  'Output/'
    ResultFolder = 'Results/'

    # Some other initializations
    MargOme         = 4
    MargPos         = Rsample
    MargPos2        = Rsample/2
    MargOme2        = 4
    MargABC         = 0.5 #in %
    MargABG         = 0.5 #in %
    chi             = 0
    LatCin          = numpy.array([LatticeConstant, LatticeConstant, LatticeConstant, 90, 90, 90])
    NumberSkipped   = 0 # nrtr in the old code
    hkls = GenerateHKLs2(nRings,Str,RingSizes,LatticeConstant,Wavelength)

    # Read Info and initialize matrices
    AllInfoFile = InputFile
    AllSpots = AllInfoRead(AllInfoFile)
    AllSpotsYZO = AllSpots[:,[0,1,2,4,8,9,10,5]]
    nrSpIds = 1
    OutFN = ResultFolder + 'OrientPosFit' + str(SpId) + '.csv'
    OrigOutFN = ResultFolder + 'OrientPosIndexer' + str(SpId) + '.csv'
    SpIdDiffFN = ResultFolder + 'NrDiffIndexerFit' + str(SpId) + '.csv'
    OutFNf = open(OutFN,'w')
    OrigOutFNf = open(OrigOutFN,'w')
    SpIdDiffFNf = open(SpIdDiffFN,'w')
    OrientsOrig = zeros((nrSpIds,10))
    PositionsOrig = zeros((nrSpIds,4))
    ErrorsOrig = zeros((nrSpIds,4))
    NrDiff = zeros((nrSpIds,4))
    OrientsFit = zeros((nrSpIds,10))
    PositionsFit = zeros((nrSpIds,4))
    StrainsFit = zeros((nrSpIds,13))
    ErrorsFin = zeros((nrSpIds,4))

    header = 'SpotID,YObsCorrPos,ZObsCorrPos,OmegaObsCorrPos,G1Obs,G2Obs,G3Obs,YExp,ZExp,OmegaExp,G1Exp,G2Exp,G3Exp,YObsCorrWedge,ZObsCorrWedge,OmegaObsCorrWedge,OmegaObs,YObs,ZObs,InternalAngle,DiffLen,DiffOmega\n'

    nSpID = 0
    timestart = time.time()
    print 'Spot ID being processed: %d, %d out of %d' % ( SpId, nSpID+1, nrSpIds)
    FileNumber = '%09d' % (SpId)
    FileName = OutputFolder + 'BestPos_' + FileNumber + '.csv'
    SpotsCompFN = OutputFolder + 'FitBest_' + FileNumber + '.csv'
    best = BestPosFileRead(FileName)
    if all(v is None for v in best):
        f = open(SpotsCompFN,'w')
        f.write(header)
        f.close()
        OutFNf.close()
        OrigOutFNf.close()
        SpIdDiffFNf.close()
        print 'This spot did not have a bestpos file.'
        sys.exit()
    TopInfo = FindTopInfoBestPosFile(FileName)
    OrientPos = TopInfo[2]
    OrientPos0 = array([float(x) for x in OrientPos])
    Orient0 = OrientPos0[1:10]
    Pos0 = OrientPos0[10:13]
    Euler0 = OrientMat2Euler(Orient0)
    IA0 = OrientPos0[0]
    rows = nonzero(best[:,2])
    Spots = best[nonzero(best[:,2])[0],:]
    spotIDS = Spots[:,13]
    nrSpotsMatched = len(spotIDS)
    spotsYZO = zeros((nrSpotsMatched,8))
    for nrSpotsDone in range(0,nrSpotsMatched):
        ThisID = int(spotIDS[nrSpotsDone])
        RowNr = nonzero(AllSpotsYZO[:,3] == ThisID)[0][0]
        ThisInfo = AllSpotsYZO[RowNr][0:8]
        spotsYZO[nrSpotsDone][0:8] = ThisInfo
    Ini = zeros((12))
    Ini[0:3] = Pos0
    Ini[3:6] = Euler0
    Ini[6:] = LatCin
    Result = CalcAngleErrors(Ini,spotsYZO,hkls,Lsd,Wavelength,OmegaRange,BoxSize,MinEta,wedge,chi)
    ErrorIni = Result['Error']
    Splist = Result['SPList']
    print 'Initial error:  %f %f %f' % (ErrorIni[0],ErrorIni[1],ErrorIni[2])
    spotsYZO = Splist
    OrientsOrig[nSpID,0] = SpId
    OrientsOrig[nSpID,1:10] = Orient0
    PositionsOrig[nSpID,0] = SpId
    PositionsOrig[nSpID,1:4] = Pos0
    ErrorsOrig[nSpID,0] = SpId
    ErrorsOrig[nSpID,1:4] = ErrorIni
    Inp = Ini
    X = Pos0
    X0 = Inp
    XLow = X - MargPos2
    XHigh = X + MargPos2
    EulerLow = Euler0 - MargOme2
    EulerHigh = Euler0 + MargOme2
    if XLow[0] < -Rsample:
        XLow[0] = -Rsample
    if XLow[1] < -Rsample:
        XLow[1] = -Rsample
    if XLow[2] < -Hbeam/2:
        XLow[2] = -Hbeam/2
    if XHigh[0] > Rsample:
        XHigh[0] = Rsample
    if XHigh[1] > Rsample:
        XHigh[1] = Rsample
    if XHigh[2] > Hbeam/2:
        XHigh[2] = Hbeam/2
    lb = zeros((12))
    ub = zeros((12))
    lb[0:3] = XLow
    lb[3:6] = EulerLow
    lb[6:12] = array([(a*(1-(MargABC/100))), (b*(1-(MargABC/100))), (c*(1-(MargABC/100))), (alph*(1-(MargABG/100))), (bet*(1-(MargABG/100))), (gamm*(1-(MargABG/100)))])
    ub[0:3] = XHigh
    ub[3:6] = EulerHigh
    ub[6:12] = array([(a*(1+(MargABC/100))), (b*(1+(MargABC/100))), (c*(1+(MargABC/100))), (alph*(1+(MargABG/100))), (bet*(1+(MargABG/100))), (gamm*(1+(MargABG/100)))])
    # Optimize position
    FitType = 'Pos'
    Result = DoFitMinuit(X0,spotsYZO,hkls,Lsd,Wavelength,OmegaRange,BoxSize,MinEta,wedge,chi,Euler0,Pos0,LatCin,FitType,lb,ub,SpId)
    X = Result
    Inp = X
    X = Inp[0:3]
    PosFit = X
    aF = LatCin[0]
    bF = LatCin[1]
    cF = LatCin[2]
    alphF = LatCin[3]
    betF = LatCin[4]
    gammF = LatCin[5]
    # Check if spots are good or not
    InpTr = Inp
    Inp = Ini
    Inp[0:3] = X
    Result = CheckForGoodSpots(Inp,AllSpotsYZO,hkls,Lsd,Wavelength,OmegaRange,BoxSize,MinEta,tol,marginPos,wedge,chi)
    SpotsYZORedef = Result['SpotsYZORedef']
    nrExp = Result['nrExp']
    nrObs = Result['nrObs']
    SpotsYZORedef = SpotsYZORedef[SpotsYZORedef[:,0]!=0,:]
    Result2 = CalcAngleErrors(Inp,SpotsYZORedef,hkls,Lsd,Wavelength,OmegaRange,BoxSize,MinEta,wedge,chi)
    ErrorInt = Result2['Error']
    SpotsYZORe = Result2['SPList']
    print 'Interim error:  %f %f %f' % (ErrorInt[0],ErrorInt[1],ErrorInt[2])
    SpotsYZORedef = SpotsYZORe
    if (nrObs/nrExp) < MinFracToAccept:
        print 'Not good enough number of spots found, skipping this spot.'
        NumberSkipped += 1
        f = open(SpotsCompFN,'w')
        f.write(header)
        f.close()
        OutFNf.close()
        OrigOutFNf.close()
        SpIdDiffFNf.close()
        sys.exit()
    # Now optimize orientation
    X0 = Euler0[0]
    X0 = numpy.append(X0,LatCin,axis=1)
    EulerLow = Euler0 - MargOme2
    EulerHigh = Euler0 + MargOme2
    lb = zeros((9))
    ub = zeros((9))
    lb[0:3] = EulerLow
    lb[3:9] = array([(aF*(1-(MargABC/100))), (bF*(1-(MargABC/100))), (cF*(1-(MargABC/100))), (alphF*(1-(MargABG/100))), (betF*(1-(MargABG/100))), (gammF*(1-(MargABG/100)))])
    ub[0:3] = EulerHigh
    ub[3:9] = array([(aF*(1+(MargABC/100))), (bF*(1+(MargABC/100))), (cF*(1+(MargABC/100))), (alphF*(1+(MargABG/100))), (betF*(1+(MargABG/100))) ,(gammF*(1+(MargABG/100)))])
    FitType = 'Orient'
    Result = DoFitMinuit(X0,SpotsYZORedef,hkls,Lsd,Wavelength,OmegaRange,BoxSize,MinEta,wedge,chi,Euler0,PosFit,LatCin,FitType,lb,ub,SpId)
    X = Result
    Inp = X
    X0 = Inp
    EulerFit =  Inp[0:3]
    LatCFit = Inp[3:9]
    OMFit = Euler2OrientMat(EulerFit)
    ResultStrain = CalcStrainTensor(LatCin,LatCFit,OMFit)
    StrainTensorGr = ResultStrain['StrainTensorGr']
    StrainTensorSample = ResultStrain['StrainTensorSample']
    aF = LatCin[0]
    bF = LatCin[1]
    cF = LatCin[2]
    alphF = LatCin[3]
    betF = LatCin[4]
    gammF = LatCin[5]
    InpTr = PosFit
    InpTr = numpy.append(InpTr,EulerFit,axis=1)
    InpTr = numpy.append(InpTr,LatCFit,axis=1)
    # Check if the spots are good or not
    ResultCheck = CheckForGoodSpots(InpTr,AllSpotsYZO,hkls,Lsd,Wavelength,OmegaRange,BoxSize,MinEta,tol,marginPos,wedge,chi)
    SpotsYZORedef = ResultCheck['SpotsYZORedef']
    nrExp = ResultCheck['nrExp']
    nrObs = ResultCheck['nrObs']
    SpotsYZORedef = SpotsYZORedef[SpotsYZORedef[:,0]!=0,:]
    Result2 = CalcAngleErrors(InpTr,SpotsYZORedef,hkls,Lsd,Wavelength,OmegaRange,BoxSize,MinEta,wedge,chi)
    ErrorInt = Result2['Error']
    SpotsYZORe = Result2['SPList']
    print 'Interim error:  %f %f %f' % (ErrorInt[0],ErrorInt[1],ErrorInt[2])
    SpotsYZORedef = SpotsYZORe
    if (nrObs/nrExp) < MinFracToAccept:
        print 'Not good enough number of spots found, skipping this spot.'
        NumberSkipped += 1
        f = open(SpotsCompFN,'w')
        f.write(header)
        f.close()
        OutFNf.close()
        OrigOutFNf.close()
        SpIdDiffFNf.close()
        sys.exit()
    # Now optimize strains
    X0 = LatCFit
    lb = zeros((6))
    ub = zeros((6))
    lb[0:6] = array([(aF*(1-(MargABC/100))), (bF*(1-(MargABC/100))), (cF*(1-(MargABC/100))), (alphF*(1-(MargABG/100))), (betF*(1-(MargABG/100))), (gammF*(1-(MargABG/100)))])
    ub[0:6] = array([(aF*(1+(MargABC/100))), (bF*(1+(MargABC/100))), (cF*(1+(MargABC/100))), (alphF*(1+(MargABG/100))), (betF*(1+(MargABG/100))) ,(gammF*(1+(MargABG/100)))])
    FitType = 'Strain'
    ResultStrains = DoFitMinuit(X0,SpotsYZORedef,hkls,Lsd,Wavelength,OmegaRange,BoxSize,MinEta,wedge,chi,EulerFit,PosFit,LatCin,FitType,lb,ub,SpId)
    X = ResultStrains
    OrientationMatrixFit = Euler2OrientMat(EulerFit)
    LatticeParameterFit = X
    PositionFit = PosFit
    OrientFit = reshape(OrientationMatrixFit,9)
    XTr = numpy.append(numpy.append(PosFit,EulerFit,axis=1),LatticeParameterFit,axis=1)
    # Check if the spots are good or not
    ResultCheck = CheckForGoodSpots(XTr,AllSpotsYZO,hkls,Lsd,Wavelength,OmegaRange,BoxSize,MinEta,tol,marginPos,wedge,chi)
    SpotsYZORedef = ResultCheck['SpotsYZORedef']
    nrExp = ResultCheck['nrExp']
    nrObs = ResultCheck['nrObs']
    SpotsYZORedef = SpotsYZORedef[SpotsYZORedef[:,0]!=0,:]
    ResultsAngleErrorsFin = CalcAngleErrors(InpTr,SpotsYZORedef,hkls,Lsd,Wavelength,OmegaRange,BoxSize,MinEta,wedge,chi)
    AngleError = ResultsAngleErrorsFin['Error']
    print 'Interim error:  %f %f %f' % (AngleError[0],AngleError[1],AngleError[2])
    SpotsCompared = ResultsAngleErrorsFin['SpotsComp']
    SpotsYZORedef = ResultsAngleErrorsFin['SPList']
    ResultsStrainTensorsFin = CalcStrainTensor(LatCin,LatticeParameterFit,OrientationMatrixFit)
    StrainTensorGr = ResultsStrainTensorsFin['StrainTensorGr']
    StrainTensorSample = ResultsStrainTensorsFin['StrainTensorSample']
    # Now optimize position again
    X0 = PosFit
    XLow = X0 - MargPos2
    XHigh = X0 + MargPos2
    if XLow[0] < -Rsample:
        XLow[0] = -Rsample
    if XLow[1] < -Rsample:
        XLow[1] = -Rsample
    if XLow[2] < -Hbeam/2:
        XLow[2] = -Hbeam/2
    if XHigh[0] > Rsample:
        XHigh[0] = Rsample
    if XHigh[1] > Rsample:
        XHigh[1] = Rsample
    if XHigh[2] > Hbeam/2:
        XHigh[2] = Hbeam/2
    lb = XLow
    ub = XHigh
    FitType = 'P2os'
    ResultPosFit = DoFitMinuit(X0,SpotsYZORedef,hkls,Lsd,Wavelength,OmegaRange,BoxSize,MinEta,wedge,chi,EulerFit,PosFit,LatticeParameterFit,FitType,lb,ub,SpId)
    PositionFit = ResultPosFit
    PosFit = PositionFit
    XTr = numpy.append(numpy.append(PosFit,EulerFit,axis=1),LatticeParameterFit,axis=1)
    ResultsAngleErrorsFin = CalcAngleErrors(XTr,SpotsYZORedef,hkls,Lsd,Wavelength,OmegaRange,BoxSize,MinEta,wedge,chi)
    AngleError = ResultsAngleErrorsFin['Error']
    SpotsCompared = ResultsAngleErrorsFin['SpotsComp']
    SpotsYZORedef = ResultsAngleErrorsFin['SPList']
    print 'Final error:    %f %f %f' % (AngleError[0],AngleError[1],AngleError[2])
    print 'Strain Tensors:'
    print 'Grain:'
    print StrainTensorGr
    print 'Sample:'
    print StrainTensorSample
    print 'Fitted lattice parameter is:'
    print LatticeParameterFit
    print 'Original orientation was:'
    print reshape(Orient0,[3,3])
    print 'Orient position was:'
    print Pos0
    print 'Average error in distance was %f.' % (ErrorIni[0])
    print 'Average error in omega was %f' % (ErrorIni[1])
    print 'Average internal angle was %f' % (ErrorIni[2])
    print 'Fitted orientation is:'
    print OrientationMatrixFit
    print 'Fitted position is:'
    print PositionFit
    print 'Average error in distance is %f.' % (AngleError[0])
    print 'Average error in omega is %f' % (AngleError[1])
    print 'Average internal angle is %f' % (AngleError[2])
    DiffPos = PositionFit - Pos0
    LenDiffPos = LA.norm(DiffPos)
    print 'Difference in position before and after fitting is %f microns.' % (LenDiffPos)
    OrientsFit[nSpID,0] = SpId
    PositionsFit[nSpID,0] = SpId
    ErrorsFin[nSpID,0] = SpId
    StrainsFit[nSpID,0] = SpId
    OrientsFit[nSpID,1:10] = OrientFit
    PositionsFit[nSpID,1:4] = PositionFit
    StrainsFit[nSpID,1:7] = append(append(StrainTensorGr[0,0:3],StrainTensorGr[1,1:3]),StrainTensorGr[2,2])
    StrainsFit[nSpID,7:13] = append(append(StrainTensorSample[0,0:3],StrainTensorSample[1,1:3]),StrainTensorSample[2,2])
    SpIDsIndexer = best[:,13]
    SpIDsFit = SpotsCompared[:,0]
    AreSame = numpy.intersect1d(SpIDsIndexer,SpIDsFit)
    NrDiff[nSpID,0:4] = array([SpId,nrExp,nrObs, (len(SpIDsFit)-len(AreSame))])
    ErrorsFin[nSpID,1:4] = AngleError
    f = open(SpotsCompFN,'w')
    f.write(header)
    f.close()
    f = open(SpotsCompFN,'a')
    numpy.savetxt(f,SpotsCompared,delimiter=',')
    f.close()
    nSpID +=1
    timeend = time.time()
    print 'Time taken for this spot: ' + str(timeend - timestart) + ' s.'
    OutRow1 = numpy.append(numpy.append(OrientsOrig,PositionsOrig,axis=1),ErrorsOrig,axis=1)
    OutRow2 = numpy.append(numpy.append(numpy.append(OrientsFit,PositionsFit,axis=1),StrainsFit,axis=1),ErrorsFin,axis=1)
    numpy.savetxt(SpIdDiffFNf,NrDiff,delimiter=',')
    numpy.savetxt(OrigOutFNf,OutRow1,delimiter=',')
    numpy.savetxt(OutFNf,OutRow2,delimiter=',')
    OutFNf.close()
    OrigOutFNf.close()
    SpIdDiffFNf.close()
    finishtime = time.time()
    print 'Total time taken for this run: ' + str(finishtime - starttime) + ' s.'

if __name__ == "__main__":
    import sys
    DoFitPosORStrains(sys.argv[1:])
