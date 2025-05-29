
#! /usr/bin/env python
""" Markov-type simulation model of ecological networks
    The size of networks are variable as are the interaction
    strengths and their number. The transition matrix is as follows:

    columns are donors [row,col]
    rows are receptors [pred,prey]
    
    transition flows are donor biomass specific flows

    """
    
# this simulation depends upon the following freely available programs:
#
# Numeric Python (v. 21.0)
#   http://www.pfdubois.com/numpy/ -- numerical computations/array handling
# Gnuplot.py (v. 1.6)
#   http://gnuplot-py.sourceforge.net/ -- a wrapper to gnuplot 
# gnuplot (v. 3.71)
#   http://www.gnuplot.org/ -- a multi-platform scientific plotting program

from numarray import *      # load Numeric 
from Numeric import *      # load Numeric 
from RandomArray import *  # load random number generators (part of Numeric)
import Gnuplot, Gnuplot.funcutils   # load Gnuplot.py


def wVector(nMean, wVectorBase):
    """ from input of nMean, create a matrix with a given random method
    Currently this is the Poisson distribution
    determine biomass size classes (g wet weight, midpoints on a log10 scale)
    declare number of compartments within each size compartment
    and then populate them
    """
    
    nInEachW = poisson(nMean, shape(wVectorBase))
    wVector = repeat(wVectorBase, nInEachW)
    wNew = size(wVector)+2
    wVector = resize(wVector, [wNew,])
    return wVector

def rSpecific(weight, rSlope, rInt):
    """ obtain specific respiration rates to be used for calculation of total R
    respiration rates are from Peters 1981 (basal metabolic rates in Watts)
    """
    
    rIntercept = where(weight <= 10**(-7.), rInt[0], rInt[1])
    rSpecific = rIntercept * ((weight/1000)**rSlope)
    # division by 1000 converts to kg
    
    """ assume: 1 W = 20.65 kcal day-1 = 7537.2 kcal year-1
                      {as 1 g wet = 1.673 kcal; see Peters 1983}
                    = 7537.2 kcal year-1 / (1.673 kcal g-1)
                    = 4505.2 g wet biomass /year
    """
    
    rSpecific = rSpecific * 4505.2  # convert from Watts to g wet year-1 ..C?
    return rSpecific

def hMaskGet(proportion, matrixSize):
    """ classifiy compartments as being autotrophic or heterotrophic groups
    no overlapping groups are implemented, binomial distribution is assumed
    """

    hMask = binomial(1, proportion, matrixSize)
    return hMask

def dMaskGet(proportion, matrixSize, hMask):
    """ classifiy compartments as being detrtivores amongst heterotrophic
    groups binomial distribution is assumed
    """
    for consumer in range(nBox):
        if hMask[consumer]: dMask = binomial(1, proportion, matrixSize)
    
    return dMask

def tMatrixMake(tMatrix, hMask, nDonorsMean):
    """ assign GENERAL AUTOTROPH & HETEROTROPH FLOWS
    
    assign magnitude of donor biomass specific flows -- assumed to be 
    a uniform random number with constraint that the sum of all internal 
    specific outflows < 1
    """

    nBox = tMatrix.shape[1]
    nConsumersPerDonor = poisson(nDonorsMean, nBox)
    for donor in range(nBox): 
        nConsumers = nConsumersPerDonor[donor]
        while nConsumers > 0:
            while 1:
                consumer = randint(0,nBox)
                if ((consumer != donor) and hMask[consumer]):
                    tMatrix[consumer,donor] = random() # check this!
                    nConsumers = nConsumers-1
                    break
    sums = sum(ceil(tMatrix, axis=0))
    for donor in range(nBox):
        if sum(donor):
            tMatrix[:,donor] = tMatrix[:,donor] / sum(donor)
    return tMatrix

if __name__ == "__main__":

    # set up graphical output handle
    g =  Gnuplot.Gnuplot()
    # g('set parametric')
    g.title('Normalised Size Spectra')
    g.xlabel('log10(mass; grams)')
    g.ylabel('log10(Number)')
#    g.size('0.5,0.5')

    # nMean--poisson mean n per sizeclass
    wVectorBase = arange(-12,2.5,.5)
    nBoxW = size(wVectorBase)
    nMean=3
    wLog10 = wVector(nMean, wVectorBase)
    w = 10**wLog10
    wX = 10**wVectorBase
    nBox = size(w)
    nDonorsMean=3 # the poisson average number of donors to each heterotroph 
    
    # determine community composition of "autotroph" vs "heterotroph"
    fractionHeterotroph = 0.25  # fraction of community that are heterotrophs
    fractionDetritivorus = 0.3  # fraction of heterotrophs (detritivorus)
    
    hMask = hMaskGet(fractionHeterotroph, nBox)
    aMask = abs(hMask - 1)
    dMask = dMaskGet(fractionDetritivorus, nBox, hMask)
    
    # initialise equal size biomass classes
    # biomass = ones(nBox)
    biomass = random(nBox)
    n = biomass / w

    # respiration intensity and rates
    rInt = (0.018, 0.140)
    rSlope = 0.75
    mineralisationRate = 0.05   # biomass-specific fixed at 5% per year
  
    mortalityRate = 0.1         # biomass-specific fixed at 10% per year

    rIntensity = rSpecific(w, rSlope, rInt)

    # losses to boundaries .. respiration and export flux
    boundaryOutflows = rIntensity  

    # boundary to autotroph flows
    # force boundary inputs to autotrophs to match their outputs (respiration)
    # and so "recipient" normalised -- assumes that these flows are dominant
    boundaryInflows = zeros(nBox)
    for consumer in range(nBox):
        if aMask[consumer]: boundaryInflows[consumer] = rIntensity[consumer]

    # recipient normalised detritivory rates
    detritusInflows = zeros(nBox)
    for consumer in range(nBox):
        if dMask[consumer]: detritusInflows[consumer] = random()
    
    # biomass flows to detritus (mortality+other losses)
    detritusOutflows = zeros(nBox) + mortalityRate

    detritusToBoundary = 0 # no detritus export / mineralisation
    boundaryToDetritus = 0 # no detritus import  

    # donor biomass-specific flows .. i.e., transfer Matrix (tMatrix)
    # initialise
    tMatrix = zeros((nBox,nBox),Float64)
    lowerLimit = 1e-10  # 1 e-8 is the lower bound to size of organisms
                        # (and therefore < 1 indivdual)
    time = 0
    nYears = 100
    tShape = tMatrix.shape
    extinct = ones(nBox)
    N = zeros(size(wVectorBase))
    
    
    while time < 365*nYears:

        time = time + 1

        nstd = standard_normal()  # the no. stdev for the error matrix
        eMatrix = nstd * standard_normal(tShape) # mean 0, std = nstd

        tMatrix = tMatrixMake(tMatrix, hMask, nDonorsMean)
        tC = tMatrix * ( 1.0 + eMatrix )  

        inputs  = matrixmultiply(tC, biomass)
        outputs = matrixmultiply(biomass, tC)
        
        # consumer biomass-specific flows 
        detritus = nstd*standard_normal(nBox)*detritusInflows - nstd*standard_normal(nBox)*detritusOutflows
        boundary = nstd*standard_normal(nBox)*boundaryInflows -  nstd*standard_normal(nBox)*boundaryOutflows
        
        dBiomass =  inputs - outputs + biomass*(detritus+boundary)
        extinct = where(biomass < lowerLimit, 0, 1)
        biomass = where(biomass < lowerLimit, 0, biomass + dBiomass/365.0 )

        rTotal = rIntensity * biomass
        richness = sum(extinct)
       
        R = sum(rTotal)
        B = sum(biomass)
        r_b = R/B
        
        bBase=zeros(nBoxW, Float)
        for i in range(nBoxW):
            for j in range(nBox):
                if wVectorBase[i] == wLog10[j]:
                    bBase[i] =  bBase[i] + biomass[j]
        
        nBase = bBase/wX
        nBase = where(nBase <= 1, 1, nBase)
        nBaseLog10 = log10(nBase)
        
        print time, richness, R, B, r_b
        g.plot(Gnuplot.Data(wVectorBase,nBaseLog10,  with='lines 2 2'))
        # g.replot()
            
