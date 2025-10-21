# define helper methods for flattening a dictionary
# of named parameters into a list of names and a list
# of parameter values
#from datetime import date

"""
Script to generate and write KD or KDUQ potential parameters for the interaction between a given nucleus and the proton and neutron

Usage:
  SampleKDUQ.py -A A -Z Z -E E 
  SampleKDUQ.py -A A -Z Z -E E --KD
  SampleKDUQ.py -h | --help

Options:
  -h --help         Show this screen.
  -A A              Integer mass of the nucleus (AMU)
  -Z Z              Charge of the nucleus
  -E E              Incidient lab energy, e.g. for 10 MeV/u, the deuteron energy is 20 MeV
  --KD              Flag for use of original KD rather than Bayesian KDUQ

"""


import numpy as np
import sys
import math
from pandas import DataFrame
from collections import OrderedDict
from json import load
from docopt import docopt

def flattenParameters(orderedDict):
    return [parameter for component in orderedDict.items() for parameter in component[1].values()]

def flattenParameterNames(orderedDict):
    return [component[0]+"_"+parameterName for component in orderedDict.items() for parameterName in component[1].keys()]

def readSample(potentialName, sampleNum):
    with open("{}/{}/parameters.json".format(potentialName, sampleNum)) as inputFile:
        return load(inputFile, object_pairs_hook=OrderedDict) # OrderedDict keeps the parameters ordered


def main(A,Z,E,UseKDUQ):
    print('Running for A = '+str(A)+', Z = '+str(Z)+', E = '+str(E))
    potentialName = "KDUQDemocratic" # name of desired potential ensemble
    nSamples = 416 # number of samples comprising KDUQDemocratic
    nflattenedSamples = nSamples
    
    # read the samples into a pandas DataFrame
    samples = [readSample(potentialName, sampleNum) for sampleNum in range(nSamples)]
    flattenedSamples = [flattenParameters(sample) for sample in samples]
    flattenedNames = flattenParameterNames(samples[0])
    
    #if KoningDelarocheInit:
        
    if not UseKDUQ:
        flattenedSamples[0][0]=59.3
        flattenedSamples[0][1]=21
        flattenedSamples[0][2]=0.024
        flattenedSamples[0][3]=0.007228
        flattenedSamples[0][4]=1.48e-6
        flattenedSamples[0][5]=1.994e-5
        flattenedSamples[0][6]=2.0e-8
        flattenedSamples[0][7]=0.007067
        flattenedSamples[0][8]=4.23e-6
        flattenedSamples[0][9]=1.729e-5
        flattenedSamples[0][10]=1.136e-8
        flattenedSamples[0][11]=7e-9
        flattenedSamples[0][12]=1.3039
        flattenedSamples[0][13]=0.4054
        flattenedSamples[0][14]=0.6778
        flattenedSamples[0][15]=1.487e-4
        flattenedSamples[0][16]=1.198
        flattenedSamples[0][17]=0.697
        flattenedSamples[0][18]=12.994
        flattenedSamples[0][19]=5.922
        flattenedSamples[0][20]=0.0030
        flattenedSamples[0][21]=0.0040
        flattenedSamples[0][22]= 1.1854
        flattenedSamples[0][23]=0.647
        flattenedSamples[0][24]=0.59
        flattenedSamples[0][25]=-3.1
        flattenedSamples[0][26]=160
        flattenedSamples[0][27]=12.195
        flattenedSamples[0][28]=0.0167
        flattenedSamples[0][29]=14.667
        flattenedSamples[0][30]=0.009629
        flattenedSamples[0][31]=73.55
        flattenedSamples[0][32]=0.0795
        flattenedSamples[0][33]=16.0
        flattenedSamples[0][34]=16.0
        flattenedSamples[0][35]=0.0180
        flattenedSamples[0][36]=0.003802
        flattenedSamples[0][37]=8.0
        flattenedSamples[0][38]=156.0
        flattenedSamples[0][39]=11.5
        flattenedSamples[0][40]=1.3424
        flattenedSamples[0][41]=0.01585
        flattenedSamples[0][42]=0.5446
        flattenedSamples[0][43]=1.656e-4
        flattenedSamples[0][44]=0.5187
        flattenedSamples[0][45]=5.205e-4
        nflattenedSamples=1
    strnucleon=[False,True]
    
    for neutron in strnucleon:
            N=A-Z
            alpha=(N-Z)/A
            if neutron:
                sign=-1
                Ef=-11.2813+0.02646*A
                string='n'
            else:
                sign=1
                Ef=-8.4075+0.01378*A
                string='p'
            DE=E-Ef
            #if KoningDelarocheInit:
            EkeV = int(E*1000.)
            if not UseKDUQ:
                filename='outputs/PotentialsKD_{string}_A{A}_Z{Z}_E{E}.pot'\
                      .format(string=string,A=A,Z=Z,E=EkeV,pulls=nflattenedSamples)
            else:
                filename='outputs/PotentialsKDUQ_{string}_A{A}_Z{Z}_E{E}.pot'\
                      .format(string=string,A=A,Z=Z,E=EkeV,pulls=nflattenedSamples)
    
            myfilepot=open(filename,'w')
    
            for i in range(0,nflattenedSamples):
                v1=flattenedSamples[i][0]+sign*flattenedSamples[i][1]*alpha-flattenedSamples[i][2]*A
                if neutron:
                    v20=flattenedSamples[i][3]
                    v2A=flattenedSamples[i][4]
                    v30=flattenedSamples[i][5]
                    v3A=flattenedSamples[i][6]
                    w10=flattenedSamples[i][27]
                    w1A=flattenedSamples[i][28]
                    aD0=flattenedSamples[i][42]
                    aDA=flattenedSamples[i][43]
                else:
                    v20=flattenedSamples[i][7]
                    v2A=flattenedSamples[i][8]
                    v30=flattenedSamples[i][9]
                    v3A=flattenedSamples[i][10]
                    w10=flattenedSamples[i][29]
                    w1A=flattenedSamples[i][30]
                    aD0=flattenedSamples[i][44]
                    aDA=flattenedSamples[i][45]
    
                v2=v20+sign*v2A*A
                v3=v30+sign*v3A*A
                v4=flattenedSamples[i][11]
                w1=w10+w1A*A
                w2=flattenedSamples[i][31]+flattenedSamples[i][32]*A
                d1=flattenedSamples[i][33]+sign*flattenedSamples[i][34]*alpha
                if((A-flattenedSamples[i][38])/flattenedSamples[i][37])> 3e2:
                    d2=flattenedSamples[i][35]
                else:
                    d2=flattenedSamples[i][35]+flattenedSamples[i][36]/(1+math.exp((A-flattenedSamples[i][38])/flattenedSamples[i][37]))
                d3=flattenedSamples[i][39]
                vso1=flattenedSamples[i][19]+flattenedSamples[i][20]*A
                vso2=flattenedSamples[i][21]
                wso1=flattenedSamples[i][25]
                wso2=flattenedSamples[i][26]                                            
    
                rV=flattenedSamples[i][12]-flattenedSamples[i][13]*pow(A,-1/3)
                aV=flattenedSamples[i][14]-flattenedSamples[i][15]*A
                rD=flattenedSamples[i][40]-flattenedSamples[i][41]*pow(A,1/3)
                aD=aD0+sign*aDA*A
                rSO=flattenedSamples[i][22]-flattenedSamples[i][23]*pow(A,-1/3)
                aSO=flattenedSamples[i][24]
                if neutron:
                    rc=1
                else:
                    rc=flattenedSamples[i][16]+flattenedSamples[i][17]*pow(A,-2/3)+flattenedSamples[i][18]*pow(A,-5/3)
    
                VC=1.73*Z*pow(A,-1/3)/rc
                                                       
                Vv=v1*(1-v2*DE+v3*pow(DE,2)-v4*pow(DE,3))
                if not neutron:
                    Vv=Vv+VC*v1*(v2-2*v3*DE+3*v4*pow(DE,2))
                Wv=w1*pow(DE,2)/(pow(DE,2)+pow(w2,2))
                WD=d1*(pow(DE,2)/(pow(DE,2)+pow(d3,2)))*math.exp(-d2*DE)
                VSO=vso1*math.exp(-vso2*DE)
                WSO=wso1*pow(DE,2)/(pow(DE,2)+pow(wso2,2))
    
    
                str_pot='{nb} {VR} {rR} {aR} {WI} {rI} {aI} {WS} {rS} {aS} {VSO} {rSO} {aSO} {WSO} {rSO} {aSO} {rC} \n'.format(nb=i,\
                            VR="%3.3f"%Vv,WI="%3.3f"%Wv,WS="%3.3f"%WD,\
                            rR="%3.3f"%rV,rI="%3.3f"%rV,rS="%3.3f"%rD,\
                            aR="%3.3f"%aV,aI="%3.3f"%aV,aS="%3.3f"%aD,\
                            VSO="%3.3f"%VSO,WSO="%3.3f"%WSO,rSO="%3.3f"%rSO,aSO="%3.3f"%aSO,rC="%3.3f"%rc)
                myfilepot.write(str_pot)
            myfilepot.close()
    
    
if __name__ == "__main__":
    args = docopt(__doc__)

    A = int(args['-A'])
    Z = int(args['-Z'])
    E = float(args['-E'])

    if args['--KD']:
        UseKDUQ = False
        print('Using original Koning Delaroche parameters')
    else:
        UseKDUQ = True
        print('Using uncertainty-parameterized Koning Delaroche parameters')

    main(A, Z, E, UseKDUQ)
 
