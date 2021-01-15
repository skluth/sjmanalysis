#!/usr/bin/env python

from array import array

import pickle as pkl

import numpy as np

from math import sqrt

import yclus


def combine( vi, vj ):
    return np.add( vi, vj )

def invariantMassSq( v ):
    return v[3]**2 - v[0]**2 - v[1]**2 - v[2]**2

def calcEvis( objects ):
    sumE= 0.0
    for obj in objects:
        sumE+= obj[3]
    return sumE
    
def distance( objects, i, j ):
    vi= objects[i]
    vj= objects[j]
    vimod= sqrt( vi[0]**2 + vi[1]**2 + vi[2]**2 )
    vjmod= sqrt( vj[0]**2 + vj[1]**2 + vj[2]**2 )
    costhetaij= ( vi[0]*vj[0] + vi[1]*vj[1] + vi[2]*vj[2] )/ ( vimod*vjmod )
    return 2.0*vi[3]*vj[3]*(1.0-costhetaij)

def calcAllYijandFindMin( objects, evis ):
    nobj= len(objects)
    minyij= 2.0*evis**2
    imin= -1
    jmin= -1
    for i in range( nobj-1 ):
        for j in range( i+1, nobj ):
            yij= distance( objects, i, j )
            if yij < minyij:
                minyij= yij
                imin= i
                jmin= j
    if imin == -1 or jmin == -1:
        raise RuntimeError( "no minimum yij found" )
    return min( imin, jmin ), max( imin, jmin ), minyij/evis**2

def JADEAlgorithm( objects ):
    # Find pairs with smallest distance, remove, add combination
    ymergeValues= []
    while( len( objects ) > 1 ):
        # print( "JADEAlgorithm: ", len(objects) )
        # find indices of pair with smallest distance and yij value
        evis= calcEvis( objects )
        i, j, yij= calcAllYijandFindMin( objects, evis )
        # print( "JADEAlgorithm: pair with smallest distance:", i, j, yij  )
        ymergeValues.append( yij )
        # Remove pair and add combined object
        obji= objects.pop( i )
        objj= objects.pop( j-1 )
        newobj= combine( obji, objj )
        objects.append( newobj )
    nynm= len(ymergeValues)
    return ymergeValues[nynm-4:nynm]

def readData( filename="da91_96.pkl", maxevt=10 ):
    with open( filename, "rb" ) as pklFile:
        data= pkl.load( pklFile )
    nevt= 0
    for key in data:
        if nevt >= maxevt:
            break
        nevt+= 1
        eventRecord= data[key]
        objects= eventRecord["objects"]
        nobjs= len(objects)
        
        fobjects= np.zeros( (4,nobjs) )
        for i in range( 4 ):
            for j in range( nobjs ):
                fobjects[i,j]= objects[j][i]
        
        ymergevalues= eventRecord["ynm"]
        print( "event number key:", nevt, key )

        JADEymergevalues= JADEAlgorithm( objects )
        print( "JADE yij values:", JADEymergevalues )
        
        print( "Ynm merge values:", ymergevalues )

        nt= len(objects)
        ierryk= yclus.ykern( 1, fobjects )
        if ierryk == 0:
            print( "YKERN yij: ", end="" )
            ptr, ierryt= yclus.ytree( False )
            if ierryt == 0:
                for i in range( 4 ):
                    print( ptr[6,i], " ", end="" )
                print( "" )
            else:
                print( "YTREE: ierr", ierryt )
        else:
            print( "YKERN ierr", ierryk )
        
    return


def main():
    readData()
    return


if __name__ == '__main__':
   main()

