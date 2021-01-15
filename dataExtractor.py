#!/usr/bin/env python

from ROOT import *

from array import array

import pickle as pkl

gROOT.LoadMacro( "libNtupleReaderDict.so" )


def writeData( filename="da91_96_200.root", outfilename="da91_96.pkl", maxevt=100 ):
    outData= {}
    ntr= LEP1NtupleReader( filename )
    ievent= 1
    while ntr.GetEvent( ievent ) and ievent < maxevt:
        ievent+= 1
        if ntr.Selection( "91.2" ):
            tlvList= []
            vtlv= ntr.GetLorentzVectors( "mt" )
            for tlv in vtlv:
                tlv.Print()
                a= array( "f" )
                for i in range( 4 ):
                    a.append( tlv[i] )
                tlvList.append( a )
            ymergeValues= array( "f" )
            for ijet in range( 2, 5 ):
                ymerge= ntr.getYmerge( "jade", "mt", ijet )
                ymergeValues.append( ymerge )
            eventRecord= { "objects": tlvList, "ynm": ymergeValues }
            outData[ievent]= eventRecord
    with open( outfilename, "wb" ) as pklFile:
        pkl.dump( outData, pklFile )
    pklFile.close()
    print( "writeData:", ievent, "events" )
    return

def main():
    writeData()
    return


if __name__ == '__main__':
   main()

