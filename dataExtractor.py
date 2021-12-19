#!/usr/bin/env python3

# run with e.g. LD_LIBRARY_PATH=/home/skluth/qcd/sjmanalysis:/home/skluth/Downloads/root/root_v6.22.06/lib:/home/skluth/qcd/hepmc/hepmc2.06.09/install/lib

from ROOT import *

gInterpreter.ProcessLine( '#include "LEP1NtupleReader.hh"' )
gROOT.LoadMacro( "libNtupleReader.so" )

from array import array

import pickle as pkl

# def writeData( filename="da91_96_200.root", outfilename="da91_96.pkl", maxevt=100 ):
def writeData( filename="mc5025_1_200.root", outfilename="mc5025_1_200.pkl", maxevt=100 ):
    # Will contain all events as "event records"
    outData= {}
    # Reads an ntuple
    ntr= LEP1NtupleReader( filename )
    ievent= 1
    # Loop over events in ntuple
    while ntr.GetEvent( ievent ) and ievent < maxevt:
        ievent+= 1
        datatlvList= []
        partontlvList= []
        hadrontlvList= []
        dataymergeValues= array( "f" )
        partonymergeValues= array( "f" )
        hadronymergeValues= array( "f" )
        # Check if selected at detector level
        txt= std.string
        if ntr.Selection( std.string("91.2") ):
            copyTlvToArrays( datatlvList, ntr, "mt" )
            copyYnm( dataymergeValues, ntr, "mt" )
        # Parton and hadron level for all events
        copyTlvToArrays( partontlvList, ntr, "parton" )
        copyYnm( partonymergeValues, ntr, "parton" )
        copyTlvToArrays( hadrontlvList, ntr, "hadron" )
        copyYnm( hadronymergeValues, ntr, "hadron" )
        # Assemble event record and add to data structure
        eventRecord= { "objects": datatlvList,
                       "ynm": dataymergeValues, 
                       "partons": partontlvList,
                       "partonynm": partonymergeValues, 
                       "hadrons": hadrontlvList,
                       "hadronynm": hadronymergeValues }
        outData[ievent]= eventRecord
    with open( outfilename, "wb" ) as pklFile:
        pkl.dump( outData, pklFile )
    pklFile.close()
    print( "writeData:", ievent, "events" )
    return

def copyYnm( ymergeValues, ntr, level ):
    for ijet in range( 2, 5 ):
        ymerge= ntr.getYmerge( "jade", level, ijet )
        ymergeValues.append( ymerge )
    return

def copyTlvToArrays( tlvList, ntr, level ):
    vtlv= ntr.GetLorentzVectors( level )
    for tlv in vtlv:
        tlv.Print()
        a= array( "f" )
        for i in range( 4 ):
            a.append( tlv[i] )
        tlvList.append( a )
    return


def main():
    writeData()
    return


if __name__ == '__main__':
   main()

