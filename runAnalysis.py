#!/usr/bin/env python

from ROOT import *
gInterpreter.ProcessLine( '#include "LEP1NtupleReader.hh"' )
gInterpreter.ProcessLine( '#include "LEP2NtupleReader.hh"' )
gInterpreter.ProcessLine( '#include "SjmConfigParser.hh"' )
gInterpreter.ProcessLine( '#include "Analysis.hh"' )
gInterpreter.ProcessLine( '#include "ObservableFactory.hh"' )
gInterpreter.ProcessLine( '#include "FilledObservable.hh"' )
gInterpreter.ProcessLine( '#include "OutputWriter.hh"' )
gROOT.LoadMacro( "libNtupleReader.so" )

def createNtupleReader( filename, ecms ):
    lep2ecms= [ "130", "136", "161", "172", "183", "189", 
                "192", "196", "200", "202", "205", "207" ]
    if ecms == "91.2":
        print( "createNtupleReader: LEP1NtupleReader" )
        result= LEP1NtupleReader( filename )
    elif ecms in lep2ecms:
        print( "createNtupleReader: LEP2NtupleReader" )
        result= LEP2NtupleReader( filename )
    else:
        raise ValueError( "createNtupleReader: wrong ecms "+ecms )
    return result

def processAnalyses( analyses, vobs, filename, ecm, maxevt ):
    print( "processAnalyses: file "+filename+", analyses:" )
    for analysis in analyses:
        print( analysis.getTag() )
    ntr= createNtupleReader( filename, ecm )
    ievent= 1
    while ntr.GetEvent( ievent ) and ievent < maxevt:
        ievent+= 1
        selections= ntr.getSelections( ecm )
        MCnonrad= ntr.MCNonRad()
        for analysis in analyses:
            cuts= analysis.getCuts()
            mccuts= analysis.getMccuts()
            if( ( cuts == "none" or selections.at( cuts ) ) and
                ( mccuts == "none" or MCnonrad ) ):
                for obs in vobs:
                    # Not all observables have all analysis variants
                    # due to filling of transfer matrices:
                    if obs.containsAnalysis( analysis ):
                        obs.fill( ntr, analysis )
    print( "processAnalyses:", ievent-1, "events" )
    return

def processUnfolding( measuredAnalyses, unfoldsource, vobs ):
    print( "processUnfolding: bin-by-bin unfolding for analyses:" )
    hadronlevel= Analysis( unfoldsource, "hadron", "none", "nonrad" )
    print( "Hadron level:", hadronlevel.getTag() )
    for measured in measuredAnalyses:
        measuredMC= Analysis( measured )
        measuredMC.setSource( unfoldsource )
        print( measured.getTag(), measuredMC.getTag() )
        unfolder= BbbUnfolder( measured, measuredMC, hadronlevel )
        for obs in vobs:
            unfolder.unfold( obs )
    return

def getFilled( vobs ):
    vfobs= vector(FilledObservable)()
    for obs in vobs:
        vfobspart= obs.getFilledObservables()
        for fobs in vfobspart:
            vfobs.push_back( fobs )        
    return vfobs

def LEP1Analysis( cfgFile="sjmconfig_91_96.cfg" ):

    # Read configs:
    sjmcp= SjmConfigParser( cfgFile )
    
    # Define analysis variations:
    measuredAnalyses= vector(Analysis)()
    measuredAnalyses.push_back( Analysis( "data", "mt", "stand" ) )
    measuredAnalyses.push_back( Analysis( "data", "mt", "costt07" ) )
    # measuredAnalyses.push_back( Analysis( "data", "mt", "nch7" ) )
    measuredAnalyses.push_back( Analysis( "data", "tc", "stand" ) )
    pyAnalyses= vector(Analysis)()
    pyAnalyses.push_back( Analysis( "py", "mt", "stand" ) )
    pyAnalyses.push_back( Analysis( "py", "mt", "costt07" ) )
    # pyAnalyses.push_back( Analysis( "py", "mt", "nch7" ) )
    pyAnalyses.push_back( Analysis( "py", "tc", "stand" ) )
    pyAnalyses.push_back( Analysis( "py", "hadron", "none", "nonrad" ) )
    hwAnalyses= vector(Analysis)()
    hwAnalyses.push_back( Analysis( "hw", "mt", "stand" ) )
    hwAnalyses.push_back( Analysis( "hw", "hadron", "none", "nonrad" ) )
    allAnalyses= vector(Analysis)( measuredAnalyses )
    allAnalyses.insert( allAnalyses.end(), pyAnalyses.begin(), pyAnalyses.end() )
    allAnalyses.insert( allAnalyses.end(), hwAnalyses.begin(), hwAnalyses.end() )

    # Define observables:
    obsnames= vector(string)()
    obsnames.push_back( "lepthrust" )
    obsnames.push_back( "partonshower" )
    obsnames.push_back( "durhamymerge23" )
    obsnames.push_back( "jadeymerge23" )
    obsnames.push_back( "durhamymergefj" )
    obsnames.push_back( "jadeymergefj" )
    obsnames.push_back( "durhamycutfj" )
    obsnames.push_back( "jadeycutfj" )
    obsnames.push_back( "antiktemin" )
    obsnames.push_back( "antiktR" )
    obsnames.push_back( "sisconeemin" )
    obsnames.push_back( "sisconeR" )
    obsfac= ObservableFactory( sjmcp )
    vobs= obsfac.createObservables( obsnames, allAnalyses )

    # Add extras for migration matrices where needed:
    pyMatrixExtras= vector(Analysis)()
    pyMatrixExtras.push_back( Analysis( "py", "hadron", "stand", "nonrad" ) )
    pyMatrixExtras.push_back( Analysis( "py", "mt", "stand", "nonrad", "hadron" ) )
    pyAnalyses.insert( pyAnalyses.end(), pyMatrixExtras.begin(), pyMatrixExtras.end() )
    hwMatrixExtras= vector(Analysis)()
    hwMatrixExtras.push_back( Analysis( "hw", "hadron", "stand", "nonrad" ) )
    hwMatrixExtras.push_back( Analysis( "hw", "mt", "stand", "nonrad", "hadron" ) )
    hwAnalyses.insert( hwAnalyses.end(), hwMatrixExtras.begin(), hwMatrixExtras.end() )
    for obs in vobs:
        name= obs.getName()
        if( name == "thrust" or name == "durhamymerge23" or
            name == "jadeymerge23" or name == "partonshower" ):
            obs.addAnalyses( pyMatrixExtras )
            obs.addAnalyses( hwMatrixExtras )

    # Fill from data and mc (PYTHIA and HERWIG) ntuples:
    # get first files only, otherwise use getFilepath
    ecm= sjmcp.getItemString( "General.energy" )
    maxevt= sjmcp.getItemInt( "General.maxevt" )
    datafilename= sjmcp.getItemString( "Data.files" )
    processAnalyses( measuredAnalyses, vobs, datafilename, ecm, maxevt )
    pyfilename= sjmcp.getItemString( "Signal.files" )
    processAnalyses( pyAnalyses, vobs, pyfilename, ecm, maxevt )
    hwfilename= sjmcp.getItemString( "AltSignal.files" )
    processAnalyses( hwAnalyses, vobs, hwfilename, ecm, maxevt )

    # Get FilledObservables for further processing:
    vfobs= getFilled( vobs )

    # Unfolding bin-by-bin:
    # PYHTHIA based:
    processUnfolding( measuredAnalyses, "py", vfobs )
    # HERWIG based for systematic:
    measuredAnalysesHw= vector(Analysis)()
    measuredAnalysesHw.push_back( Analysis( "data", "mt", "stand" ) )
    processUnfolding( measuredAnalysesHw, "hw", vfobs )
    # MC detector level with MC as cross check for PYTHIA and HERWIG:
    measuredPyAnalyses= vector(Analysis)()
    measuredPyAnalyses.push_back( Analysis( "py", "mt", "stand" ) )
    measuredPyAnalyses.push_back( Analysis( "py", "mt", "costt07" ) )
    # measuredPyAnalyses.push_back( Analysis( "py", "mt", "nch7" ) )
    measuredPyAnalyses.push_back( Analysis( "py", "tc", "stand" ) )
    processUnfolding( measuredPyAnalyses, "py", vfobs )
    measuredHwAnalyses= vector(Analysis)()
    measuredHwAnalyses.push_back( Analysis( "hw", "mt", "stand" ) )
    processUnfolding( measuredHwAnalyses, "hw", vfobs )

    # Normalise and calculate stat errors, print
    # Normalisation only during postprocessing
    for fobs in vfobs:
        # vfobs[i]->finalise();
        #fobs.print()
        pass

    # Write root objects (TH1D or TGraphErrors, and TH2D):
    outfileName= sjmcp.getItemString( "General.outfile" )
    writer= OutputWriter( outfileName )
    # writer= OutputWriter( "LEP1Analysis.root" )
    writer.write( vfobs )

    # The End:
    return


def main():
    LEP1Analysis()
    return


if __name__ == '__main__':
   main()

