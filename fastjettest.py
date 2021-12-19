#!/usr/bin/env python3

"""Simple example to try out fastjet from python with ROOT / cppyy

"""

# Import ROOT, needed fastjet headers and libs
import ROOT
ROOT.gInterpreter.AddIncludePath( "/home/iwsatlas1/skluth/qcd/fastjet/fastjet-3.4.0/install/include/" )
ROOT.gInterpreter.ProcessLine( '#include "fastjet/JetDefinition.hh"' )
ROOT.gInterpreter.ProcessLine( '#include "fastjet/ClusterSequence.hh"' )
ROOT.gInterpreter.ProcessLine( '#include "fastjet/PseudoJet.hh"' )
ROOT.gInterpreter.ProcessLine( '#include "fastjet/JadePlugin.hh"' )
ROOT.gInterpreter.ProcessLine( '#include "fastjet/SISConeSphericalPlugin.hh"' )
ROOT.gInterpreter.ProcessLine( '#include "fastjet/PxConePlugin.hh"' )
ROOT.gSystem.Load( "libfastjet.so" )
ROOT.gSystem.Load( "libfastjetplugins.so" )
ROOT.gSystem.Load( "libfastjettools.so" )
ROOT.gSystem.Load( "/usr/lib/x86_64-linux-gnu/libgfortran.so.5" )

# Make fastjet symbols available
from ROOT.fastjet import JetDefinition, ClusterSequence, PseudoJet
from ROOT.fastjet import ee_kt_algorithm, JadePlugin, SISConeSphericalPlugin, PxConePlugin, sorted_by_E

# Pass a recombiner subclass to fastjet::JetDefinition:
ROOT.gInterpreter.ProcessLine( '#include "EEE0Recombiner.hh"' )
from ROOT import EEE0Recombiner

def main( algo="jade", njet=3 ):

    # get the e+e- event
    filename= "single-ee-event.dat"
    event= read_event( filename )

    # Event 4-vector sum
    pjetsum= PseudoJet( 0.0, 0.0, 0.0, 0.0 )
    for pjet in event:
        pjetsum+= pjet
    print( "Event 4-vector sum" )
    print( f"{'px':>10s} {'py':>10s} {'pz':>10s} {'E':>10s}" )
    print( f"{pjetsum[0]:10.3f} {pjetsum[1]:10.3f} {pjetsum[2]:10.3f} {pjetsum[3]:10.3f}" )
    Evis= pjetsum[3]
    
    # set up jet definition
    if "durham" in algo:
        jetDef= JetDefinition( ee_kt_algorithm )
    elif "jade" in algo:
        plugin= JadePlugin()
        jetDef= JetDefinition( plugin )
        rc= EEE0Recombiner()
        jetDef.set_recombiner( rc )
    elif "eesiscone" in algo:
        R= 0.6
        Emin= 0.0
        plugin= SISConeSphericalPlugin( R, 0.75 )
        jetDef= JetDefinition( plugin )
    elif "pxcone" in algo:
        R= 0.6
        Emin= 0.01*Evis
        plugin= PxConePlugin( R, Emin, 0.75, True, 1 )
        jetDef= JetDefinition( plugin )
    print( "jet definition is:", jetDef.description() )
    
    # cluster it
    clusseq= ClusterSequence( event, jetDef )

    # print out some information about the event and clustering
    print( "Event has {0} particles".format( len(event) ) )
    if "durham" in algo or "jade" in algo:
        clusseqjets= clusseq.exclusive_jets( njet )
    elif "eesiscone" or "pxcone" in algo:
        clusseqjets= clusseq.inclusive_jets( Emin )
    jets= sorted_by_E( clusseqjets )
        
    # print the jets
    print_jets( jets )
    
    # print internal information about the jets
    for ijet, jet in enumerate(jets):
        print( f"Number of constituents of jet {ijet} is {len( jet.constituents() )}" ) 

    # The End
    return

# Print jet 4-vectors
def print_jets( jets ):
    print( f"{'jet #':>5s} {'px':>10s} {'py':>10s} {'pz':>10s} {'E':>10s}" )
    for ijet, jet in enumerate( jets ):
        print( f"{ijet:5d} {jet.px():10.3f} {jet.py():10.3f} {jet.pz():10.3f} {jet.E():10.3f}" )
    return

# Read fastjet example event
def read_event( filename ):
    f= open( filename, 'r' )
    event= ROOT.std.vector(PseudoJet)()
    for line in f:
        p= line.split()
        event.push_back( PseudoJet( float(p[0]), float(p[1]), float(p[2]), float(p[3]) ) )
    return event
    
if __name__ == '__main__':
    main()
