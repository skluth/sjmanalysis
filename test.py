
from ROOT import TFile, gROOT, gPad, TVectorD, TObject

from ROOT import TGraphErrors, TH1D, TLegend, TCanvas, TLatex, SetOwnership

# From aholzner.wordpress.com/tag/pyroot/
# T* subclasses to handle (disable!) python ownership
oldTGraphErrorsinit= TGraphErrors.__init__
class GarbageCollectionResistentTGraphErrors( TGraphErrors ):
  def __init__( self, *args ):
    oldTGraphErrorsinit( self, *args )
    SetOwnership( self, False )
TGraphErrors= GarbageCollectionResistentTGraphErrors

oldTH1Dinit= TH1D.__init__
class GarbageCollectionResistentTH1D( TH1D ):
  def __init__( self, *args ):
    oldTH1Dinit( self, *args )
    SetOwnership( self, False )
TH1D= GarbageCollectionResistentTH1D

oldTLegendinit= TLegend.__init__
class GarbageCollectionResistentTLegend( TLegend ):
  def __init__( self, *args ):
    oldTLegendinit( self, *args )
    SetOwnership( self, False )
TLegend= GarbageCollectionResistentTLegend

oldTCanvasinit= TCanvas.__init__
class GarbageCollectionResistentTCanvas( TCanvas ):
  def __init__( self, *args ):
    oldTCanvasinit( self, *args )
    SetOwnership( self, False )
TCanvas= GarbageCollectionResistentTCanvas

oldTLatexinit= TLatex.__init__
class GarbageCollectionResistentTLatex( TLatex ):
  def __init__( self, *args ):
    oldTLatexinit( self, *args )
    SetOwnership( self, False )
TLatex= GarbageCollectionResistentTLatex

gROOT.LoadMacro( "libAnalysisDict.so" )
from ROOT import Analysis, TH1DAnalysisObject, TGEAnalysisObject

from array import array

import numpy as np

# Read numbers columnwise from ascii txt files into arrays indexed by column number:
def ascii2arrays( filename ):
    lines= [ line.rstrip( '\n' ) for line in open( filename ) ]
    arrays= dict()
    for line in lines:
        tokens= line.split()
        for itoken in range( len( tokens ) ):
            if not itoken in arrays:
                arrays[itoken]= array( "d" )
            arrays[itoken].append( float( tokens[itoken] ) )
    return arrays

# Factory method to create AnalysisObject instances
def getAnalysisObjectFromFile( tfile, obs, analysis ):
    ao= None
    key= obs+" "+analysis.getTag()+";1"
    obj= tfile.Get( key )
    if not obj:
        raise RuntimeError( "getAnalysisObjectFromFile: AnalysisObject with key "+key+" not in file "+tfile.GetName() )
    if obj.ClassName() == "TH1D":
        errobj= tfile.Get( "errm "+obs+" "+analysis.getTag() )
        if errobj:
            ao= TH1DAnalysisObject( obj, errobj )
        else:
            ao= TH1DAnalysisObject( obj )
    elif obj.ClassName() == "TGraphErrors":
        ao= TGEAnalysisObject( obj )
    else:
        raise RuntimeError( "getAnalysisObjectFromFile: can't handle class name"+obj.ClassName() )
    return ao

# Interface for analyses
class AnalysisObservable:

    def __init__( self, name ):
        self.obs= name
        self.aostand=None
        self.points=None
        self.values=None
        self.sterrs=None
        self.syerrs=None
        self.variationsDelta=None
        self.events= dict()
        self.rawEvents= dict()
        return

    def setupStandardAnalysis( self, standardAnalysis, tfile ):
        self.aostand= getAnalysisObjectFromFile( tfile, self.obs, standardAnalysis )
        self.points= array( "d", self.aostand.getPoints() )
        self.values= array( "d", self.aostand.getValues() )
        self.sterrs= array( "d", self.aostand.getErrors() )
        self.events["stand"]= self.aostand.getNEvents()
        return
    
    def subtractVariations( self, analysisVariations, tfile ):
        self.variationsDelta= dict()
        for key in analysisVariations.keys():
            ao= getAnalysisObjectFromFile( tfile, self.obs, analysisVariations[key] )
            variationData= array( "d", ao.getValues() )
            self.variationsDelta[key]= np.subtract( variationData, self.values )
            self.events[key]= ao.getNEvents()
        return
       
    def calcSystSumSq( self, keys ):
        self.syerrs= 0.0
        for key in keys:
            self.syerrs+= np.square( self.variationsDelta[key] )
        self.syerrs= np.sqrt( self.syerrs )
        return
    
    def printResults( self, width=7, precision=3, pointwidth=4, pointprec=2, opt="?" ):
        print( "Results for", self.obs )
        print( self.aostand.getPointLabel( pointwidth ), end="" )
        fmt= "{:>"+str(width)+"}"
        for key in [ "val", "stat", "sys" ]:
            print( fmt.format( key ), end="" )
        if "d" in opt:
            for key in sorted( self.variationsDelta.keys() ):
                print( fmt.format( key ), end="" )
        print()
        if "m" in opt:
            sterrs= self.aostand.getErrors( "m" )
        else:
            sterrs= self.sterrs            
        fmt="{:"+str(width)+"."+str(precision)+"f}"
        for i in range( len( self.values ) ):
            if self.obs.find( "EEC" ) >= 0 and i < len( self.values )-1:
                rad2grad= 180.0/3.14159
                leftedge= self.points[i]*rad2grad
                rightedge= self.points[i+1]*rad2grad
                print( "{0:3.0f} {1:3.0f}  ".format( leftedge, rightedge ), end="" )
            else:
                print( self.aostand.getPointStr( i, pointwidth, pointprec ), end="" )
            print( fmt.format( self.values[i] ), end="" )
            print( fmt.format( sterrs[i] ), end="" )
            print( fmt.format( self.syerrs[i] ), end="" )
            if "d" in opt:
                for key in sorted( self.variationsDelta.keys() ):
                    print( fmt.format( self.variationsDelta[key][i] ), end="" )
            print()
        return

    def printErrors( self, width=7, precision=4 ):
        from math import sqrt
        errorMatrix= self.aostand.getErrorMatrix()
        fmt="{:"+str(width)+"."+str(precision)+"f}"
        for i in range( len( self.sterrs )-1 ):
            binw= self.points[i+1]-self.points[i]
            diagError= sqrt( errorMatrix(i,i) )/binw
            print( fmt.format( self.sterrs[i] ), fmt.format( diagError ) )
        return

    def plot( self, plotoptions, opt="?" ):
        vx= array( "d", self.aostand.getPointsCenter() )
        values= self.values
        sterrs= self.sterrs
        if "m" in opt:
            print( "AnalysisObservable::plot: use errors from error matrix" )
            sterrs= array( "d", self.aostand.getErrors( "m" ) )
        syerrs= self.syerrs
        npoints= len(vx)
        if "xshift" in plotoptions:
            for i in range(npoints):
                vx[i]+= plotoptions["xshift"]
        vex= array( "d", npoints*[0.0] )
        tgest= TGraphErrors( npoints, vx, values, vex, sterrs )
        #SetOwnership( tgest, False )
        toterrs= np.sqrt( np.add( np.square( sterrs ),  np.square( syerrs ) ) )
        tgesy= TGraphErrors( npoints, vx, values, vex, toterrs )
        #SetOwnership( tgesy, False )
        tgesy.SetMarkerStyle( plotoptions["markerStyle"] )
        tgesy.SetMarkerSize( plotoptions["markerSize"] )
        drawas= plotoptions["drawas"] if "drawas" in plotoptions else "p"
        tgesy.SetName( self.obs )
        if "fillcolor" in plotoptions:
            tgesy.SetFillColor(plotoptions["fillcolor"])
            tgest.SetFillColor(plotoptions["fillcolor"])
        if "s" in opt:
            tgesy.Draw( "psame" )
        else:
            if "title" in plotoptions:
                tgesy.SetTitle( plotoptions["title"] )
            else:
                tgesy.SetTitle( self.obs )
            tgesy.SetMinimum( plotoptions["ymin"] )
            tgesy.SetMaximum( plotoptions["ymax"] )
            xaxis= tgesy.GetXaxis()
            xaxis.SetLimits( plotoptions["xmin"], plotoptions["xmax"] )
            if "xlabel" in plotoptions:
                xaxis.SetTitle( plotoptions["xlabel"] )
            if "ylabel" in plotoptions:
                tgesy.GetYaxis().SetTitle( plotoptions["ylabel"] )
            tgesy.Draw( "a"+drawas )
        optlogx= plotoptions["logx"] if "logx" in plotoptions else 0
        gPad.SetLogx( optlogx )
        optlogy= plotoptions["logy"] if "logy" in plotoptions else 0
        gPad.SetLogy( optlogy )
        tgest.Draw( "same"+drawas )
        return tgest, tgesy
    
    def maxAbsErrorSq( self, errorKey1, errorKey2 ):
        return np.square( np.maximum( np.absolute( self.variationsDelta[errorKey1] ),
                                      np.absolute( self.variationsDelta[errorKey2] ) ) )
    
    def printEvents( self ):
        for key in sorted( self.events.keys() ):
            print( key, self.events[key] )
        return
    def printRawEvents( self ):
        for key in sorted( self.rawEvents.keys() ):
            print( key, self.rawEvents[key] )
        return

    def readRawEvents( self, standardAnalysis, analysisVariations, tfile, srclist=[] ):
        allAnalyses= analysisVariations.copy()
        allAnalyses["stand"]= standardAnalysis
        for source in [ "data", "py" ]+srclist:
            for key in allAnalyses.keys():
                analysis= allAnalyses[key]
                rawAnalysis= Analysis( source, analysis.getReco(), analysis.getCuts() )
                ao= getAnalysisObjectFromFile( tfile, self.obs, rawAnalysis )
                self.rawEvents[rawAnalysis.getTag()]= ao.getNEvents()
        hwRawAnalysis= Analysis( "hw", "mt", "stand" )
        ao= getAnalysisObjectFromFile( tfile, self.obs, hwRawAnalysis )
        self.rawEvents[hwRawAnalysis.getTag()]= ao.getNEvents()
        return
    
    
# LEP1 Analysis:
class LEP1AnalysisObservable( AnalysisObservable ):

    def __init__( self, obs ):
        AnalysisObservable.__init__( self, obs )
        return

    def setupFromFile( self, tfile, unf="bbb" ):
        standardAnalysis= Analysis( "data mt stand none none none py " + unf )
        analysisVariations= {
            "tc": Analysis( "data tc stand none none none py " + unf ),
            "costt07": Analysis( "data mt costt07 none none none py " + unf ),
            "hw": Analysis( "data mt stand none none none hw " + unf ) }
        self.setupStandardAnalysis( standardAnalysis, tfile )
        self.subtractVariations( analysisVariations, tfile )
        self.calcSystSumSq( analysisVariations.keys() )
        self.readRawEvents( standardAnalysis, analysisVariations, tfile )
        return

# LEP1.5 Analysis:
class LEP15AnalysisObservable( AnalysisObservable ):

    def __init__( self, obs ):
        AnalysisObservable.__init__( self, obs )
        return

    def setupFromFile( self, tfile, unf="bbb" ):
        standardAnalysis= Analysis( "data mt stand none none none py " + unf )
        analysisVariations= {
            "tc": Analysis( "data tc stand none none none py " + unf ),
            "costt07": Analysis( "data mt costt07 none none none py " + unf ),
            "hw": Analysis( "data mt stand none none none hw " + unf ),
            "sprold": Analysis( "data mt sprold none none none py " + unf  ) }
        self.setupStandardAnalysis( standardAnalysis, tfile )
        self.subtractVariations( analysisVariations, tfile )
        self.calcSystSumSq( analysisVariations.keys() )
        self.readRawEvents( standardAnalysis, analysisVariations, tfile )        
        return
    
    def clone( self, values, sterrs, variationsDelta ):
        aocloned= LEP15AnalysisObservable( self.obs )
        aocloned.aostand= self.aostand
        aocloned.points= self.points
        aocloned.values= values
        aocloned.sterrs= sterrs
        aocloned.variationsDelta= variationsDelta
        aocloned.calcSystSumSq( variationsDelta.keys() )
        return aocloned
 
# LEP2 Analysis
class LEP2AnalysisObservable( AnalysisObservable ):

    def __init__( self, obs ):
        AnalysisObservable.__init__( self, obs )
        return

    def setupFromFile( self, tfile, unf="bbb" ):
        standardAnalysis= Analysis( "data mt stand none none llqq:qqqq:eeqq py " + unf )
        self.setupStandardAnalysis( standardAnalysis, tfile )
        analysisVariations= {
            "tc": Analysis( "data tc stand none none llqq:qqqq:eeqq py " + unf  ),
            "costt07": Analysis( "data mt costt07 none none llqq:qqqq:eeqq py " + unf  ), 
            "sprold": Analysis( "data mt sprold none none llqq:qqqq:eeqq py " + unf  ),
            "hw": Analysis( "data mt stand none none llqq:qqqq:eeqq hw " + unf  ),
            "wqqlnhi": Analysis( "data mt wqqlnhi none none llqq:qqqq:eeqq py " + unf  ),
            "wqqlnlo": Analysis( "data mt wqqlnlo none none llqq:qqqq:eeqq py " + unf  ),
            "wqqqqhi": Analysis( "data mt wqqqqhi none none llqq:qqqq:eeqq py " + unf  ),
            "wqqqqlo": Analysis( "data mt wqqqqlo none none llqq:qqqq:eeqq py " + unf  ),
            "bkghi": Analysis( "data mt stand none none llqq:qqqq:eeqq:hi py " + unf  ), 
            "bkglo": Analysis( "data mt stand none none llqq:qqqq:eeqq:lo py " + unf  ) }
        self.subtractVariations( analysisVariations, tfile )
        self.calcSyst()
        self.readRawEvents( standardAnalysis, analysisVariations, tfile,
                                [ "llqq", "qqqq", "eeqq" ] )
        return

    def calcSyst( self ):
        self.calcSystSumSq( [ "tc", "costt07", "hw", "sprold" ] )
        syerrbkg= self.maxAbsErrorSq( "wqqlnhi", "wqqlnlo" )
        syerrbkg+= self.maxAbsErrorSq( "wqqqqhi", "wqqqqlo" )
        syerrbkg+= self.maxAbsErrorSq( "bkghi", "bkglo" )
        self.syerrs= np.sqrt( np.square( self.syerrs ) + syerrbkg )
        return

    def clone( self, values, sterrs, variationsDelta ):
        aocloned= LEP2AnalysisObservable( self.obs )
        aocloned.aostand= self.aostand
        aocloned.points= self.points
        aocloned.values= values
        aocloned.sterrs= sterrs
        aocloned.variationsDelta= variationsDelta
        aocloned.calcSyst()
        return aocloned
   
# Factory method to create AnalysisObservable objects:
def createAnalysisObservable( tfile, obs="thrust", unf="bbb" ):
    filename= tfile.GetName()
    ao= None
    print( "createAnalysisObservable: create for", obs, "from", filename, end=" " )
    if "sjm91" in filename:
        print( "LEP1AnalysisObservable" )
        ao= LEP1AnalysisObservable( obs )
    elif( "sjm130" in filename or "sjm136" in filename ):
        print( "LEP15AnalysisObservable" )
        ao= LEP15AnalysisObservable( obs )
    elif( "sjm161" in filename or "sjm172" in filename or "sjm183" in filename or
          "sjm189" in filename or "sjm192" in filename or "sjm196" in filename or
          "sjm200" in filename or "sjm202" in filename or "sjm205" in filename or
          "sjm207" in filename ):
        print( "LEP2AnalysisObservable" )
        ao= LEP2AnalysisObservable( obs )
    else:
        print( "no matching AnalysisObservable" )
    ao.setupFromFile( tfile, unf )
    return ao

# Error weighted average of results of input observables:
def combineAnalysisObservables( aobs ):
    firstao= aobs[0]
    for ao in aobs:
        if ao.obs != firstao.obs:
            raise ValueError( "Observables don't match: "+firstao.obs+" "+ao.obs )
    wgts= dict()
    nvalues= len(firstao.values)
    sumwgts= array( "d", nvalues*[ 0.0 ] )
    for ao in aobs:
        wgts[ao]= np.divide( 1.0, np.square( ao.sterrs ) )
        sumwgts= np.add( wgts[ao], sumwgts )
    values= array( "d", nvalues*[ 0.0 ] )
    for ao in aobs:
        values= np.add( np.multiply( ao.values, wgts[ao] ), values )
    values= np.divide( values, sumwgts )
    sterrs= np.divide( 1.0, np.sqrt( sumwgts ) )
    variationsDelta= dict()
    for key in firstao.variationsDelta.keys():
        deltas= array( "d", nvalues*[ 0.0 ] )
        for ao in aobs:
            deltas= np.add( np.multiply( ao.variationsDelta[key], wgts[ao] ), deltas )
        variationsDelta[key]= np.divide( deltas, sumwgts )
    aocombined= firstao.clone( values, sterrs, variationsDelta )        
    return aocombined

# Create combined observable from file list:
def createCombineAnalysisObservables( filenames, obs="thrust" ):
    if len(filenames) == 1:
        f= TFile( filenames[0] )
        aocomb= createAnalysisObservable( f, obs )
    else:
        print( "createCombineAnalysisObservables: combine from", end="" )
        aobs= list()
        for filename in filenames:
            print( filename, end="" )
        print()
        for filename in filenames:
            f= TFile( filename )
            ao= createAnalysisObservable( f, obs )
            aobs.append( ao )
        aocomb= combineAnalysisObservables( aobs )
    return aocomb

# Extract ecm from file names, Match digits between "sjm" and "_" or "."
def ecmFromFilename( filename ):
    import re
    result= re.match( "^sjm(\\d+)[_|.]", filename )
    ecm= result.group( 1 )
    return ecm

# Plot all groomed observables at combined ecms into pdf:
def plotAllGroomedAveraged():
    canv= TCanvas( "canv", "All groomed shapes", 1200, 800 )
    canv.Divide( 3, 2 )
    observables= [ "grthrust" , "grcpar" ]
    filenameslists= [ [ "sjm91_all.root" ],
     [ "sjm130.root", "sjm136.root" ],
     [ "sjm161.root", "sjm172.root", "sjm183.root", "sjm189.root" ],
     [  "sjm192.root", "sjm196.root", "sjm200.root","sjm202.root", "sjm205.root", "sjm207.root" ] ]
    ecms= [ "91", "133", "177", "197" ]
    for obs in observables:
        iecm= 0
        for filenames in filenameslists:
            postfix=""
            if filenames == filenameslists[0] and obs == observables[0]:
                postfix= "("
            elif filenames == filenameslists[-1] and obs == observables[-1]:
                postfix= ")"
            ecm= ecms[iecm]
            plotGroomed( obs, filenames, ecm, logy=1, canv=canv )
            title= "Title: "+obs+" "+ecm+" GeV"
            print( title )
            canv.Print( "plots_averaged.pdf"+postfix, title )
            iecm= iecm+1
    return
   
# Plot all groomed observables into pdf:
def plotAllGroomed():
    filenames= [ "sjm91_all.root",
                     "sjm130.root",
                     "sjm136.root",
                     "sjm161.root",
                     "sjm172.root",
                     "sjm183.root",
                     "sjm189.root",
                     "sjm192.root",
                     "sjm196.root",
                     "sjm200.root",
                     "sjm202.root",
                     "sjm205.root",
                     "sjm207.root" ]
    canv= TCanvas( "canv", "All groomed shapes", 1200, 800 )
    canv.Divide( 3, 2 )
    observables= [ "grthrust" , "grcpar" ]
    for obs in observables:
        for filename in filenames:
            postfix=""
            if filename == filenames[0] and obs == observables[0]:
                postfix= "("
            elif filename == filenames[-1] and obs == observables[-1]:
                postfix= ")"
            ecm= ecmFromFilename( filename )
            plotGroomed( obs, [ filename ], ecm, logy=1, canv=canv )
            title= "Title: "+obs+" "+ecm+" GeV"
            print( title )
            canv.Print( "plots.pdf"+postfix, title )    
    return

# Plot groomed observables: 
def plotGroomed( obs="grthrust", filenames=[ "sjm136_test.root" ], ecm="136",
                 logy=1, canv=None ):
    thplotoptions= { "xmin": 0.0, "xmax": 0.5, "ymin": 0.005, "ymax": 50.0,
                     "markerStyle": 20, "markerSize": 0.5,
                     "title": "groomed Thrust",
                     "xlabel": "1-T_{gr}", "ylabel": "1/\\sigma d\\sigma/d(1-T_{gr})",
                     "logy":logy }
    cpplotoptions= { "xmin": 0.0, "xmax": 1.0, "ymin": 0.03, "ymax": 30.0,
                     "markerStyle": 20, "markerSize": 0.5,
                     "title": "groomed C-parameter",
                     "xlabel": "C_{gr}", "ylabel": '1/\\sigma d\\sigma/d(C_{gr})',
                     "logy":logy }
    plotopts= { "grthrust": thplotoptions, "grcpar": cpplotoptions }
    if canv == None:
        canv= TCanvas( "canv", obs+" "+ecm, 1200, 800 )
    icanv= 0
    for beta in [ "0.0", "1.0" ]:
        for zcut in [ "0.05", "0.10", "0.15" ]:
            icanv= icanv+1
            canv.cd( icanv )
            gPad.SetLeftMargin( 0.15 )
            gPad.SetRightMargin( 0.025 )
            key= obs + "_" + beta + "_" + zcut
            print( key )
            aogr= createCombineAnalysisObservables( filenames, key )
            aogr.plot( plotopts[obs] )
            tl= TLegend( 0.4, 0.8, 0.85, 0.85 )
            tl.SetTextSize( 0.05 )
            tl.SetBorderSize( 0 )
            tl.AddEntry( key, "OPAL "+ecm+" GeV", "ep" )
            tl.Draw( "same" )
            txt= TLatex( 0.6, 0.7, "#beta="+beta+ " z_{cut}="+zcut )
            txt.SetNDC( True )
            txt.SetTextSize( 0.035 )
            txt.Draw()            
    return

# Check jet rates add up to one:
def checkJetrates( filename="sjm91_all.root", obs="durhamycut" ):
    f= TFile( filename )
    valuesmap= dict()
    for rate in [ "R2", "R3", "R4", "R5", "R6" ]:
        ao= createAnalysisObservable( f, obs+rate )
        valuesmap[rate]= ao.values
    valuessum= valuesmap["R2"]
    for rate in [ "R3", "R4", "R5", "R6" ]:
        valuessum= np.add( valuessum, valuesmap[rate] )
    print( valuessum )
    return

# Compare y23 to M.T. Ford:
def compareY23ds():
    canv= TCanvas( "canv", "y_{23}(D) comparison 91 - 189", 1000, 1200 )
    canv.Divide( 2, 3 )
    canv.cd( 1 )
    compareY23d( "sjm91_all.root" )
    canv.cd( 2 )
    compareY23d( "sjm133.root" )
    canv.cd( 3 )
    compareY23d( "sjm161.root" )
    canv.cd( 4 )
    compareY23d( "sjm172.root" )
    canv.cd( 5 )
    compareY23d( "sjm183.root" )
    canv.cd( 6 )
    compareY23d( "sjm189.root" )
    canv2= TCanvas( "canv2", "y_{23}(D) comparison 192 - 207", 1000, 1200 )
    canv2.Divide( 2, 3 )
    canv2.cd( 1 )
    compareY23d( "sjm192.root" )
    canv2.cd( 2 )
    compareY23d( "sjm196.root" )
    canv2.cd( 3 )
    compareY23d( "sjm200.root" )
    canv2.cd( 4 )
    compareY23d( "sjm202.root" )
    canv2.cd( 5 )
    compareY23d( "sjm205.root" )
    canv2.cd( 6 )
    compareY23d( "sjm207.root" )
    return
def compareY23d( filename="sjm91_all.root", mtffilename=None, opt="m" ):
    if mtffilename == None:
        ecm= ecmFromFilename( filename )
        mtffilename= "mtford-y23d"+ecm+".txt" 
    arrays= ascii2arrays( mtffilename )
    mtfordpointsl= arrays[0]
    mtfordpointsr= arrays[1]
    mtfordpoints= np.divide( np.add( arrays[0], arrays[1] ), 2.0 )
    mtfordvalues= arrays[2]
    mtfordsterrs= arrays[3]
    mtfordsyerrs= arrays[4]
    mtforderrs= np.sqrt( np.add( np.square( mtfordsterrs ),  np.square( mtfordsyerrs ) ) )
    if filename=="sjm133.root":
        f1= TFile( "sjm130.root" )
        ao1= createAnalysisObservable( f1, "durhamymerge23" )
        f2= TFile( "sjm136.root" )
        ao2= createAnalysisObservable( f2, "durhamymerge23" )
        ao= combineAnalysisObservables( [ ao1, ao2 ] )
    else:
        f= TFile( filename )
        ao= createAnalysisObservable( f, "durhamymerge23" )
    npoints= len( mtfordpoints )
    vex= array( "d", npoints*[0.0] )
    tgest= TGraphErrors( npoints, mtfordpoints, mtfordvalues, vex, mtfordsterrs )
    tgetot= TGraphErrors( npoints, mtfordpoints, mtfordvalues, vex, mtforderrs )
    plotoptions= { "xmin": 0.0003, "xmax": 0.5, "ymin": 0.5, "ymax": 500.0,
                   "markerStyle": 20, "markerSize": 0.75,
                   "title": "Durham y23 "+filename,
                   "xlabel": "y_{23}", "ylabel": "1/\\sigma d\\sigma/dy_{23}",
                   "logx":1, "logy":1 }
    ao.plot( plotoptions, opt )
    tgetot.SetMarkerStyle( 24 )
    tgetot.SetMarkerSize( 1.25 )
    tgetot.SetName( "mtford" )
    tgetot.Draw( "psame" )
    tgest.Draw( "psame" )
    tl= TLegend( 0.7, 0.9, 0.7, 0.9 )
    tl.AddEntry( "mtford", "M.T. Ford thesis", "ep" )
    tl.AddEntry( "durhamymerge23", "sjmanalysis", "ep" )
    tl.Draw()
    return

# Compare thrust to M.T. Ford:
def compareThrusts():
    canv= TCanvas( "canv", "Thrust comparison to M.T. Ford", 1000, 1200 )
    canv.Divide( 2, 3 )
    canv.cd( 1 )
    compareThrust( "sjm91_all.root" )
    canv.cd( 2 )
    compareThrust( "sjm133.root" )
    canv.cd( 3 )
    compareThrust( "sjm161.root" )
    canv.cd( 4 )
    compareThrust( "sjm172.root" )
    canv.cd( 5 )
    compareThrust( "sjm183.root" )
    canv.cd( 6 )
    compareThrust( "sjm189.root" )
    canv.Print( "thrustplots.pdf(", "Title: 91 - 189 GeV" ) 
    canv.cd( 1 )
    compareThrust( "sjm192.root" )
    canv.cd( 2 )
    compareThrust( "sjm196.root" )
    canv.cd( 3 )
    compareThrust( "sjm200.root" )
    canv.cd( 4 )
    compareThrust( "sjm202.root" )
    canv.cd( 5 )
    compareThrust( "sjm205.root" )
    canv.cd( 6 )
    compareThrust( "sjm207.root" )
    canv.Print( "thrustplots.pdf)", "Title: 192 - 207 GeV" ) 
    return
def compareThrust( filename="sjm91_all.root", mtffilename=None ):
    if mtffilename == None:
        ecm= ecmFromFilename( filename )
        mtffilename= "mtford-thrust"+ecm+".txt" 
    arrays= ascii2arrays( mtffilename )
    mtfordvalues= arrays[2]
    mtfordsterrs= arrays[3]
    mtfordsyerrs= arrays[4]
    mtforderrs= np.sqrt( np.add( np.square( mtfordsterrs ),  np.square( mtfordsyerrs ) ) )
    if filename=="sjm133.root":
        aothrust= createCombineAnalysisObservables( ( "sjm130.root", "sjm136.root" ),
                                                    "lepthrust" )
    else:
        f= TFile( filename )
        aothrust= createAnalysisObservable( f, "lepthrust" )    
    vx= array( "d", aothrust.aostand.getPointsCenter() )
    npoints= len(vx)-1
    vex= array( "d", npoints*[0.0] )
    tgethrustst= TGraphErrors( npoints, vx, mtfordvalues, vex, mtfordsterrs )
    tgethrusttot= TGraphErrors( npoints, vx, mtfordvalues, vex, mtforderrs )
    plotoptions= { "xmin": 0.0, "xmax": 0.5, "ymin": 0.2, "ymax": 30,
                   "markerStyle": 20, "markerSize": 0.8,
                   "title": "Thrust "+filename,
                   "xlabel": "1-T", "ylabel": "1/\\sigma d\\sigma/d(1-T)",
                   "logy": 1 }
    aothrust.plot( plotoptions )
    tgethrusttot.SetMarkerStyle( 24 )
    tgethrusttot.SetMarkerSize( 1.25 )
    tgethrusttot.SetName( "mtford" )
    tgethrusttot.Draw( "psame" )
    tgethrustst.Draw( "psame" )
    tl= TLegend( 0.6, 0.75, 0.85, 0.9 )
    tl.AddEntry( "mtford", "M.T. Ford thesis", "ep" )
    tl.AddEntry( "thrust", "sjmanalysis", "ep" )
    tl.Draw()
    return

# Compare PCONE OPAL results for given variant and jetrate:
def comparePxcone( filename="sjm91_all.root", optKind="emin", optRate="R2" ):

    pr097vals= dict()
    pr097st= dict()
    pr097sy= dict()
    arrays= ascii2arrays( "pr097-pxcone"+optKind+".txt" )
    pr097pts= arrays[0]
    pr097vals["R2"]= arrays[1]
    pr097st["R2"]= arrays[2]
    pr097sy["R2"]= arrays[3]
    pr097vals["R3"]= arrays[4]
    pr097st["R3"]= arrays[5]
    pr097sy["R3"]= arrays[6]
    pr097vals["R4"]= arrays[7]
    pr097st["R4"]= arrays[8]
    pr097sy["R4"]= arrays[9]        
    npr097pts= len( pr097pts )
    vexpr097= array( "d", npr097pts*[0.0] )
    pr097vals= np.divide( pr097vals[optRate], 100.0 )
    pr097st= np.divide( pr097st[optRate], 100.0 )
    pr097sy= np.divide( pr097sy[optRate], 100.0 )    
    pr097tot= np.sqrt( np.add( np.square( pr097st ),  np.square( pr097sy ) ) )

    pr408vals= dict()
    pr408st= dict()
    pr408sy= dict()
    arrays= ascii2arrays( "pr408-pxcone"+optKind+"91.txt" )
    pr408pts= arrays[0]
    pr408vals["R2"]= arrays[1]
    pr408st["R2"]= arrays[2]
    pr408sy["R2"]= arrays[3]
    pr408vals["R3"]= arrays[4]
    pr408st["R3"]= arrays[5]
    pr408sy["R3"]= arrays[6]
    pr408vals["R4"]= arrays[7]
    pr408st["R4"]= arrays[8]
    pr408sy["R4"]= arrays[9]        
    npr408pts= len( pr408pts )
    vexpr408= array( "d", npr408pts*[0.0] )
    pr408vals= pr408vals[optRate]
    pr408st= np.divide( pr408st[optRate], 100.0 )
    pr408sy= np.divide( pr408sy[optRate] , 100.0 )
    pr408tot= np.sqrt( np.add( np.square( pr408st ),  np.square( pr408sy ) ) )
    
    f= TFile( filename )
    aopxcone= createAnalysisObservable( f, "pxcone"+optKind+optRate )
    xmax= { "R": 1.7, "emin": 27.0 }
    ymax= { "R2": 1.1, "R3": 0.35, "R4": 0.18 }
    xlabel= { "R": "R [rad.]", "emin": "E_{min} [GeV]" }
    ylabel= { "R2": "2-jet rate", "R3": "3-jet rate", "R4": "4-jet rate" }
    plotoptions= { "xmin": 0.0, "xmax": xmax[optKind], "ymin": 0.0, "ymax": ymax[optRate],
                       "markerStyle": 20, "markerSize": 0.8,
                       "xlabel": xlabel[optKind], "ylabel": ylabel[optRate],
                       "title": "Cone "+optKind+" "+filename }
    aopxcone.plot( plotoptions )
    xshift= { "R": 0.02, "emin": 0.2 }
    pr097pts= np.add( pr097pts, -xshift[optKind] )
    tgepr097= TGraphErrors( npr097pts, pr097pts, pr097vals, vexpr097, pr097tot )
    tgepr097.SetMarkerStyle( 24 )
    tgepr097.SetMarkerSize( 1.0 )
    tgepr097.SetName( "pr097" )
    tgepr097.Draw( "psame" )
    pr408pts= np.add( pr408pts, xshift[optKind] )
    tgepr408= TGraphErrors( npr408pts, pr408pts, pr408vals, vexpr408, pr408tot )
    tgepr408.SetMarkerStyle( 29 )
    tgepr408.SetMarkerSize( 1.0 )
    tgepr408.SetName( "pr408" )
    tgepr408.Draw( "psame" )    
    tl= TLegend( 0.7, 0.5, 0.9, 0.7 )
    tl.AddEntry( "pr097", "OPAL PR097", "ep" )
    tl.AddEntry( "pr408", "OPAL PR408", "ep" )
    tl.AddEntry( "pxcone"+optKind+optRate, filename, "ep" )
    tl.Draw()
    return

# Compare OPAL PXCONE results:
def comparePxcones( filename="sjm91_all.root" ):
    canv= TCanvas( "canv", "PXCONE comparison", 1000, 1200 )
    canv.Divide(2,3)
    canv.cd(1)
    comparePxcone( filename, "R", "R2" )
    canv.cd(2)
    comparePxcone( filename, "emin", "R2" )
    canv.cd(3)
    comparePxcone( filename, "R", "R3" )
    canv.cd(4)
    comparePxcone( filename, "emin", "R3" )
    canv.cd(5)
    comparePxcone( filename, "R", "R4" )
    canv.cd(6)
    comparePxcone( filename, "emin", "R4" )
    canv.SaveAs( "comparePxcones.pdf" )
    return

# Compare antikt, siscone and PXCONE jets in same plot        
def compareConejets( filename="sjm91_all.root", optKind="R", optR="R3" ):
    f= TFile( filename )
    canv= TCanvas( "canv", "Cone jets", 800, 800 )
    algantikt= "antikt"+optKind
    algsiscone= "siscone"+optKind
    # PXCONE with "2" on same points as other cone jet algs
    algpxcone= "pxcone"+optKind+"2"
    aktao= createAnalysisObservable( f, algantikt+optR )
    ymax= { "R2": 1.0, "R3": 0.5, "R4": 0.3, "R5": 0.3, "R6": 0.3 }
    xmax= { "R": 1.0, "emin": 0.15 }
    xlabel= { "R": "R [rad]", "emin": "Emin fraction" }
    plotoptions= { "xmin": 0.0, "xmax": xmax[optKind], "ymin": 0.0, "ymax": ymax[optR],
                   "markerStyle": 20, "markerSize": 1.0,
                   "title": "Cone "+optKind+" "+optR+" "+filename,
                   "xlabel": xlabel[optKind], "ylabel": optR }
    akttgest, akttgesy= aktao.plot( plotoptions )
    sisao= createAnalysisObservable( f, algsiscone+optR )
    plotoptions["markerStyle"]= 21
    plotoptions["xshift"]= xmax[optKind]/100.0
    sistgest, sistgesy= sisao.plot( plotoptions, "s" )
    pxao= createAnalysisObservable( f, algpxcone+optR )
    plotoptions["markerStyle"]= 22
    plotoptions["xshift"]= -xmax[optKind]/100.0
    pxtgest, pxtgesy= pxao.plot( plotoptions, "s" )
    l= TLegend( 0.7, 0.7, 0.9, 0.9 )
    l.AddEntry( algantikt+optR, "anti-k_t "+optR, "ep" )
    l.AddEntry( algsiscone+optR, "SISCone "+optR, "ep" )
    l.AddEntry( algpxcone+optR, "PXCONE "+optR, "ep" )
    l.Draw()
    return

# Compare Andrii's Durham jet rates 
def compareAllDurhamjetrates():
    canv= TCanvas( "canv", "Durham jetrates comparison", 1000, 1200 )
    canv.Divide(2,3)
    canv.cd( 1 )
    compareDurhamjetrates( "sjm91_all.root",
                               "/home/iwsatlas1/skluth/Downloads/JRTMC/share/NEW/data.dat",
                               "donkers-durhamjets91.txt" )
    canv.cd( 2 )
    compareDurhamjetrates( "sjm130.root",
                               "/home/iwsatlas1/skluth/Downloads/JRTMC/share/NEW2/data.dat",
                               None )
    canv.cd( 3 )
    compareDurhamjetrates( "sjm136.root",
                               "/home/iwsatlas1/skluth/Downloads/JRTMC/share/NEW3/data.dat",
                               None )
    canv.cd( 4 )
    compareDurhamjetrates( "sjm161.root",
                               "/home/iwsatlas1/skluth/Downloads/JRTMC/share/NEW4/data.dat",
                               "donkers-durhamjets161.txt" )
    canv.cd( 5 )
    compareDurhamjetrates( "sjm189.root",
                               "/home/iwsatlas1/skluth/Downloads/JRTMC/share/NEW7/data.dat",
                              "donkers-durhamjets189.txt" )
    canv.cd( 6 )
    compareDurhamjetrates( "sjm192.root",
                               "/home/iwsatlas1/skluth/Downloads/JRTMC/share/NEW8/data.dat",
                               "donkers-durhamjets192.txt" )
    return
def compareDurhamjetrates( filename="sjm91_all.root",
                           datafilename="/home/iwsatlas1/skluth/Downloads/JRTMC/share/NEW/data.dat",
                           donkersfilename="donkers-durhamjets91.txt" ):
    f= TFile( filename )
    R2ao= createAnalysisObservable( f, "durhamycutfjR2" )
    R3ao= createAnalysisObservable( f, "durhamycutfjR3" )
    plotoptions= { "xmin": 0.0005, "xmax": 0.5, "ymin": 0.0, "ymax": 1.05, "markerStyle": 20,
                       "markerSize": 0.75, "title": "Durham R2 and R3 "+filename,
                       "xlabel": "y_{cut}", "ylabel": "Jet rates", "logx": 1 }
    R2tgest, R2tgesy= R2ao.plot( plotoptions )
    plotoptions["markerStyle"]= 21
    R3tgest, R3tgesy= R3ao.plot( plotoptions, "s" )

    arrays= ascii2arrays( datafilename )
    ycutpoints= arrays[0]
    R2values= np.divide( arrays[1], 100.0 )
    R2sterrs= np.divide( arrays[2], 100.0 )
    R2syerrs= np.divide( arrays[3], 100.0 )
    R3values= np.divide( arrays[4], 100.0 )
    R3sterrs= np.divide( arrays[5], 100.0 )
    R3syerrs= np.divide( arrays[6], 100.0 )
    R2errs= np.sqrt( np.add( np.square( R2sterrs ), np.square( R2syerrs ) ) )
    R3errs= np.sqrt( np.add( np.square( R3sterrs ), np.square( R3syerrs ) ) )
    n= len(ycutpoints)
    xerrs= array( "d", n*[0.0] )

    R2datatge= TGraphErrors( n, ycutpoints, R2values, xerrs, R2errs )
    R2datatge.SetMarkerStyle( 24 )
    R2datatge.SetMarkerSize( 0.75 )    
    R2datatge.SetName( "R2datatge" )
    R2datatge.Draw( "psame" )
    R3datatge= TGraphErrors( n, ycutpoints, R3values, xerrs, R3errs )
    R3datatge.SetMarkerStyle( 25 )
    R3datatge.SetMarkerSize( 0.75 )    
    R3datatge.SetName( "R3datatge" )
    R3datatge.Draw( "psame" )

    legend= TLegend( 0.6, 0.6, 0.9, 0.9 )
    R2tgesy.SetName( "R2tgesy" )
    legend.AddEntry( "R2tgesy", "OPAL R2", "pe" )
    R3tgesy.SetName( "R3tgesy" )
    legend.AddEntry( "R3tgesy", "OPAL R3", "pe" )
    legend.AddEntry( "R2datatge", "Andrii R2", "pe" )
    legend.AddEntry( "R3datatge", "Andrii R3", "pe" )

    if donkersfilename:
        dkarrays= ascii2arrays( donkersfilename )
        dkycutpoints= np.power( 10.0, dkarrays[0] )
        dkR2values= dkarrays[1]
        dkR2sterrs= np.divide( dkarrays[2], 100.0 )
        dkR2syerrs= np.divide( dkarrays[3], 100.0 )
        dkR3values= dkarrays[4]
        dkR3sterrs= np.divide( dkarrays[5], 100.0 )
        dkR3syerrs= np.divide( dkarrays[6], 100.0 )
        dkR2errs= np.sqrt( np.add( np.square( dkR2sterrs ), np.square( dkR2syerrs ) ) )
        dkR3errs= np.sqrt( np.add( np.square( dkR3sterrs ), np.square( dkR3syerrs ) ) )
        dkn= len( dkycutpoints )
        dkxerrs= array( "d", dkn*[0.0] )
        dkR2datatge= TGraphErrors( dkn, dkycutpoints, dkR2values, dkxerrs, dkR2errs )
        dkR2datatge.SetMarkerStyle( 26 )
        dkR2datatge.SetMarkerSize( 0.75 )    
        dkR2datatge.SetName( "dkR2datatge" )
        dkR2datatge.Draw( "psame" )
        dkR3datatge= TGraphErrors( dkn, dkycutpoints, dkR3values, dkxerrs, dkR3errs )
        dkR3datatge.SetMarkerStyle( 27 )
        dkR3datatge.SetMarkerSize( 0.75 )    
        dkR3datatge.SetName( "dkR3datatge" );
        dkR3datatge.Draw( "psame" )
        legend.AddEntry( "dkR2datatge", "Donkers R2", "pe" )
        legend.AddEntry( "dkR3datatge", "Donkers R3", "pe" )
        
    legend.Draw()
    return

# Compare EEC from various sources with own measurements 
def compareEEC( filename="sjm91_all.root", datafilename="../EECMC/share/OPAL/data.dat" ):
    f= TFile( filename )
    ao= createAnalysisObservable( f, "EEC" )
    tokens= datafilename.split( "/" )
    exp= tokens[3]
    plotoptions= { "xmin": 0.0, "xmax": 3.14159, "ymin": 0.05, "ymax": 5.0,
                   "markerStyle": 20, "markerSize": 0.5, "drawas": "3", "fillcolor": 6,
                   "title": "EEC "+exp,
                   "xlabel": "\\chi\\ [rad.]", "ylabel": "1/\\sigma d\\Sigma/d\\chi",
                   "logy": 1 }
    tgest, tgesy= ao.plot( plotoptions )
    lines= [ line.rstrip( '\n' ) for line in open( datafilename ) ]
    n= len( lines )
    points= TVectorD( n )
    values= TVectorD( n )
    errors= TVectorD( n )
    perrs= TVectorD(n)
    grad2rad= 3.14159/180.0
    for i in range( n ):
        line= (lines[i]).split()
        points[i]= float(line[0])*grad2rad
        values[i]= float(line[3])
        errors[i]= float(line[4])
        perrs[i]= 0.0
    datatge= TGraphErrors( points, values, perrs, errors )
    datatge.SetMarkerStyle( 20 )
    datatge.SetMarkerSize( 0.5 )    
    datatge.Draw( "psame" )
    legend= TLegend( 0.2, 0.7, 0.5, 0.85 )
    datatge.SetName( "datatge" );
    tgesy.SetName( "tgesy" )
    legend.AddEntry( "datatge", exp+" data", "pe" )
    legend.AddEntry( "tgesy", "OPAL "+filename, "f" )
    legend.Draw()
    return 

def compareEECs( filename="sjm91_all.root" ):
    canv= TCanvas( "canv", "EEC comparison", 1000, 1200 )
    canv.Divide(2,3)
    canv.cd(1)
    compareEEC( filename, datafilename="../EECMC/share/OPAL/data.dat" )
    canv.cd(2)
    compareEEC( filename, datafilename="../EECMC/share/OPAL2/data.dat" )
    canv.cd(3)
    compareEEC( filename, datafilename="../EECMC/share/OPAL3/data.dat" )
    canv.cd(4)
    compareEEC( filename, datafilename="../EECMC/share/DELPHI/data.dat" )
    canv.cd(5)
    compareEEC( filename, datafilename="../EECMC/share/SLD/data.dat" )
    canv.cd(6)
    compareEEC( filename, datafilename="../EECMC/share/L3/data.dat" )
    canv.SaveAs( "compareEECs.pdf" )
    return


def testMigrationMatrix( obs="lepthrust", filename="sjm91_96-98.root" ):
    hdetstr= obs+" py mt stand"
    hhstr= obs+" py hadron stand"
    hhnrstr= obs+" py hadron none nonrad"
    mstr= "migr "+obs+" py mt stand hadron"
    f= TFile( filename )
    hdet= f.Get( hdetstr )
    hdet.Print()
    m= f.Get( mstr )
    m.Print()
    hh= f.Get( hhstr )
    hh.Print()
    hhnr= f.Get( hhnrstr )
    hhnr.Print()
    nbin= hdet.GetNbinsX()    
    import numpy as np
    valuesd= np.array( nbin*[0.0] )
    valuesh= np.array( nbin*[0.0] )
    valueshnr= np.array( nbin*[0.0] )
    cacc= np.array( nbin*[0.0] )
    R= np.array( np.zeros( (nbin,nbin) ) )
    for i in range( nbin ):
        valuesd[i]= hdet.GetBinContent( i+1 )*hdet.GetEntries()*hdet.GetBinWidth( i+1 )
        valuesh[i]= hh.GetBinContent( i+1 )*hh.GetEntries()*hh.GetBinWidth( i+1 )
        valueshnr[i]= hhnr.GetBinContent( i+1 )*hhnr.GetEntries()*hhnr.GetBinWidth( i+1 )
        if valuesh[i] != 0.0:
            cacc[i]= valueshnr[i]/valuesh[i]
        else:
            cacc[i]= 0.0
        for j in range( nbin ):
            R[j,i]= m.GetBinContent( i+1, j+1 )

    fmt71= "{:7.1f}"
    fmt73= "{:7.3f}"
    for i in range( nbin ):
        print( fmt71.format( valueshnr[i] ), end=" " )
        print( fmt71.format( valuesh[i] ), end=" " )
        for j in range( nbin ):
            print( fmt73.format( R[i,j] ), end=" " )
        print()
    print( "               ", end=" " )
    for i in range( nbin ):
        print( fmt71.format( valuesd[i] ), end=" " )
    print()

    for i in range( nbin ):
        sumcol= sum( R[:,i] )
        if sumcol != 0.0:
            R[:,i]/= sumcol
    C= np.diag( cacc )
    CR= np.dot( C, R )
    valuesc= np.dot( CR, valuesd )
    print( valueshnr )
    print( valuesc )

    return

