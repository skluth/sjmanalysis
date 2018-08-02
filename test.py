
from ROOT import TFile, gROOT, gPad, TVectorD, TObject

from ROOT import TGraphErrors, TH1D, TLegend, TCanvas
TGraphErrors.__init__._creates= False
TH1D.__init__._creates= False
TLegend.__init__._creates= False
TCanvas.__init__._creates= False

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
        return

    def setupStandardAnalysis( self, standardAnalysis, tfile ):
        self.aostand= getAnalysisObjectFromFile( tfile, self.obs, standardAnalysis )
        self.points= array( "d", self.aostand.getPoints() )
        self.values= array( "d", self.aostand.getValues() )
        self.sterrs= array( "d", self.aostand.getErrors() )
        return
    
    def subtractVariations( self, analysisVariations, tfile ):
        self.variationsDelta= dict()
        for key in analysisVariations.keys():
            variationData= array( "d", getAnalysisObjectFromFile( tfile, self.obs, analysisVariations[key] ).getValues() )
            self.variationsDelta[key]= np.subtract( variationData, self.values )
        return
       
    def calcSystSumSq( self, keys ):
        self.syerrs= 0.0
        for key in keys:
            self.syerrs+= np.square( self.variationsDelta[key] )
        self.syerrs= np.sqrt( self.syerrs )
        return
    
    def printResults( self, width=7, precision=3, opt="" ):
        print "Results for", self.obs
        print self.aostand.getPointLabel(),
        fmt= "{:>"+str(width)+"}"
        for key in [ "val", "stat", "sys" ]:
            print fmt.format( key ),
        if "d" in opt:
            for key in sorted( self.variationsDelta.keys() ):
                print fmt.format( key ),
        print
        if "m" in opt:
            sterrs= self.aostand.getErrors( "m" )
        else:
            sterrs= self.sterrs            
        fmt="{:"+str(width)+"."+str(precision)+"f}"
        for i in range(len(self.values)):
            if( self.obs.find( "EEC" ) >= 0 and not
                self.aostand.getPointStr(i).Index( "Integral" ) >= 0 ):
                rad2grad= 180.0/3.14159
                leftedge= self.points[i]*rad2grad
                rightedge= self.points[i+1]*rad2grad
                print "{0:3.0f} {1:3.0f}  ".format( leftedge, rightedge ),
            else:
                print self.aostand.getPointStr(i),
            print fmt.format( self.values[i] ),
            print fmt.format( sterrs[i] ),
            print fmt.format( self.syerrs[i] ),
            if "d" in opt:
                for key in sorted( self.variationsDelta.keys() ):
                    print fmt.format( self.variationsDelta[key][i] ),
            print
        return

    def printErrors( self, width=7, precision=4 ):
        from math import sqrt
        errorMatrix= self.aostand.getErrorMatrix()
        fmt="{:"+str(width)+"."+str(precision)+"f}"
        for i in range( len( self.sterrs )-1 ):
            binw= self.points[i+1]-self.points[i]
            diagError= sqrt( errorMatrix(i,i) )/binw
            print fmt.format( self.sterrs[i] ), fmt.format( diagError )
        return

    def plot( self, plotoptions, opt="?" ):
        vx= array( "d", self.aostand.getPointsCenter() )
        npoints= len(vx)
        if "xshift" in plotoptions:
            for i in range(npoints):
                vx[i]+= plotoptions["xshift"]
        vex= array( "d", npoints*[0.0] )
        tgest= TGraphErrors( npoints, vx, self.values, vex, self.sterrs )
        tgesy= TGraphErrors( npoints, vx, self.values, vex, self.syerrs )
        tgesy.SetMarkerStyle( plotoptions["markerStyle"] )
        tgesy.SetMarkerSize( plotoptions["markerSize"] )
        drawas= plotoptions["drawas"] if "drawas" in plotoptions else "p"
        tgesy.SetName( self.obs )
        if "fillcolor" in plotoptions:
            tgesy.SetFillColor(plotoptions["fillcolor"])
            tgest.SetFillColor(plotoptions["fillcolor"])            
        if opt.find( "s" ) >= 0:
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

    
# LEP1 Analysis
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
        return

# LEP1.5 Analysis
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
    print "createAnalysisObservable: create for", obs, "from", filename,
    if "sjm91" in filename:
        print "LEP1AnalysisObservable"
        ao= LEP1AnalysisObservable( obs )
    elif( "sjm130" in filename or "sjm136" in filename ):
        print "LEP15AnalysisObservable"
        ao= LEP15AnalysisObservable( obs )
    elif( "sjm161" in filename or "sjm172" in filename or "sjm183" in filename or
          "sjm189" in filename or "sjm192" in filename or "sjm196" in filename or
          "sjm200" in filename or "sjm202" in filename or "sjm205" in filename or
          "sjm207" in filename ):
        print "LEP2AnalysisObservable"
        ao= LEP2AnalysisObservable( obs )
    else:
        print "no matching AnalysisObservable"
    ao.setupFromFile( tfile, unf )
    return ao

# Error weighted average of results of input observables:
def combineAnalysisObservables( *aobs ):
    for ao in aobs:
        if ao.obs != aobs[0].obs:
            raise ValueError( "Observables don't match: "+ao[0].obs+" "+ao.obs )
    wgts= dict()
    sumwgts= array( "d", len(aobs[0].values)*[ 0.0 ] )
    for ao in aobs:
        wgts[ao]= np.divide( 1.0, np.square( ao.sterrs ) )
        sumwgts= np.add( wgts[ao], sumwgts )
    values= array( "d", len(aobs[0].values)*[ 0.0 ] )
    for ao in aobs:
        values= np.add( np.multiply( ao.values, wgts[ao] ), values )
    values= np.divide( values, sumwgts )
    sterrs= np.divide( 1.0, np.sqrt( sumwgts ) )
    variationsDelta= dict()
    for key in aobs[0].variationsDelta.keys():
        deltas= array( "d", len(aobs[0].values)*[ 0.0 ] )
        for ao in aobs:
            deltas= np.add( np.multiply( ao.variationsDelta[key], wgts[ao] ), deltas )
        variationsDelta[key]= np.divide( deltas, sumwgts )
    aocombined= aobs[0].clone( values, sterrs, variationsDelta )        
    return aocombined

# Check jet rates add up to one:
def checkJetrates( filename="sjm91_all_test.root", obs="durhamycut" ):
    f= TFile( filename )
    valuesmap= dict()
    for rate in [ "R2", "R3", "R4", "R5", "R6" ]:
        ao= createAnalysisObservable( f, obs+rate )
        valuesmap[rate]= ao.values
    valuessum= valuesmap["R2"]
    for rate in [ "R3", "R4", "R5", "R6" ]:
        valuessum= np.add( valuessum, valuesmap[rate] )
    print valuessum
    return

def compareThrusts():
    from ROOT import TCanvas
    canv= TCanvas( "canv", "Thrust comparison", 1000, 1200 )
    canv.Divide( 2, 2 )
    canv.cd( 1 )
    compareThrust( "sjm91_all_test.root", "mtford-thrust91.txt" )
    canv.cd( 2 )
    compareThrust( "sjm133.root", "mtford-thrust133.txt" )
    canv.cd( 3 )
    compareThrust( "sjm189.root", "mtford-thrust189.txt" )
    canv.cd( 4 )
    compareThrust( "sjm192.root", "mtford-thrust192.txt" )
    return

def compareThrust( filename="sjm91_all_test.root", mtffilename="mtford-thrust91.txt" ):
    arrays= ascii2arrays( mtffilename )
    mtfordvalues= arrays[2]
    mtfordsterrs= arrays[3]
    mtfordsyerrs= arrays[4]
    mtforderrs= np.sqrt( np.add( np.square( mtfordsterrs ),  np.square( mtfordsyerrs ) ) )
    if filename=="sjm133.root":
        f1= TFile( "sjm130.root" )
        aothrust1= createAnalysisObservable( f1, "thrust" )
        f2= TFile( "sjm136.root" )
        aothrust2= createAnalysisObservable( f2, "thrust" )
        aothrust= combineAnalysisObservables( aothrust1, aothrust2 )
    else:
        f= TFile( filename )
        aothrust= createAnalysisObservable( f, "thrust" )    
    vx= array( "d", aothrust.aostand.getPointsCenter() )
    npoints= len(vx)-1
    vex= array( "d", npoints*[0.0] )
    tgethrustst= TGraphErrors( npoints, vx, mtfordvalues, vex, mtfordsterrs )
    tgethrusttot= TGraphErrors( npoints, vx, mtfordvalues, vex, mtforderrs )
    plotoptions= { "xmin": 0.0, "xmax": 0.5, "ymin": 0.2, "ymax": 30, "markerStyle": 20,
                       "markerSize": 0.8, "title": "Thrust "+filename, "logy": 1,
                       "xlabel": "1-T", "ylabel": "1/\sigma d\sigma/d(1-T)" }
    aothrust.plot( plotoptions )
    tgethrusttot.SetMarkerStyle( 24 )
    tgethrusttot.SetMarkerSize( 1.25 )
    tgethrusttot.SetName( "mtford" )
    tgethrusttot.Draw( "psame" )
    tgethrustst.Draw( "psame" )
    tl= TLegend( 0.7, 0.9, 0.7, 0.9 )
    tl.AddEntry( "mtford", "M.T. Ford thesis", "ep" )
    tl.AddEntry( "thrust", "sjmanalysis", "ep" )
    tl.Draw()
    return

# Compare PCONE OPAL results for given variant and jetrate:
def comparePxcone( filename="sjm91_all_test.root", optKind="emin", optRate="R2" ):

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
def comparePxcones( filename="sjm91_all_test.root" ):
    from ROOT import TCanvas
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
def compareConejets( filename="sjm91_all_test.root", optKind="R", optR="R3" ):
    f= TFile( filename )
    algantikt= "antikt"+optKind
    algsiscone= "siscone"+optKind
    algpxcone= "pxcone"+optKind+"2"
    aktao= createAnalysisObservable( f, algantikt+optR )
    ymax= { "R2":1.0, "R3":0.5, "R4":0.3, "R5":0.3, "R6":0.3 }
    xmax= { "R":1.0, "emin":0.15 }
    plotoptions= { "xmin": 0.0, "xmax": xmax[optKind], "ymin": 0.0, "ymax": ymax[optR],
                       "markerStyle": 20, "markerSize": 0.8,
                       "title": "Cone "+optKind+" "+optR+" "+filename }
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
    from ROOT import TCanvas
    canv= TCanvas( "canv", "Durham jetrates comparison", 1000, 1200 )
    canv.Divide(2,3)
    canv.cd( 1 )
    compareDurhamjetrates( "sjm91_all_test.root",
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

def compareDurhamjetrates( filename="sjm91_all_test.root",
                               datafilename="/home/iwsatlas1/skluth/Downloads/JRTMC/share/NEW/data.dat",
                               donkersfilename="donkers-durhamjets91.txt" ):
    f= TFile( filename )
    R2ao= createAnalysisObservable( f, "durhamycutfjR2" )
    R3ao= createAnalysisObservable( f, "durhamycutfjR3" )
    plotoptions= { "xmin": 0.0, "xmax": 4.0, "ymin": 0.0, "ymax": 1.05, "markerStyle": 20,
                       "markerSize": 0.75, "title": "Durham R2 and R3 "+filename,
                       "xlabel": "-log10(y_{cut})", "ylabel": "Jet rates" }
    R2tgest, R2tgesy= R2ao.plot( plotoptions )
    plotoptions["markerStyle"]= 21
    R3tgest, R3tgesy= R3ao.plot( plotoptions, "s" )

    arrays= ascii2arrays( datafilename )
    ycutpoints= - np.log10( arrays[0] )
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
        dkycutpoints= np.multiply( dkarrays[0], -1.0 )
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
def compareEEC( filename="sjm91_all_test.root", datafilename="../EECMC/share/OPAL/data.dat" ):
    f= TFile( filename )
    ao= createAnalysisObservable( f, "EEC" )
    tokens= datafilename.split( "/" )
    exp= tokens[3]
    plotoptions= { "xmin": 0.0, "xmax": 3.14159, "ymin": 0.05, "ymax": 5.0, "markerStyle": 20, "markerSize": 0.5, "drawas": "3", "fillcolor": 6, "title": "EEC "+exp, "xlabel": "\chi\ [rad.]", "ylabel": "1/\sigma d\Sigma/d\chi", "logy": 1 }
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

def compareEECs( filename="sjm91_all_test.root" ):
    from ROOT import TCanvas
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

def testMigrationMatrix( obs="thrust", filename="sjm91_96_test.root" ):
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

    width, precision= 7, 3
    fmt= "{:"+str(width)+"."+str(precision)+"f}"
    for i in range( nbin ):
        print fmt.format( valueshnr[i] ),
        print fmt.format( valuesh[i] ),
        for j in range( nbin ):
            print fmt.format( R[i,j] ),
        print
    print "               ",
    for i in range( nbin ):
        print fmt.format( valuesd[i] ),
    print

    for i in range( nbin ):
        sumcol= sum( R[:,i] )
        if sumcol != 0.0:
            R[:,i]/= sumcol
    C= np.diag( cacc )
    CR= np.dot( C, R )
    valuesc= np.dot( CR, valuesd )
    print valueshnr
    print valuesc

    return

