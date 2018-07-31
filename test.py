
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


# Factory method to create AnalysisObject instances
def getAnalysisObjectFromFile( tfile, obs, analysis ):
    ao= None
    key= obs+" "+analysis.getTag()+";1"
    obj= tfile.Get( key )
    if not obj:
        # print "AnalysisObject with key", key, "not in file", tfile.GetName()
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

    def __init__( self, obs, tfile, unf="bbb" ):
        AnalysisObservable.__init__( self, obs )
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

    def __init__( self, obs, tfile, unf="bbb" ):
        AnalysisObservable.__init__( self, obs )
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

# LEP2 Analysis
class LEP2AnalysisObservable( AnalysisObservable ):

    def __init__( self, obs, tfile, unf="bbb" ):
        AnalysisObservable.__init__( self, obs )
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
        self.calcSystSumSq( [ "tc", "costt07", "hw", "sprold" ] )
        syerrbkg= self.maxAbsErrorSq( "wqqlnhi", "wqqlnlo" )
        syerrbkg+= self.maxAbsErrorSq( "wqqqqhi", "wqqqqlo" )
        syerrbkg+= self.maxAbsErrorSq( "bkghi", "bkglo" )
        self.syerrs= np.sqrt( np.square( self.syerrs ) + syerrbkg )
        return

   
# Factory method to create AnalysisObservable objects:
def createAnalysisObservable( tfile, obs="thrust", unf="bbb" ):
    filename= tfile.GetName()
    ao= None
    print "createAnalysisObservable: create for", obs, "from", filename,
    if "sjm91" in filename:
        print "LEP1AnalysisObservable"
        ao= LEP1AnalysisObservable( obs, tfile, unf )
    elif( "sjm161" in filename or "sjm172" in filename or "sjm183" in filename or
          "sjm189" in filename or "sjm192" in filename or "sjm196" in filename or
          "sjm200" in filename or "sjm202" in filename or "sjm205" in filename or
          "sjm207" in filename ):
        print "LEP2AnalysisObservable"
        ao= LEP2AnalysisObservable( obs, tfile, unf )
    else:
        print "no matching AnalysisObservable"
    return ao

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

def compareThrust( filename="sjm91_all_test.root" ):
    if "91" in filename:
        mtfordvalues= array( "d", [ 1.273, 12.26, 18.38, 13.86, 9.80, 6.502, 4.133, 2.649,
                                        1.705, 0.913, 0.3704 ] )
        mtfordsterrs= array( "d", [ 0.019, 0.05, 0.07, 0.06, 0.05, 0.027, 0.022, 0.014, 0.012, 0.006, 0.0033 ] )
        mtfordsyerrs= array( "d", [ 0.043, 0.40, 0.28, 0.15, 0.18, 0.086, 0.037, 0.058, 0.075, 0.020, 0.0090 ] )
    elif "189" in filename:
        mtfordvalues= array( "d", [ 9.36, 21.93, 15.35, 9.82, 7.15, 5.38, 2.89, 2.44, 1.36, 0.72, 0.479 ] )
        mtfordsterrs= array( "d", [ 0.55, 0.80, 0.64, 0.54, 0.47, 0.28, 0.24, 0.16, 0.14, 0.08, 0.065 ] )
        mtfordsyerrs= array( "d", [ 0.58, 0.91, 0.89, 0.68, 0.51, 0.28, 0.16, 0.08, 0.14, 0.11, 0.096 ] )
    mtforderrs= np.sqrt( np.add( np.square( mtfordsterrs ),  np.square( mtfordsyerrs ) ) )
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
    pr097pts= dict()
    pr097pts["emin"]= array( "d", [ 3.0, 5.0, 7.0, 9.0, 12.0, 15.0, 18.0, 21.0, 25.0 ] )
    pr097pts["R"]= array( "d", [ 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5 ] )
    npr097pts= len( pr097pts[optKind] )
    vexpr097= array( "d", npr097pts*[0.0] )
    pr097vals= dict()
    pr097vals["R"]= dict()
    pr097vals["emin"]= dict()
    pr097st= dict()
    pr097st["R"]= dict()
    pr097st["emin"]= dict()
    pr097sy= dict()
    pr097sy["R"]= dict()
    pr097sy["emin"]= dict()
    pr097vals["emin"]["R2"]= array( "d", [ 61.22, 69.55, 74.69, 78.79, 83.92, 88.30, 92.34, 95.82, 98.04 ] )
    pr097st["emin"]["R2"]= array( "d", [ 0.14, 0.10, 0.12, 0.09, 0.08, 0.08, 0.08, 0.06, 0.06 ] )
    pr097sy["emin"]["R2"]= array( "d", [ 1.0, 1.0, 0.96, 1.0, 0.91, 1.04, 1.01, 1.12, 1.44 ] )    
    pr097vals["R"]["R2"]= array( "d", [ 62.33, 68.74, 74.69, 81.03, 87.71, 94.87, 99.47 ] )
    pr097st["R"]["R2"]= array( "d", [ 0.10, 0.11, 0.12, 0.07, 0.07, 0.04, 0.02 ] )
    pr097sy["R"]["R2"]= array( "d", [ 1.63, 1.21, 0.96, 0.48, 0.25, 0.18, 0.10 ] )
    pr097vals["emin"]["R3"]= array( "d", [ 28.84, 24.98, 21.90, 19.08, 15.10, 11.29, 7.46, 3.87, 0.61 ] )
    pr097st["emin"]["R3"]= array( "d", [ 0.13, 0.09, 0.12, 0.10, 0.08, 0.08, 0.08, 0.05, 0.02 ] )
    pr097sy["emin"]["R3"]= array( "d", [ 0.81, 0.77, 0.66, 0.70, 0.54, 0.59, 0.47, 0.24, 0.08 ] )
    pr097vals["R"]["R3"]= array( "d", [ 28.36, 25.30, 21.90, 17.57, 11.99, 5.08, 0.54 ] )
    pr097st["R"]["R3"]= array( "d", [ 0.12, 0.12, 0.12, 0.08, 0.08, 0.04, 0.02 ])
    pr097sy["R"]["R3"]= array( "d", [ 0.50, 0.72, 0.66, 0.35, 0.16, 0.05, 0.02 ] )
    pr097vals["emin"]["R4"]= array( "d", [ 8.26, 4.86, 3.15, 1.97, 0.87, 0.27, 0.03, 0.0, 0.0 ] )
    pr097st["emin"]["R4"]= array( "d", [ 0.07, 0.07, 0.05, 0.05, 0.03, 0.01, 0.01, 0.0, 0.0 ] )
    pr097sy["emin"]["R4"]= array( "d", [ 0.46, 0.28, 0.21, 0.14, 0.08, 0.04, 0.01, 0.0, 0.0 ] )
    pr097vals["R"]["R4"]= array( "d", [ 7.94, 5.19, 3.15, 1.35, 0.29, 0.03, 0.0  ] )
    pr097st["R"]["R4"]= array( "d", [ 0.08, 0.07, 0.05, 0.03, 0.02, 0.01, 0.0 ])
    pr097sy["R"]["R4"]= array( "d", [ 0.42, 0.38, 0.21, 0.12, 0.02, 0.01, 0.0 ] )
    pr097vals= np.divide( pr097vals[optKind][optRate], 100.0 )
    pr097st= np.divide( pr097st[optKind][optRate], 100.0 )
    pr097sy= np.divide( pr097sy[optKind][optRate], 100.0 )
    pr097tot= np.sqrt( np.add( np.square( pr097st ),  np.square( pr097sy ) ) )
    pr408pts= dict()
    pr408pts["emin"]= array( "d", [ 2.0, 6.0, 10.0, 14.0, 18.0, 22.0, 25.50 ] )
    pr408pts["R"]= array( "d", [ 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5 ] )
    npr408pts= len( pr408pts[optKind] )
    vexpr408= array( "d", npr408pts*[0.0] )
    pr408vals= dict()
    pr408vals["R"]= dict()
    pr408vals["emin"]= dict()
    pr408st= dict()
    pr408st["R"]= dict()
    pr408st["emin"]= dict()
    pr408sy= dict()
    pr408sy["R"]= dict()
    pr408sy["emin"]= dict()
    pr408vals["emin"]["R2"]= array( "d", [ 0.5201, 0.7175, 0.8007, 0.8659, 0.9214, 0.9701, 0.9957 ] )
    pr408st["emin"]["R2"]= array( "d", [ 0.0018, 0.0017, 0.0015, 0.0013, 0.0010, 0.0006, 0.0003 ] )
    pr408sy["emin"]["R2"]= array( "d", [ 1.0402, 1.4350, 1.6015, 1.7318, 1.8428, 1.9402, 1.9914 ] )    
    pr408vals["R"]["R2"]= array( "d", [ 0.6142, 0.6794, 0.7422, 0.8049, 0.8740, 0.9474, 0.9944 ] )
    pr408st["R"]["R2"]= array( "d",   [ 0.0018, 0.0017, 0.0016, 0.0015, 0.0012, 0.0008, 0.0003 ] )
    pr408sy["R"]["R2"]= array( "d",   [ 1.2284, 1.3587, 1.4845, 1.6097, 1.7479, 1.8947, 1.9888 ] )
    pr408vals["emin"]["R3"]= array( "d", [ 0.3235, 0.2396, 0.1827, 0.1294, 0.0778, 0.0295, 0.0044 ] )
    pr408st["emin"]["R3"]= array( "d",   [ 0.0018, 0.0016, 0.0015, 0.0013, 0.0010, 0.0006, 0.0002 ] )
    pr408sy["emin"]["R3"]= array( "d",   [ 0.6470, 0.4793, 0.3654, 0.2589, 0.1556, 0.0591, 0.0088 ] )
    pr408vals["R"]["R3"]= array( "d",    [ 0.2870, 0.2604, 0.2237, 0.1818, 0.1233, 0.0522, 0.0055 ] )
    pr408st["R"]["R3"]= array( "d",      [ 0.0017, 0.0017, 0.0016, 0.0015, 0.0013, 0.0009, 0.0003 ] )
    pr408sy["R"]["R3"]= array( "d",      [ 0.5740, 0.5209, 0.4475, 0.3637, 0.2466, 0.1044, 0.0110 ] )
    pr408vals["emin"]["R4"]= array( "d", [ 0.1570, 0.0431, 0.0165, 0.0043, 0.0003, 0.0, 0.0 ] )
    pr408st["emin"]["R4"]= array( "d",   [ 0.0015, 0.0008, 0.0005, 0.0003, 0.0001, 0.0, 0.0 ] )
    pr408sy["emin"]["R4"]= array( "d",   [ 0.3141, 0.0863, 0.0330, 0.0085, 0.0007, 0.0, 0.0 ] )
    pr408vals["R"]["R4"]= array( "d",    [ 0.0987, 0.0601, 0.0342, 0.0134, 0.0028, 0.0004, 0.0 ] )
    pr408st["R"]["R4"]= array( "d",      [ 0.0012, 0.0010, 0.0008, 0.0005, 0.0002, 0.0001, 0.0 ] )
    pr408sy["R"]["R4"]= array( "d",      [ 0.1975, 0.1201, 0.0685, 0.0268, 0.0057, 0.0009, 0.0 ] )
    pr408vals= pr408vals[optKind][optRate]
    pr408st= pr408st[optKind][optRate]
    pr408sy= np.divide( pr408sy[optKind][optRate], 100.0 )
    pr408tot= np.sqrt( np.add( np.square( pr408st ),  np.square( pr408sy ) ) )
    f= TFile( filename )
    aopxcone= createAnalysisObservable( f, "pxcone"+optKind+optRate )
    xmax= { "R": 1.7, "emin": 27.0 }
    ymax= { "R2": 1.1, "R3": 0.35, "R4": 0.18 }
    xlabel= { "R": "R [rad.]", "emin": "E_min [GeV]" }
    ylabel= { "R2": "2-jet rate", "R3": "3-jet rate", "R4": "4-jet rate" }
    plotoptions= { "xmin": 0.0, "xmax": xmax[optKind], "ymin": 0.0, "ymax": ymax[optRate],
                       "markerStyle": 20, "markerSize": 0.8,
                       "xlabel": xlabel[optKind], "ylabel": ylabel[optRate],
                       "title": "Cone "+optKind+" "+filename }
    aopxcone.plot( plotoptions )
    xshift= { "R": 0.02, "emin": 0.2 }
    pr097pts= np.add( pr097pts[optKind], -xshift[optKind] )
    tgepr097= TGraphErrors( npr097pts, pr097pts, pr097vals, vexpr097, pr097tot )
    tgepr097.SetMarkerStyle( 24 )
    tgepr097.SetMarkerSize( 1.0 )
    tgepr097.SetName( "pr097" )
    tgepr097.Draw( "psame" )
    pr408pts= np.add( pr408pts[optKind], xshift[optKind] )
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
def compareDurhamjetrates( filename="sjm91_all_test.root",
                            datafilename="/home/iwsatlas1/skluth/Downloads/JRTMC/share/NEW/data.dat" ):
    f= TFile( filename )
    R2ao= createAnalysisObservable( f, "durhamycutfjR2" )
    R3ao= createAnalysisObservable( f, "durhamycutfjR3" )
    plotoptions= { "xmin": 0.0, "xmax": 4.0, "ymin": 0.0, "ymax": 1.05, "markerStyle": 20,
                       "markerSize": 0.75, "title": "Durham R2 and R3 "+filename,
                       "xlabel": "-log10(y_{cut})", "ylabel": "Jet rates" }
    R2tgest, R2tgesy= R2ao.plot( plotoptions )
    plotoptions["markerStyle"]= 21
    R3tgest, R3tgesy= R3ao.plot( plotoptions, "s" )
    lines= [ line.rstrip( '\n' ) for line in open( datafilename ) ]
    n= len( lines )
    ycutpoints= array( "d", n*[0.0] )
    R2values= array( "d", n*[0.0] )
    R2sterrs= array( "d", n*[0.0] )
    R2syerrs= array( "d", n*[0.0] )
    R3values= array( "d", n*[0.0] )
    R3sterrs= array( "d", n*[0.0] )
    R3syerrs= array( "d", n*[0.0] )
    xerrs= array( "d", n*[0.0] )
    from math import log10
    for i in range( n ):
        line= (lines[i]).split()
        ycutpoints[i]= - log10( float( line[0] ) )
        R2values[i]= float( line[1] )/100.0
        R2sterrs[i]= float( line[2] )/100.0
        R2syerrs[i]= float( line[3] )/100.0
        R3values[i]= float( line[4] )/100.0
        R3sterrs[i]= float( line[5] )/100.0
        R3syerrs[i]= float( line[6] )/100.0
    R2datatge= TGraphErrors( n, ycutpoints, R2values, xerrs, R2syerrs )
    R2datatge.SetMarkerStyle( 24 )
    R2datatge.SetMarkerSize( 0.75 )    
    R2datatge.SetName( "R2datatge" );
    R2datatge.Draw( "psame" )
    R3datatge= TGraphErrors( n, ycutpoints, R3values, xerrs, R3syerrs )
    R3datatge.SetMarkerStyle( 25 )
    R3datatge.SetMarkerSize( 0.75 )    
    R3datatge.SetName( "R3datatge" );
    R3datatge.Draw( "psame" )
    legend= TLegend( 0.6, 0.6, 0.9, 0.9 )
    R2tgesy.SetName( "R2tgesy" )
    legend.AddEntry( "R2tgesy", "OPAL R2", "pe" )
    R3tgesy.SetName( "R3tgesy" )
    legend.AddEntry( "R3tgesy", "OPAL R3", "pe" )
    legend.AddEntry( "R2datatge", "Andrii R2", "pe" )
    legend.AddEntry( "R3datatge", "Andrii R3", "pe" )
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

