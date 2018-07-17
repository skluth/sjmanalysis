
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

# Compare antikt, siscone and PXCONE jets in same plot        
def compareConejets( filename="sjm91_96_test.root", optKind="R", optR="R3" ):
    f= TFile( filename )
    algantikt= "antikt"+optKind
    algsiscone= "siscone"+optKind
    algpxcone= "pxcone"+optKind+"2"
    aktao= createAnalysisObservable( f, algantikt+optR )
    ymax= { "R2":1.0, "R3":0.5, "R4":0.3, "R5":0.3, "R6":0.3 }
    xmax= { "R":1.0, "emin":0.15 }
    plotoptions= { "xmin": 0.0, "xmax": xmax[optKind], "ymin":0.0, "ymax":ymax[optR], "markerStyle": 20, "markerSize": 0.8, "title":"Cone "+optKind+" "+optR+" "+filename }
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

# Compare EEC from various sources with own measurements 
def compareEEC( filename="sjm91_96.root", datafilename="../EECMC/share/OPAL/data.dat" ):

    f= TFile( filename )
    ao= createAnalysisObservable( f, "EEC" )
    tokens= datafilename.split( "/" )
    exp= tokens[3]
    plotoptions= { "xmin": 0.0, "xmax": 3.14159, "ymin": 0.05, "ymax": 5.0, "markerStyle": 20, "markerSize": 0.5, "drawas": "3", "fillcolor": 6, "title": "EEC "+exp, "xlabel": "\chi\ [rad.]", "ylabel": "1/\sigma d\Sigma/d\chi" }
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
    legend.AddEntry( "tgesy", "OPAL 96", "f" )
    legend.Draw()
    
    return 

def compareEECs( filename="sjm91_96.root" ):
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

