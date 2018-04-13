
from ROOT import TFile, gROOT, gPad, TVectorD, TObject

from ROOT import TGraphErrors, TH1D, TLegend, TCanvas
TGraphErrors.__init__._creates= False
TH1D.__init__._creates= False
TLegend.__init__._creates= False
TCanvas.__init__._creates= False

gROOT.LoadMacro( "libAnalysisDict.so" )

from ROOT import Analysis, TH1DAnalysisObject, TGEAnalysisObject

from array import array


# Factory method to create AnalysisObject instances
def getAnalysisObjectFromFile( tfile, obs, analysis ):
    obj= tfile.Get( obs+" "+analysis.getTag()+";1" )
    if obj.ClassName() == "TH1D":
        errobj= tfile.Get( "errm "+obs+" "+analysis.getTag() )
        if errobj:
            ao= TH1DAnalysisObject( obj, errobj )
        else:
            ao= TH1DAnalysisObject( obj )
    elif obj.ClassName() == "TGraphErrors":
        ao= TGEAnalysisObject( obj )
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
        return
        
    def printResults( self, width=6, precision=3 ):
        print "Results for", self.obs
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
            print fmt.format( self.sterrs[i] ),
            print fmt.format( self.syerrs[i] )
        return

    def plot( self, plotoptions, opt="?" ):
        vx= array( "d", self.aostand.getPointsCenter() )
        npoints= len(vx)
        vex= array( "d", npoints*[0.0] )
        tgest= TGraphErrors( npoints, vx, self.values, vex, self.sterrs )
        tgesy= TGraphErrors( npoints, vx, self.values, vex, self.syerrs )
        tgesy.SetMarkerStyle( plotoptions["markerStyle"] )
        tgesy.SetMarkerSize( plotoptions["markerSize"] )
        drawas= plotoptions["drawas"] if "drawas" in plotoptions else "p"
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
        gPad.SetLogy()
        tgest.Draw( "same"+drawas )
        return tgest, tgesy

# LEP1 Analysis
class LEP1AnalysisObservable( AnalysisObservable ):

    def __init__( self, obs, tfile, unf="bbb" ):
        AnalysisObservable.__init__( self, obs )
        standardAnalysis= Analysis( "data mt stand none none none py " + unf )
        tcAnalysis= Analysis( "data tc stand none none none py " + unf )
        costt07Analysis= Analysis( "data mt costt07 none none none py " + unf )
        nch7Analysis= Analysis( "data mt nch7 none none none py " + unf )
        hwAnalysis= Analysis( "data mt stand none none none hw " + unf )
        self.aostand= getAnalysisObjectFromFile( tfile, obs, standardAnalysis )
        self.points= array( "d", self.aostand.getPoints() )
        self.values= array( "d", self.aostand.getValues() )
        self.sterrs= array( "d", self.aostand.getErrors() )
        tctvd= array( "d", getAnalysisObjectFromFile( tfile, obs, tcAnalysis ).getValues() )
        costt07tvd= array( "d", getAnalysisObjectFromFile( tfile, obs, costt07Analysis ).getValues() )
        nch7tvd= array( "d", getAnalysisObjectFromFile( tfile, obs, nch7Analysis ).getValues() )
        hwtvd= array( "d", getAnalysisObjectFromFile( tfile, obs, hwAnalysis ).getValues() )
        from numpy import square, sqrt, subtract
        tcdelta= subtract( tctvd, self.values )
        costt07delta= subtract( costt07tvd, self.values )
        nch7delta= subtract( nch7tvd, self.values )
        hwdelta= subtract( hwtvd, self.values )
        self.syerrs= square(tcdelta)+square(costt07delta)+square(nch7delta)+square(hwdelta)
        self.syerrs= sqrt(self.syerrs)
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
            
    
# LEP2 Analysis
class LEP2AnalysisObservable( AnalysisObservable ):

    def __init__( self, obs, tfile ):
        AnalysisObservable.__init__( self, obs )
        standardAnalysis= Analysis( "data mt stand none none llqq:qqqq:eeqq py bbb" )
        tcAnalysis= Analysis( "data tc stand none none llqq:qqqq:eeqq py bbb" )
        costt07Analysis= Analysis( "data mt costt07 none none llqq:qqqq:eeqq py bbb" )

        hwAnalysis= Analysis( "data mt stand none none llqq:qqqq:eeqq hw bbb" )
        self.aostand= getAnalysisObjectFromFile( tfile, obs, standardAnalysis )
        self.points= array( "d", self.aostand.getPoints() )
        self.values= array( "d", self.aostand.getValues() )
        self.sterrs= array( "d", self.aostand.getErrors() )
        tctvd= array( "d", getAnalysisObjectFromFile( tfile, obs, tcAnalysis ).getValues() )
        costt07tvd= array( "d", getAnalysisObjectFromFile( tfile, obs, costt07Analysis ).getValues() )

        hwtvd= array( "d", getAnalysisObjectFromFile( tfile, obs, hwAnalysis ).getValues() )
        from numpy import square, sqrt, subtract
        tcdelta= subtract( tctvd, self.values )
        costt07delta= subtract( costt07tvd, self.values )

        hwdelta= subtract( hwtvd, self.values )
        self.syerrs= square(tcdelta)+square(costt07delta)+square(hwdelta)
        self.syerrs= sqrt(self.syerrs)
        return



# Compare antikt, siscone and PXCONE jets in same plot        
def compareConeRjets( filename="sjm91_96_test.root" ):
    f= TFile( filename )
    ao= LEP1AnalysisObservable( "antiktRR3", f )
    plotoptions= { "xmin": 0.0, "xmax": 1.5, "ymin":0.0, "ymax":0.5, "markerStyle": 20, "markerSize": 0.75 }
    ao.plot( plotoptions )
    ao= LEP1AnalysisObservable( "sisconeRR3", f )
    plotoptions["markerStyle"]= 21
    ao.plot( plotoptions, "s" )
    ao= LEP1AnalysisObservable( "pxconeR2R3", f )
    plotoptions["markerStyle"]= 22
    ao.plot( plotoptions, "s" )    
    return
def compareConeEminjets( filename="sjm91_96_test.root" ):
    f= TFile( filename )
    ao= LEP1AnalysisObservable( "antikteminR3", f )
    plotoptions= { "xmin": 0.0, "xmax": 0.2, "ymin":0.0, "ymax":0.5, "markerStyle": 20, "markerSize": 0.75 }
    ao.plot( plotoptions )
    ao= LEP1AnalysisObservable( "sisconeeminR3", f )
    plotoptions["markerStyle"]= 21
    ao.plot( plotoptions, "s" )
    ao= LEP1AnalysisObservable( "pxconeemin2R3", f )
    plotoptions["markerStyle"]= 22
    ao.plot( plotoptions, "s" )    
    return

# Compare EEC from various sources with own measurements 
def compareEEC( filename="sjm91_96.root", datafilename="../EECMC/share/OPAL/data.dat" ):

    f= TFile( filename )
    ao= LEP1AnalysisObservable( "EEC", f )
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

