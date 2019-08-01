
from ROOT import TChain, TFile
    

def checkSelection( cfgfile="sjmconfig_192.cfg" ):

    import SjmConfigParser

    configfiles= [ cfgfile, "sjmGeneralOptions.cfg" ]
    config= SjmConfigParser.readConfig( configfiles )

    ecms= config.get( "General", "energy" )
    print cfgfile
    
    testSample( config, "Signal", ecms )
    testSample( config, "AltSignal", ecms )

    testSample( config, "BkgWWllqq", ecms )
    testSample( config, "BkgWWqqqq", ecms )
    testSample( config, "BkgWWeeqq", ecms )
    
    return

    
def testSample( config, sample, ecms ):
    filenames= config.get( sample, "files" ).split()
    print sample, config.get( sample, "name" )
    chain= TChain( "h10" )
    for filename in filenames:
        print filename
        chain.Add( filename )
    testChain( chain, ecms )
    return

        
def testChain( chain, ecms ):

    ecms+= ".0"
    
    nevtStand= chain.Draw( "Tdtc", "Icjst==3&&Iebst==3&&Il2mh==1&&Ntkd02>6&&fabs(Tvectc[2])<0.9&&"+ecms+"-sqrt( fabs( Pspr[3]*Pspr[3]-Pspr[2]*Pspr[2]-Pspr[1]*Pspr[1]-Pspr[0]*Pspr[0] ) )<10.0&&Lwqqln<0.5&&Lwqqqq<0.25" )
    nevtCostt07= chain.Draw( "Tdtc", "Icjst==3&&Iebst==3&&Il2mh==1&&Ntkd02>6&&fabs(Tvectc[2])<0.7&&"+ecms+"-sqrt( fabs( Pspr[3]*Pspr[3]-Pspr[2]*Pspr[2]-Pspr[1]*Pspr[1]-Pspr[0]*Pspr[0] ) )<10.0&&Lwqqln<0.5&&Lwqqqq<0.25" )
    nevtSprold= chain.Draw( "Tdtc", "Icjst==3&&Iebst==3&&Il2mh==1&&Ntkd02>6&&fabs(Tvectc[2])<0.9&&"+ecms+"-sqrt( fabs( Pspri[3]*Pspri[3]-Pspri[2]*Pspri[2]-Pspri[1]*Pspri[1]-Pspri[0]*Pspri[0] ) )<10.0&&Lwqqln<0.5&&Lwqqqq<0.25" )

    nevtqqlnhi= chain.Draw( "Tdtc", "Icjst==3&&Iebst==3&&Il2mh==1&&Ntkd02>6&&fabs(Tvectc[2])<0.9&&"+ecms+"-sqrt( fabs( Pspr[3]*Pspr[3]-Pspr[2]*Pspr[2]-Pspr[1]*Pspr[1]-Pspr[0]*Pspr[0] ) )<10.0&&Lwqqln<0.75&&Lwqqqq<0.25" )
    nevtqqlnlo= chain.Draw( "Tdtc", "Icjst==3&&Iebst==3&&Il2mh==1&&Ntkd02>6&&fabs(Tvectc[2])<0.9&&"+ecms+"-sqrt( fabs( Pspr[3]*Pspr[3]-Pspr[2]*Pspr[2]-Pspr[1]*Pspr[1]-Pspr[0]*Pspr[0] ) )<10.0&&Lwqqln<0.25&&Lwqqqq<0.25" )
    nevtqqqqhi= chain.Draw( "Tdtc", "Icjst==3&&Iebst==3&&Il2mh==1&&Ntkd02>6&&fabs(Tvectc[2])<0.9&&"+ecms+"-sqrt( fabs( Pspr[3]*Pspr[3]-Pspr[2]*Pspr[2]-Pspr[1]*Pspr[1]-Pspr[0]*Pspr[0] ) )<10.0&&Lwqqln<0.5&&Lwqqqq<0.4" )
    nevtqqqqlo= chain.Draw( "Tdtc", "Icjst==3&&Iebst==3&&Il2mh==1&&Ntkd02>6&&fabs(Tvectc[2])<0.9&&"+ecms+"-sqrt( fabs( Pspr[3]*Pspr[3]-Pspr[2]*Pspr[2]-Pspr[1]*Pspr[1]-Pspr[0]*Pspr[0] ) )<10.0&&Lwqqln<0.5&&Lwqqqq<0.1" )
    print "Standard", nevtStand
    print "Costt07", nevtCostt07
    print "Sprold", nevtSprold
    print "Wqqlnhi", nevtqqlnhi
    print "Wqqlnlo", nevtqqlnlo
    print "Wqqqqhi", nevtqqqqhi
    print "Wqqqqlo", nevtqqqqlo

    return


def countEvents( obs="thrust" ):

    selectionKeys= [ "stand", "costt07", "sprold", "wqqlnhi", "wqqlnlo", "wqqqqhi", "wqqqqlo" ]
    ecms= [ "91", "130", "136", "161", "172", "183", "189", "192", "196",
                "200", "202", "205", "207" ]

    print "Number of events for ecm (using "+obs+")"
    print " ecm",
    for selectionKey in selectionKeys:
        print "{:>8s}".format( selectionKey ),
    print

    for ecm in ecms:
        filename= "sjm"+ecm+".root"
        if ecm == "91":
            filename= "sjm91_all.root"
        f= TFile( filename )
        counters= dict()
        for selectionKey in selectionKeys:
            hkey= obs+" data mt "+selectionKey
            hist= f.Get( hkey )
            if hist:
                counters[selectionKey]= hist.GetEntries()
            else:
                counters[selectionKey]= 0
        print "{:>4s}".format( ecm ),
        for selectionKey in selectionKeys:
            print "{:8.1f}".format( counters.get( selectionKey ) ),
        print
        
    return

def cutFlow( ecms= [ "130" ] ):

    for ecm in ecms:
        filename= "sjm"+ecm+"_test.root"
        f= TFile( filename )
        keyList= f.GetListOfKeys()
        countersKeys= [ key.GetTitle() for key in keyList if "Cutflow_" in key.GetTitle() ]
        for counterskey in countersKeys:
            print counterskey
            hist= f.Get( counterskey )
            if hist:
                xaxis= hist.GetXaxis()
                counterKeys= xaxis.GetLabels()
                for counterKey in counterKeys:
                    ibin= xaxis.FindBin( counterKey.String().Data() )
                    counter= hist.GetBinContent( ibin )
                    print counterKey, counter

    return

    
