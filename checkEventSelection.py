
from ROOT import TChain, TFile

def checkSelection( cfgfile="sjmconfig_192.cfg" ):
    import SjmConfigParser
    configfiles= [ cfgfile, "sjmGeneralOptions.cfg" ]
    config= SjmConfigParser.readConfig( configfiles )
    ecms= config.get( "General", "energy" )
    path= config.get( "General", "path" )
    print( cfgfile )
    testSample( config, "Signal", ecms, path )
    testSample( config, "AltSignal", ecms, path )
    testSample( config, "BkgWWllqq", ecms, path )
    testSample( config, "BkgWWqqqq", ecms, path )
    testSample( config, "BkgWWeeqq", ecms, path )
    return

def testSample( config, sample, ecms, path ):
    filenames= config.get( sample, "files" ).split()
    print( sample, config.get( sample, "name" ) )
    chain= TChain( "h10" )
    for filename in filenames:
        print( filename )
        chain.Add( path+"/"+filename )
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
    
    print( "Standard", nevtStand )
    print( "Costt07", nevtCostt07 )
    print( "Sprold", nevtSprold )
    print( "Wqqlnhi", nevtqqlnhi )
    print( "Wqqlnlo", nevtqqlnlo )
    print( "Wqqqqhi", nevtqqqqhi )
    print( "Wqqqqlo", nevtqqqqlo )

    return


def countEvents( obs="lepthrust" ):

    selectionKeys= [ "stand", "costt07", "sprold", "wqqlnhi", "wqqlnlo", "wqqqqhi", "wqqqqlo" ]
    ecms= [ "91", "130", "136", "161", "172", "183", "189", "192", "196",
                "200", "202", "205", "207" ]

    print( "Number of events for ecm (using "+obs+")" )
    print( " ecm", end="" )
    for selectionKey in selectionKeys:
        print( "{:>10s}".format( selectionKey ), end="" )
    print()

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
        print( "{:>4s}".format( ecm ), end="" )
        for selectionKey in selectionKeys:
            print( "{:10.1f}".format( counters.get( selectionKey ) ), end="" )
        print()
        
    return

def unpackAlphanumericHisto( hist ):
    values= dict()
    xaxis= hist.GetXaxis()
    if xaxis.IsAlphanumeric():
        keys= xaxis.GetLabels()
        for key in keys:
            strKey= key.String().Data()
            ibin= xaxis.FindBin( strKey )
            value= hist.GetBinContent( ibin )
            values[strKey]= value
    else:
        print( "unpackAlphanumericHisto: not alphanumeric ", hist.GetName() )
    return values

def getCutflow( ecm="161" ):
    sampleCounters= dict()
    filename= "sjm"+ecm+".root"
    if ecm == "91":
        filename= "sjm"+ecm+"_all.root"
    f= TFile( filename )
    eventCountsHisto= f.Get( "Eventcounts" )
    eventCounts= unpackAlphanumericHisto( eventCountsHisto )
    keyList= f.GetListOfKeys()
    MCsamples= [ "Signal", "AltSignal", "BkgWWqqqq", "BkgWWllqq", "BkgWWeeqq" ]
    samples= [ "Data" ] + MCsamples
    countersKeys= [ key.GetTitle() for key in keyList if key.GetTitle() in samples ]
    for counterskey in countersKeys:
        hist= f.Get( counterskey )
        if hist:
            counters= unpackAlphanumericHisto( hist )
            counters["events"]= eventCounts[counterskey]
            sampleCounters[counterskey]= counters
    import SjmConfigParser
    cfgFilename= "sjmconfig_"+ecm+".cfg"
    config= SjmConfigParser.readConfig( [ cfgFilename ] )
    lumi= getFromConfig( config, "Data", "lumi" )
    sampleCounters["Data"]["lumi"]= lumi
    for section in MCsamples:
        xsec= getFromConfig( config, section, "xsec" )
        if xsec:
            sampleCounters[section]["xsec"]= xsec
    return sampleCounters

def getFromConfig( config, section, key ):
    import configparser
    result= None
    try:
        result= config.getfloat( section, key )
    except( configparser.NoSectionError, configparser.NoOptionError ):
        # print( "getFromConfig: section or key error", section, key )
        pass
    return result

def printCutflow( ecm="161" ):
    sampleCounters= getCutflow( ecm )
    dataCounter= sampleCounters["Data"]
    lumi= dataCounter["lumi"]
    MCsamples= [ "Signal", "AltSignal", "BkgWWqqqq", "BkgWWllqq", "BkgWWeeqq" ]
    print( "Cutflow for ecm=", ecm )
    print( "cut    ", end="" )
    for sample in [ "Data" ] + MCsamples:
        if sample in sampleCounters:
            print( "{:>10}".format( sample ), end="" )
    print( "{:>10}".format( "MC sum" ) )
    cuts= [ "l2mh", "nch7", "costt09", "sprime", "wqqqq", "wqqln" ]
    if ecm == "91":
        cuts= [ "tkmh", "nch7", "costt09" ]
    for cut in cuts:
        print( "{:7}".format( cut ), end="" )
        print( "{:10.1f}".format( dataCounter[cut] ), end="" )
        mcsum= 0.0
        for sample in MCsamples:
            if sample in sampleCounters:
                sampleCounter= sampleCounters[sample]
                wgt= lumi*sampleCounter["xsec"]/sampleCounter["events"]
                value= sampleCounter[cut] * wgt
                if sample != "AltSignal":
                    mcsum+= value
                print( "{:10.1f}".format( value ), end="" )
        print( "{:10.1f}".format( mcsum ) )
    return

def printCutflows():
    for ecm in [ 130, 136, 161, 172, 183, 189, 192, 196, 200, 202, 205, 207 ]:
        printCutflow( str(ecm) )
    return

