

import re

# Get lumi and efficiency values from logfile
def getValues( logfileName="lu161.log", SWokonly=False, debug=True ):
    effSWpattern= "Eff\(SW ok && MHAD ok\)"
    ieffSW= 5
    iSWRL= 5
    if SWokonly:
        effSWpattern= "Eff\(SW ok\)"
        ieffSW= 2
        iSWRL= 4
    effSWvalues= dict()
    SWRLlumiValues= dict()
    energyranges= dict()
    print "Get values from", logfileName
    with open( logfileName ) as logfile:
        for line in logfile:
            if re.search( "QLPRNT, energy bin  ", line ):
                if debug:
                    print line,
                tokens= line.split()
                energybin= tokens[4]
                energybin= energybin.strip( ":" )
                energyrange= ( int(tokens[5]), int(tokens[7]) )
                energyranges[energybin]= energyrange
                if debug:
                    print "energybin:", energybin
                effSW= 0.0
                SWRLlumi= 0.0
                SWRLlumiSyst= 0.0
                SWRLlumiStat= 0.0
                for i in range( 30 ):
                    nextline= next( logfile )
                    if re.search( "ROCROS sum  \d+ not filled", nextline ):
                        if debug:
                            print "empty bin, break"
                        break
                    # if re.search( "Eff\(SW ok && MHAD ok\)", nextline ):
                    if re.search( effSWpattern, nextline ):
                        tokens= nextline.split()
                        effSW= ( float(tokens[ieffSW]), float(tokens[ieffSW+2]) )
                    if re.search( "SW RL integrated lumi", nextline ):
                        tokens= nextline.split()
                        SWRLlumi= float(tokens[iSWRL])
                    if re.search( "Syst error on SW RL", nextline ):
                        tokens= nextline.split()
                        SWRLlumiSyst= float(tokens[iSWRL+2])
                    if re.search( "Stat error on SW RL", nextline ):
                        tokens= nextline.split()
                        SWRLlumiStat= float(tokens[iSWRL+2])
                if debug:
                    print "Eff(SW ok):", effSW, "energybin:", energybin
                    print "SW RL lumi:", SWRLlumi
                effSWvalues[energybin]= effSW
                SWRLlumiValues[energybin]= ( SWRLlumi, SWRLlumiStat, SWRLlumiSyst )
    if debug:
        for key in sorted( energyranges.keys() ):
            print key, energyranges[key], effSWvalues[key], SWRLlumiValues[key]
    return energyranges, effSWvalues, SWRLlumiValues

# Get lumi and efficiency values for energy bin and calculate correction
def calcLumi( ecm=161, year=None, SWokonly=False, debug=False ):
    logfilename= "lu"+str(ecm)+".log"
    if year:
        logfilename= "lu"+year+".log"
    energyranges, effSWvalues, SWRLlumiValues= getValues( logfilename, SWokonly, debug )
    energybin= None
    ecmMeV= ecm*1000
    for key in sorted( energyranges.keys() ):
        energyrange= energyranges[key]
        if energyrange[0] < ecmMeV and ecmMeV < energyrange[1]:
            energybin= key
            break
    SWRLlumi= SWRLlumiValues[energybin][0]/1000.0
    effSW= effSWvalues[energybin][0]
    lumi= SWRLlumi/effSW
    return ( SWRLlumi, effSW, lumi )

# Lumi calculation for all energy bins
def calcLumis( SWokonly=False ):

    results= dict()
    results["130_95"]= calcLumi( 130, "133_95", SWokonly )
    results["130_97"]= calcLumi( 130, "133_97", SWokonly )
    results["136_95"]= calcLumi( 136, "133_95", SWokonly )
    results["136_97"]= calcLumi( 136, "133_97", SWokonly )
    results["161"]= calcLumi( 161, SWokonly=SWokonly )
    results["172"]= calcLumi( 172, SWokonly=SWokonly )
    results["183"]= calcLumi( 183, SWokonly=SWokonly )
    results["189"]= calcLumi( 189, SWokonly=SWokonly )
    results["192"]= calcLumi( 192, "1999", SWokonly )
    results["196"]= calcLumi( 196, "1999", SWokonly )
    results["200_99"]= calcLumi( 200, "1999", SWokonly )
    results["200_2k"]= calcLumi( 200, "2000", SWokonly )
    results["202_99"]= calcLumi( 202, "1999", SWokonly )
    results["202_2k"]= calcLumi( 202, "2000", SWokonly )
    results["205"]= calcLumi( 205, "2000", SWokonly )
    results["207"]= calcLumi( 207, "2000", SWokonly )

    print "Lumi values per year and energy bin",
    if SWokonly:
        print "(SW ok only)"
    else:
        print
    print "ecm/year SW lumi  eff   corrected"
    for key in sorted( results.keys() ):
        SWRLlumi, effSW, lumi= results[key]
        print "{0:6s}".format( key ), "{0:7.2f}".format( SWRLlumi ), "{0:7.4f}".format( effSW ), "{0:7.2f}".format( lumi )

    sums= dict()
    for ecm in [ 161, 172, 183, 189, 192, 196, 205, 207 ]:
        key= str(ecm)
        SWRLlumi, effSW, lumi= results[key]
        sums[key]= lumi
    sums["130"]= calcLumisum( results, [ "130_95", "130_97" ] )
    sums["136"]= calcLumisum( results, [ "136_95", "136_97" ] )
    sums["133"]= calcLumisum( results, [ "130_95", "130_97", "136_95", "136_97" ] )
    sums["200"]= calcLumisum( results, [ "200_99", "200_2k" ] )
    sums["202"]= calcLumisum( results, [ "202_99", "202_2k" ] )

    print "Lumi values per energy bin"
    for key in sorted( sums.keys() ):
        print "{0:6s}".format( key ), "{0:7.2f}".format( sums[key] )

    return

def calcLumisum( results, keys ):
    lumisum= 0.0
    for key in keys:
        SWRLlumi, effSW, lumi= results[key]
        lumisum+= lumi
    return lumisum
