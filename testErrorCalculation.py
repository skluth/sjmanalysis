
from ROOT import TH1D, TH2D, TF1, TMath, TCanvas, SetOwnership, TRandom3, gRandom, gStyle

TCanvas.__init__._creates= False

import numpy as np

from math import sqrt

binEdges= np.array( [ 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.09, 0.12, 0.15, 0.22, 0.30, 0.50 ] )
nbin= len(binEdges)-1

def inBinRange( value ):
    return binEdges[0] < value and value < binEdges[nbin]

def testWeighted( mnevent=1000, nsample=1000, resolutionFactor=0.02, nmom=1, width=7, precision=4 ):

    # Prepare
    gRandom= TRandom3()
    fun= TF1( "fun", "1000*x*exp(-x*20)", 0.0 , 0.5 )
    histh= TH1D( "histh", "1000*x*exp(-x*20) h", nbin, binEdges )
    histd= TH1D( "histd", "1000*x*exp(-x*20) d", nbin, binEdges )
    histm= TH2D( "histm", "migration h vs d", nbin, binEdges, nbin, binEdges )
    histhw= TH1D( "histhw", "1000*x*exp(-x*20) w h", nbin, binEdges )
    histdw= TH1D( "histdw", "1000*x*exp(-x*20) w d", nbin, binEdges )
    histmwh= TH2D( "histmwh", "migration h vs d w h", nbin, binEdges, nbin, binEdges )
    histmwd= TH2D( "histmwd", "migration h vs d w d", nbin, binEdges, nbin, binEdges )
    histhw2= TH1D( "histhw2", "1000*x*exp(-x*20) w2 h", nbin, binEdges )
    histdw2= TH1D( "histdw2", "1000*x*exp(-x*20) w2 d", nbin, binEdges )
    histmw2h= TH2D( "histmw2h", "migration h vs d w2 h", nbin, binEdges, nbin, binEdges )
    histmw2d= TH2D( "histmw2d", "migration h vs d w2 d", nbin, binEdges, nbin, binEdges )
    SetOwnership( histh, False )
    SetOwnership( histd, False )
    SetOwnership( histm, False )
    SetOwnership( histhw, False )
    SetOwnership( histdw, False )
    SetOwnership( histmwh, False )
    SetOwnership( histmwd, False )
    SetOwnership( histhw2, False )
    SetOwnership( histdw2, False )
    SetOwnership( histmw2h, False )
    SetOwnership( histmw2d, False )

    # Generate and fill histos
    for i in range( mnevent ):
        value= fun.GetRandom()
        histh.Fill( value )
        resValue= resolution( value, resolutionFactor )
        histd.Fill( resValue )
        histm.Fill( resValue, value )
        histhw.Fill( value, value**nmom )
        histdw.Fill( resValue, resValue**nmom )
        histmwh.Fill( resValue, value, value**nmom )
        histmwd.Fill( resValue, value, resValue**nmom )
        nmom2= nmom*2
        histhw2.Fill( value, value**nmom2 )
        histdw2.Fill( resValue, resValue**nmom2 )
        histmw2h.Fill( resValue, value, value**nmom2 )
        histmw2d.Fill( resValue, value, resValue**nmom2 )
        
    # Print raw migration matrices
    print "Migration of events hadron detector level"
    migrMatrix= getMatrix( histm )
    valuesh= getArray( histh )
    valuesd= getArray( histd )
    printMatrixVectors( migrMatrix, valuesd, valuesh )
    print "Migration weighted hadron level"
    migrMatrixwh= getMatrix( histmwh )
    valueshw= getArray( histhw )
    valuesdw= getArray( histdw )
    printMatrixVectors( migrMatrixwh, valuesdw, valueshw )
    print "Migration weighted detector level"
    migrMatrixwd= getMatrix( histmwd )
    printMatrixVectors( migrMatrixwd, valuesdw, valueshw )
    
    # Calculate "unfolding" matrix detector->hadron level
    # Normalised weighted det-level migration
    migrMatrixwdNorm= normaliseColumns( migrMatrixwd )
    # print "Migration weighted detector level normalised column-wise"
    # printMatrix( migrMatrixwdNorm )
    # Ratio of weighted had- to det-level migration matrices
    hdratioMatrix= divideChecked( migrMatrixwh, migrMatrixwd )
    # print "Migration weighted ratio hadron/detector"
    # printMatrix( hdratioMatrix )
    # Multiply normalised weighted det-level migration by had- to det-level ratio
    # of weighted migrations
    migrMatrixwRec= hdratioMatrix*migrMatrixwdNorm
    # print "Migration weighted detector->hadron corrected normalised column-wise"
    # printMatrix( migrMatrixwRec )
    # Efficiency at hadron level correction
    histmwy= histmwh.ProjectionY( "Migr weighted h y proj", 1, nbin )
    cacc= divideHistos( histhw, histmwy )
    # print "Effciciency correction hadron level", cacc
    # Final detector->hadron correction matrix and closure test
    CR= np.dot( np.diag( cacc ), migrMatrixwRec )
    valueshwRec= np.dot( CR, valuesdw )
    print "Closure test weighted distributions reconstructed vs hadron level"
    print valueshwRec
    print valueshw
    # Same for moment**2
    migrMatrixw2h= getMatrix( histmw2h )
    migrMatrixw2d= getMatrix( histmw2d )
    migrMatrixw2dNorm= normaliseColumns( migrMatrixw2d )
    hdratioMatrix2= divideChecked( migrMatrixw2h, migrMatrixw2d )
    migrMatrixw2Rec= hdratioMatrix2*migrMatrixw2dNorm
    histmw2y= histmw2h.ProjectionY( "Migr w2 h y proj", 1, nbin )
    cacc2= divideHistos( histhw2, histmw2y )
    CR2= np.dot( np.diag( cacc2 ), migrMatrixw2Rec )
    valuesdw2= getArray( histdw2 )
    valueshw2Rec= np.dot( CR2, valuesdw2 )
    
    # Calculate folding matrix hadron->detector level 
    # Normalised (row-wise) weighted had-level migration
    migrMatrixwhNorm= normaliseRows( migrMatrixwh )
    #print "Migration weighted had level normalised row-wise"
    #printMatrix( migrMatrixwhNorm )
    # Ratio of weighted had- to det-level migration matrices
    hdratioMatrix2= divideChecked( migrMatrixwd, migrMatrixwh )
    #print "Migration weighted ratio detector/hadron"
    #printMatrix( hdratioMatrix2 )
    # Multiply normalised weighted det-level migration by had- to det-level ratio
    # of weighted migrations
    migrMatrixwh2d= hdratioMatrix2*migrMatrixwhNorm
    #print "Migration weighted detector->hadron corrected normalised column-wise"
    #printMatrix( migrMatrixwh2d )
    # Efficiency at hadron level correction
    caccd= np.divide( 1.0, cacc )
    # print "Effciciency correction detector level", caccd
    # Final hadron-> correction matrix and closure test
    CRh2d= np.dot( migrMatrixwh2d.transpose(), np.diag( caccd ) )
    valuesdwRec= np.dot( CRh2d, valueshw )
    print "Closure test weighted distributions folded vs detector level"
    print valuesdwRec
    print valuesdw
    
    # Put results with propagated errors into histo for normalisation
    histhrec= TH1D( "histhrec", "1000*x*exp(-x*20) h rec", nbin, binEdges )
    for i in range( nbin ):
        histhrec.SetBinContent( i+1, valueshwRec[i] )
    histhrec.SetEntries( mnevent )
    errorsdw= getArrayErrors( histdw )
    diagErrorMatrixdw= np.diag( errorsdw**2 )
    errorMatrixhw= np.dot( CR, np.dot( diagErrorMatrixdw, CR.transpose() ) )
    errorMatrixhwNorm= normalise( histhrec, errorMatrixhw )

    # Print results
    fmt="{:"+str(width)+"."+str(precision)+"f}"
    print "Result corrected normalised weighted distribution"
    for i in range( nbin ):
        valuehwRecNorm= histhrec.GetBinContent( i+1 )
        errorhwRecNorm= sqrt( errorMatrixhwNorm[i,i] )
        histhrec.SetBinError( i+1, errorhwRecNorm )
        print fmt.format( valuehwRecNorm ), fmt.format( errorhwRecNorm )
    print "Correlation matrix"
    corr= cov2corr( errorMatrixhwNorm )
    printMatrix( corr )

    # Integral of corrected weighted is 1st moment of unweighted,
    # error from diagonal bin errors only
    sumw= calcIntegral( histhrec )
    errors= getArrayErrors( histhrec )
    binw= binEdges[1:]-binEdges[:nbin]
    errorsbw= errors*binw
    errSumw= sqrt( sum(errorsbw**2) )

    # Results from unnormalised distribution and errors from covariance or mom*2 
    sumwunn= np.sum( valueshwRec ) / float(mnevent)
    errSumwunn= sqrt( np.sum( errorMatrixhw ) ) / float(mnevent)
    sumw2unn= np.sum( valueshw2Rec ) / float(mnevent)
    # "Error of mean is s.d./sqrt(N)", i.e. from moment*2
    errSumw2unn= sqrt( ( sumw2unn - sumwunn**2 ) / float(mnevent) )

    # # Error of integral given by Sum V_ij*binw_i*binw_j
    # # but does not work in single closure sample
    # errSumw= 0.0
    # for i in range( nbin-1 ):
    #     for j in range( nbin-1 ):
    #         errSumw+= errorMatrixhwNorm[i,j]*binw[i]*binw[j]
    
    fmt="{:"+str(width)+"."+str(precision)+"f}"
    print "Moment", nmom, "from integrated corrected weighted distribution:",
    print fmt.format( sumw ), "+/-", fmt.format( errSumw ), "(diag.)"
    print "Moment", nmom, "from unnormalised weighted distribution:",
    print fmt.format( sumwunn ), "+/-", fmt.format( errSumwunn ), "(covm)",
    print fmt.format( errSumw2unn ), "(s.d.)"
    print "Mean of had-level histo from ROOT:",
    print fmt.format( histh.GetMean() ), "+/-", fmt.format( histh.GetMeanError() )

    # Plots
    canv= TCanvas( "canv", "Weighted correction", 800, 1000 )
    canv.Divide( 2, 3 )
    canv.cd( 1 )
    normalise( histhw )
    histhw.Draw()
    canv.cd( 2 )
    normalise( histdw )
    histdw.Draw()
    canv.cd( 3 )
    histm.Draw( "box" )
    canv.cd( 4 )
    histmwh.Draw( "box" )
    canv.cd( 5 )
    histmwd.Draw( "box" )
    canv.cd( 6 )
    histh.Draw()
        
    return

def divideHistos( hnom, hdenom ):
    result= np.array( nbin*[0.0] )
    for i in range( nbin ):
        result[i]= hnom.GetBinContent( i+1 ) / hdenom.GetBinContent( i+1 )
    return result

def divideChecked( nom, denom ):
    shapeNom= nom.shape
    shapeDenom= denom.shape
    if shapeNom != shapeDenom:
        raise "divideChecked: shapes not equal"
    result= np.array( np.zeros( shapeNom ) )
    if len(shapeNom) == 1:
        for i in range( shapeNom[0] ):
            result[i]= divideOrZero( nom[i], denom[i] )
    if len(shapeNom) == 2:
        for i in range( shapeNom[0] ):
            for j in range( shapeNom[1] ):
                result[i,j]= divideOrZero( nom[i,j], denom[i,j] )
    return result

def divideOrZero( nom, denom ):
    if denom == 0.0:
        return 0.0
    else:
        return nom/denom

def printMatrixVectors( m, vx, vy, width=7, precision=3 ):
    fmt="{:"+str(width)+"."+str(precision)+"f}"
    for i in range( nbin ):
        print fmt.format( vy[i] ),
        for j in range( nbin ):
            print fmt.format( m[i,j] ),
        print
    print "       ",
    for i in range( nbin ):
        print fmt.format( vx[i] ),
    print
    return

def calcIntegral( hist ):
    nbin= hist.GetNbinsX()
    sumw= 0.0
    for i in range( nbin ):
        sumw+= hist.GetBinContent( i+1 )*hist.GetBinWidth( i+1 )
    return sumw

def getMatrix( hist2d ):
    matrix= np.array( np.zeros( (nbin,nbin) ) )
    for i in range(nbin):
        for j in range(nbin):
            matrix[j,i]= hist2d.GetBinContent( i+1, j+1 )
    return matrix

def getArray( hist ):
    arr= np.array( nbin*[0.0] )
    for i in range(nbin):
        arr[i]= hist.GetBinContent( i+1 )
    return arr
def getArrayErrors( hist ):
    arr= np.array( nbin*[0.0] )
    for i in range(nbin):
        arr[i]= hist.GetBinError( i+1 )
    return arr

# Test sampling
def testSampling( mnevent=1000, nsamples=1000, opt="p" ):
    print "Number of samples:", nsamples
    if "p" in opt:
        print "Variable number of events from Poisson mu=", mnevent
    else:
        print "Fixed number of events n=", mnevent
    if "w1" in opt:
        print "Weighted distribution w1"
    elif "w2" in opt:
        print "Weighted distribution w2"
    gRandom= TRandom3()
    fun= TF1( "fun", "1000*x*exp(-x*20)", 0.0 , 0.5 )
    hists= dict()
    for isample in range( nsamples ):
        hist= TH1D( "histh"+str(isample), "1000*x*exp(-x*20) h", nbin, binEdges )
        if "p" in opt:
            nevent= gRandom.Poisson( mnevent )
        else:
            nevent= mnevent
        for i in range( int( nevent ) ):
            value= fun.GetRandom()
            nmom= nmomFromOpt( opt )
            # fillweight= weight( value, opt )
            fillweight= value**nmom
            hist.Fill( value, fillweight )
        hists[isample]= hist
    errorMatrixSample= sampleErrorMatrix( hists )
    print "Error matrix"
    printMatrix( errorMatrixSample, 7, 4 )
    corr= cov2corr( errorMatrixSample )
    print "Correlation matrix"
    printMatrix( corr, 6, 3 )
    return

# Sampling tests of error calculation
def sample( nsample=1000, acc=0.8, resolutionFactor=0.02, mnevent=1000, opt="sme",
            width=7, precision=4 ):

    if "s" in opt:
        print "Sample from truth and gaussian smearing"
    elif "d" in opt:
        print "Sample from reference detector level"
    if "m" in opt:
        print "Use matrix correction"
    elif "b" in opt:
        print "Use bin-by-bin correction"
    if "e" in opt:
        print "Use x*exp(-x) function"
    elif "g" in opt:
        print "Use Gaus function"
    if "u" in opt:
        print "Use expected errors for sampled detector level"
    import re
    mo= re.search( "w[0-5]", opt )
    if mo:
        print "Use weighted distribution", mo.group()

    # Initialise ROOT random numbers
    gRandom= TRandom3()

    # Function to sample
    if "e" in opt:
        fun= TF1( "fun", "1000*x*exp(-x*20)", 0.0 , 0.5 )
    elif "g" in opt:
        fun= TF1( "fun", "TMath::Gaus(x,0.25,0.15)", 0.0 , 0.5 )
    fun.SetNormalized(True)
    
    # Correction factors, matrix and reference histo
    cf, CR, cf2, CR2, eventAcc, histref, histdet, histref2, histdet2= generate( fun, 1000000, acc, resolutionFactor, opt )

    # Create and fill sample histos
    if "s" in opt:
        histds, histd2s, hisths= fillSampleHistosSmear( fun, nsample, mnevent,
                                                        resolutionFactor, opt )
    elif "d" in opt:
        histds, histd2s, hisths= fillSampleHistosDraw( histref, histdet, histdet2, nsample,
                                                       mnevent, fun.GetTitle() )

    # Error matrices after filling
    errorMatrixhSample= sampleErrorMatrix( hisths )
    errorMatrixdSample= sampleErrorMatrix( histds )
    
    # Get results with migration matrix or bin-by-bin correction
    if "b" in opt:
        histcs, errorMatricesCorr= bbbCorrection( histds, cf, eventAcc,
                                                  histdet, opt )
        histc2s, errorMatricesCorr2= bbbCorrection( histd2s, cf2, eventAcc,
                                                    histdet2, opt )
    elif "m" in opt:
        histcs, errorMatricesCorr= matrixCorrection( histds, cf, CR, eventAcc,
                                                     histdet, opt )
        histc2s, errorMatricesCorr2= matrixCorrection( histd2s, cf2, CR2, eventAcc,
                                                       histdet2, opt )

    # Error matrices after correction
    errorMatrixCorrAvg= averageMatrix( errorMatricesCorr )
    errorMatrixCorrSample= sampleErrorMatrix( histcs )

    print "Error matrix hadron"
    printMatrix( errorMatrixhSample, 7, 4 )

    corr= cov2corr( errorMatrixhSample )
    print "Correlation matrix hadron"
    printMatrix( corr )
    
    print "Error matrix detector"
    printMatrix( errorMatrixdSample, 7, 4 )

    corr= cov2corr( errorMatrixdSample )
    print "Correlation matrix detector"
    printMatrix( corr )

    print "Error matrix after correction calculated avg."
    printMatrix( errorMatrixCorrAvg, 7, 4 )
    print "Error matrix after correction sampled"
    printMatrix( errorMatrixCorrSample, 7, 4 )

    # Normalise
    errorMatricesNorm= dict()
    for isample in range( nsample ):
        histc= histcs[isample]
        # Take error matrix from sampling or calculated
        errorMatrix= errorMatricesCorr[isample]
        # errorMatrix= errorMatrixCorrAvg
        # errorMatrix= errorMatrixCorrSample
        errorMatricesNorm[isample]= normalise( histc, errorMatrix )
        histh= hisths[isample]
        normalise( histh )
    
    # Sample averages        
    mvalues= np.array( nbin*[0.0] )
    merrors= np.array( nbin*[0.0] )
    mvaluesh= np.array( nbin*[0.0] )
    merrorsh= np.array( nbin*[0.0] )
    merrors2= np.array( nbin*[0.0] )
    mmomentc= 0.0
    mmomentc2= 0.0
    mmomerrc= 0.0
    mmomenth= 0.0
    mmomenth2= 0.0
    for isample in range( nsample ):
        histc= histcs[isample]
        histc2= histc2s[isample]
        histh= hisths[isample]
        for i in range( nbin ):
            valueh= histh.GetBinContent( i+1 )
            mvaluesh[i]+= valueh
            merrorsh[i]+= valueh**2
            valuei= histc.GetBinContent( i+1 )
            mvalues[i]+= valuei
            merrors[i]+= valuei**2
            error= histc.GetBinError( i+1 )
            merrors2[i]+= error
        momentc= calcIntegral( histc )
        momentc2= histc2.Integral( 1, nbin )/histc2.GetEntries()
        mmomentc+= momentc
        mmomentc2+= momentc**2
        mmomerrc+= np.sqrt( np.sum( errorMatricesCorr[isample] ) )/histc.GetEntries()
        # mmomerrc+= sqrt( ( momentc2 - momentc**2 ) / histc2.GetEntries() )
        momenth= calcIntegral( histh )
        mmomenth+= momenth
        mmomenth2+= momenth**2
    mvaluesh/= float(nsample)
    merrorsh/= float(nsample)
    mvalues/= float(nsample)
    merrors/= float(nsample)
    merrors-= mvalues**2
    merrors= np.sqrt( merrors )
    merrorsh-= mvaluesh**2
    merrorsh= np.sqrt( merrorsh )
    merrors2/= float(nsample)
    mmomentc/= float(nsample)
    mmomentc2/= float(nsample)
    mmomerrc/= float(nsample)
    mmomenth/= float(nsample)
    mmomenth2/= float(nsample)

    errorMatrixNormSample= sampleErrorMatrix( histcs )
    errorMatrixNormAvg= averageMatrix( errorMatricesNorm )
    
    # Reference from histref
    normalise( histref )
        
    # Pulls
    pulls= np.array( nbin*[0.0] )
    pullsstd= np.array( nbin*[0.0] )
    pullmom= 0.0
    pullmom2= 0.0
    for isample in range( nsample ):
        histc= histcs[isample]
        # errorMatrix= errorMatricesNorm[isample]
        # errorMatrix= errorMatrixNormAvg
        errorMatrix= errorMatrixNormSample
        for i in range( nbin ):
            if errorMatrix[i,i] != 0.0:
                pull= ( ( histc.GetBinContent( i+1 ) - histref.GetBinContent( i+1 ) ) /
                        sqrt( errorMatrix[i,i] ) )
                pulls[i]+= pull
                pullsstd[i]+= pull**2
        # Error from covariance matrix or from sampling
        pullmerr= np.sqrt( np.sum( errorMatricesCorr[isample] ) )/histc.GetEntries()
        # pullmerr= sqrt( mmomentc2 - mmomentc**2 )
        histc2= histc2s[isample]
        moment= calcIntegral( histc )
        nevent= histc2.GetEntries()
        moment2= histc2.Integral( 1, nbin )/nevent
        # pullmerr= sqrt( ( moment2 - moment**2 ) / nevent )
        pullm= ( moment - calcIntegral( histref ) ) / pullmerr
        pullmom+= pullm
        pullmom2+= pullm**2
    pulls/= float(nsample)
    pullsstd/= float(nsample)
    pullsstd-= pulls**2
    pullsstd= np.sqrt( pullsstd )
    pullmom/= float(nsample)
    pullmom2/= float(nsample)
            
    # Chi^2s w.r.t. reference
    hchi2diag= TH1D( "hchi2diag", "P(Chi^2) diag. errors", 20, 0.0, 1.0 )
    hchi2matrix= TH1D( "hchi2matrix", "P(Chi^2) error matrix calc", 20, 0.0, 1.0 )
    hchi2matrixs= TH1D( "hchi2matrixs", "P(Chi^2) error matrix sample", 20, 0.0, 1.0 )
    hchi2matrixa= TH1D( "hchi2matrixa", "P(Chi^2) error matrix calc avg", 20, 0.0, 1.0 )
    SetOwnership( hchi2diag, False )
    SetOwnership( hchi2matrix, False )
    SetOwnership( hchi2matrixs, False )
    SetOwnership( hchi2matrixa, False )
    for isample in range( nsample ):
        histc= histcs[isample]
        chi2diag= 0.0
        for i in range( 1, nbin+1 ):
            error= histc.GetBinError( i )
            if error > 0.0:
                chi2diag+= ( ( histc.GetBinContent( i ) - histref.GetBinContent( i ) ) /
                                 error )**2                
        hchi2diag.Fill( TMath.Prob( chi2diag, nbin ) )
        chi2matrix= chi2Matrix( errorMatricesNorm[isample], histc, histref )
        hchi2matrix.Fill( TMath.Prob( chi2matrix, nbin-1 ) )
        chi2matrix= chi2Matrix( errorMatrixNormSample, histc, histref  )
        hchi2matrixs.Fill( TMath.Prob( chi2matrix, nbin-1 ) )
        chi2matrix= chi2Matrix( errorMatrixNormAvg, histc, histref  )
        hchi2matrixa.Fill( TMath.Prob( chi2matrix, nbin-1 ) )

    # Printing
    fmt="{:"+str(width)+"."+str(precision)+"f}"
    print "Sample averages"
    print "Avg. rec. moment:", fmt.format(mmomentc), "+/-", fmt.format(mmomerrc), "(avg.)",
    print fmt.format(sqrt( mmomentc2 - mmomentc**2 )), "(sample)"
    print "Avg. had. moment:", fmt.format(mmomenth), "+/-", fmt.format(sqrt( mmomenth2 - mmomenth**2 ))
    print "Ref. moment, pull:", fmt.format(histref.Integral("width")), fmt.format(pullmom), fmt.format(sqrt(pullmom2-pullmom**2))
    print "val.corr.err.sam.err.avg.M.sam.  M.calc. val.had.err.had.ref.    pull    pull.std."
    for i in range( nbin ):
        print fmt.format( mvalues[i] ),
        print fmt.format( merrors[i] ),
        print fmt.format( merrors2[i] ),
        print fmt.format( sqrt( errorMatrixNormSample[i][i] ) ),
        print fmt.format( sqrt( errorMatrixNormAvg[i][i] ) ),
        print fmt.format( mvaluesh[i] ),
        print fmt.format( merrorsh[i] ),
        print fmt.format( histref.GetBinContent( i+1 ) ),
        print fmt.format( pulls[i] ),
        print fmt.format( pullsstd[i] )
    corr= cov2corr( errorMatrixNormSample )
    print "Correlation from sampled error matrix"
    printMatrix( corr )
    corr= cov2corr( errorMatrixNormAvg )
    print "Correlation from average calculated error matrix"
    printMatrix( corr )

    # Plots
    gStyle.SetErrorX( 0.0 )
    canv1= TCanvas( "canv1", "Chi^2", 800, 1000 )
    canv1.Divide( 2, 3 )
    canv1.cd( 1 )
    hchi2diag.SetMarkerStyle( 20 )
    hchi2diag.SetMarkerSize( 0.5 )
    hchi2diag.Draw()
    canv1.cd( 2 )
    hchi2matrix.SetMarkerStyle( 20 )
    hchi2matrix.SetMarkerSize( 0.5 )
    hchi2matrix.Draw()
    canv1.cd( 3 )
    hchi2matrixs.SetMarkerStyle( 20 )
    hchi2matrixs.SetMarkerSize( 0.5 )
    hchi2matrixs.Draw()
    canv1.cd( 4 )
    hchi2matrixa.SetMarkerStyle( 20 )
    hchi2matrixa.SetMarkerSize( 0.5 )
    hchi2matrixa.Draw()
    
    canv1.cd( 5 )
    SetOwnership( histref, False )
    histref.SetMarkerStyle( 20 )
    histref.SetMarkerSize( 0.5 )
    histref.Draw( "hist" )
    histres= TH1D( "histres", "1000*x*exp(-x*20) avg. result", nbin, binEdges )
    SetOwnership( histres, False )
    for i in range( nbin ):
        histres.SetBinContent( i+1, mvalues[i] )
        histres.SetBinError( i+1, sqrt( errorMatrixNormAvg[i][i] ) )
    histres.SetMarkerStyle( 20 )
    histres.SetMarkerSize( 0.5 )
    histres.Draw( "samepE1" )
    canv1.cd( 6 )
    histpulls= TH1D( "histpulls", "pulls corrected - reference ", nbin, binEdges )
    SetOwnership( histpulls, False )
    for i in range( nbin ):
        histpulls.SetBinContent( i+1, pulls[i] )
        histpulls.SetBinError( i+1, pullsstd[i] )
    histpulls.SetMarkerStyle( 20 )
    histpulls.SetMarkerSize( 0.5 )
    histpulls.Draw( "pE1" )
    
    # Fin
    return

# Sample by generating from reference histos
def fillSampleHistosDraw( histref, histdet, histdet2, nsample, mnevent, title ):
    histds= dict()
    histd2s= dict()
    hisths= dict()
    for isample in range( nsample ):
        # title= fun.GetTitle()
        histh= TH1D( "histh"+str(isample), title+" h", nbin, binEdges )
        histd= TH1D( "histd"+str(isample), title+" d", nbin, binEdges )
        histd2= TH1D( "histd2"+str(isample), title+" d 2", nbin, binEdges )
        nevent= gRandom.Poisson( mnevent )
        for ievent in range( nevent ):
            histh.Fill( histref.GetRandom() )
            histd.Fill( histdet.GetRandom() )
            histd2.Fill( histdet2.GetRandom() )
        hisths[isample]= histh
        histds[isample]= histd
        histd2s[isample]= histd2
    return histds, histd2s, hisths
        
# Sample by generating from "true" function and smearing
def fillSampleHistosSmear( fun, nsample, mnevent, resolutionFactor, opt ):
    histds= dict()
    histd2s= dict()
    hisths= dict()
    for isample in range( nsample ):
        title= fun.GetTitle()
        histh= TH1D( "histh"+str(isample), title+" h", nbin, binEdges )
        histd= TH1D( "histd"+str(isample), title+" d", nbin, binEdges )
        histd2= TH1D( "histd2"+str(isample), title+" d 2", nbin, binEdges )
        nevent= gRandom.Poisson( mnevent )
        neventsFilled= 0.0
        for i in range( nevent ):
            nmom= nmomFromOpt( opt )
            value= fun.GetRandom()
            hWeight= value**nmom
            histh.Fill( value, hWeight )
            resValue= resolution( value, resolutionFactor )
            dWeight= resValue**nmom
            dWeight2= dWeight**2
            histd.Fill( resValue, dWeight )
            histd2.Fill( resValue, dWeight2 )
            if inBinRange( resValue ):
                neventsFilled+= 1.0
        hisths[isample]= histh
        histd.SetEntries( neventsFilled )
        histds[isample]= histd
        histd2.SetEntries( neventsFilled )
        histd2s[isample]= histd2
    return histds, histd2s, hisths

def nmomFromOpt( opt ):
    import re
    mo= re.search( "w[1-5]", opt )
    result= 0
    if mo:
        result= int( mo.group()[1] )
    return result
    
# Average matrix from collection
def averageMatrix( errorMatrices ):
    merrorMatrix= np.zeros( (nbin,nbin) )
    nsample= len( errorMatrices )
    for isample in range( nsample ):
        merrorMatrix+= errorMatrices[isample]
    merrorMatrix/= float( nsample )
    return merrorMatrix

# Sample error matrix from collection of histograms
def sampleErrorMatrix( hists ):
    merrorMatrix= np.zeros( (nbin,nbin) )
    mvalues= np.array( nbin*[0.0] )
    nsample= len( hists )
    for isample in range( nsample ):
        hist= hists[isample]
        for i in range( nbin ):
            valuei= hist.GetBinContent( i+1 )
            mvalues[i]+= valuei
            for j in range( nbin ):
                 valuej= hist.GetBinContent( j+1 )
                 merrorMatrix[i,j]+= valuei*valuej
    mvalues/= float( nsample )
    merrorMatrix/= float( nsample )
    for i in range( nbin ):
            for j in range( nbin ):
                merrorMatrix[i,j]-= mvalues[i]*mvalues[j]
    return merrorMatrix

# Correct using reverse migration matrix incl. efficiency
def matrixCorrection( histds, cf, CR, eventAcc, histdet, opt ):
    histcs= dict()
    errorMatrices= dict()
    nsample= len( histds )
    for isample in range( nsample ):
        histd= histds[isample]
        title= histd.GetTitle()
        histc= TH1D( "histc"+str(isample), title+" corr", nbin, binEdges )
        SetOwnership( histc, False )
        valuesd= getArray( histd )
        errorsd= getArrayErrors( histd )
        valuesh= CR.dot( valuesd )
        if "u" in opt:
            # Use expected error
            scf= histd.GetEntries()/histdet.GetEntries()
            for i in range( nbin ):
                errorsd[i]= histdet.GetBinError( i+1 )*sqrt( scf )        
        # Error propagation from correction matrix; it is not
        # symmetric so order of indices matters
        diagErrorMatrixd= np.diag( errorsd**2 )
        errorMatrixh= np.dot( CR, np.dot( diagErrorMatrixd, CR.transpose() ) )
        # Set corrected values and bin-by-bin error in histo
        errorsh= errorsd*cf
        for i in range( nbin ):
            histc.SetBinContent( i+1, valuesh[i] )
            histc.SetBinError( i+1, errorsh[i] )
        # Number of entries from filled events and event acceptance
        nevent= histd.GetEntries()
        neventCorrected= nevent*eventAcc
        histc.SetEntries( neventCorrected )
        histcs[isample]= histc
        errorMatrices[isample]= errorMatrixh
    # Fin
    return histcs, errorMatrices

# Bin-by-bin correction
def bbbCorrection( histds, cf, eventAcc, histdet, opt ):
    histcs= dict()
    errorMatrices= dict()
    nsample= len( histds )
    for isample in range( nsample ):
        histd= histds[isample]
        title= histd.GetTitle()
        histc= TH1D( "histc", title+" corr", nbin, binEdges )
        valuesh= np.array( nbin*[0.0] )
        errorsh= np.array( nbin*[0.0] )
        for i in range( nbin ):
            valueh= histd.GetBinContent( i+1 )*cf[i]
            if "u" in opt:
                # Expected error
                scf= histd.GetEntries()/histdet.GetEntries()
                errorh= histdet.GetBinError( i+1 )*sqrt( scf )
            else:
                errorh= histd.GetBinError( i+1 )
            errorh*= cf[i]
            histc.SetBinContent( i+1, valueh )
            histc.SetBinError( i+1, errorh )
            valuesh[i]= valueh
            errorsh[i]= errorh
        nevent= histd.GetEntries()
        neventCorrected= nevent*eventAcc
        histc.SetEntries( neventCorrected )
        histcs[isample]= histc
        errorMatrices[isample]= np.diag( errorsh**2 )
    return histcs, errorMatrices

# Calculate chi^2 with error matrix
def chi2Matrix( errorMatrix, histc, histt, first=0, last=nbin-1 ):
    from scipy import linalg
    errorMatrixInv= linalg.inv( errorMatrix[first:last,first:last] )
    chi2= 0.0
    for i in range( last-first ):
        for j in range( last-first ):
            chi2+= ( errorMatrixInv[i][j]*
                         ( histc.GetBinContent( first+i+1 ) -
                               histt.GetBinContent( first+i+1 ) )*
                         ( histc.GetBinContent( first+j+1 ) -
                               histt.GetBinContent( first+j+1 ) ) )
    return chi2

# Print a matrix formatted
def printMatrix( m, width=7, precision=3 ):
    fmt="{:"+str(width)+"."+str(precision)+"f}"
    ( imax, jmax )= m.shape
    for i in range( imax ):
        for j in range( jmax ):
            print fmt.format( m[i,j] ),
        print
    return

# Convert covariance to correlation
def cov2corr( m ):
    corr= np.zeros( (nbin,nbin) )
    for i in range( nbin ):
        for j in range( nbin ):
            corr[i,j]= m[i,j]/sqrt( m[i,i]*m[j,j] )
    return corr

# Error matrix for normalised histo from diagonal errors
def calcErrorMatrix( values, errors, nevent ):
    errorMatrix= np.zeros( (nbin,nbin) )
    sumw= sum( values )
    for i in range( nbin ):
        for j in range( nbin ):
            cov= 0.0
            for k in range( nbin ):
                deltaik= 0.0
                if i == k:
                    deltaik= 1.0
                deltajk= 0.0
                if j == k:
                    deltajk= 1.0
                cov+= ( errors[k]**2*
                            ( sumw*deltaik - values[i] )*
                            ( sumw*deltajk - values[j] ) )
            cov/= (float(nevent)*sumw)**2
            errorMatrix[i,j]= cov
    return errorMatrix

# Generate correction matrix and factors
def generate( fun, n=1000, acc=0.8, resolutionFactor=0.01, opt="" ):

    print "Generate events for corrections"
    
    # Histos at hadron and detector level, migrations
    title= fun.GetTitle()
    histh= TH1D( "histh", title+" h", nbin, binEdges )
    histd= TH1D( "histd", title+" d", nbin, binEdges )
    histmd= TH2D( "histmd", "migration h vs d d", nbin, binEdges, nbin, binEdges )
    histmh= TH2D( "histmh", "migration h vs d h", nbin, binEdges, nbin, binEdges )
    histh2= TH1D( "histh2", title+" h", nbin, binEdges )
    histd2= TH1D( "histd2", title+" d", nbin, binEdges )
    histmd2= TH2D( "histmd2", "migration h vs d d", nbin, binEdges, nbin, binEdges )
    histmh2= TH2D( "histmh2", "migration h vs d h", nbin, binEdges, nbin, binEdges )
    neventsd= 0.0
    for i in range( n ):
        nmom= nmomFromOpt( opt )
        hValue= fun.GetRandom()
        hWeight= hValue**nmom
        hWeight2= hWeight**2
        histh.Fill( hValue, hWeight )
        histh2.Fill( hValue, hWeight2 )
        if acceptance( acc ):
            dValue= resolution( hValue, resolutionFactor )
            dWeight= dValue**nmom
            dWeight2= dWeight**2
            histd.Fill( dValue, dWeight )
            histmd.Fill( dValue, hValue, dWeight )
            histmh.Fill( dValue, hValue, hWeight )
            histd2.Fill( dValue, dWeight2 )
            histmd2.Fill( dValue, hValue, dWeight2 )
            histmh2.Fill( dValue, hValue, hWeight2 )
            if inBinRange( dValue ):
                neventsd+= 1.0

    # Prepare for expected errors calculation, exclude
    # over- and underflows
    histd.SetEntries( neventsd )
    histd2.SetEntries( neventsd )

    # Event acceptance
    eventAcc= histh.GetEntries()/neventsd
    print "Generated events before/after cuts:", histh.GetEntries(), neventsd
    print "Event acceptance:", eventAcc
    
    # Bin-by-bin correction factors
    cf= divideHistos( histh, histd )
    cf2= divideHistos( histh2, histd2 )

    # Migration matrix normalised column-wise as correction
    revMigr= getMatrix( histmd )
    revMigrNorm= normaliseColumns( revMigr )
    revMigr2= getMatrix( histmd2 )
    revMigrNorm2= normaliseColumns( revMigr2 )

    # h->d correction of migration matrix if weighted
    if "w" in opt:
        revMigrh= getMatrix( histmh )
        hdRatioMatrix= divideChecked( revMigrh, revMigr )
        revMigrNormCorr= hdRatioMatrix*revMigrNorm
        revMigrh2= getMatrix( histmh2 )
        hdRatioMatrix2= divideChecked( revMigrh2, revMigr2 )
        revMigrNormCorr2= hdRatioMatrix2*revMigrNorm2
    else:
        revMigrNormCorr= revMigrNorm
        revMigrNormCorr2= revMigrNorm2

    # Acceptance correction
    histmhy= histmh.ProjectionY( "revMigr y proj", 1, nbin )
    cacc= divideHistos( histh, histmhy )
    histmhy2= histmh2.ProjectionY( "revMigr2 y proj", 1, nbin )
    cacc2= divideHistos( histh2, histmhy2 )

    # Combine matrix and acceptance correction
    CR= np.dot( np.diag( cacc ), revMigrNormCorr )
    CR2= np.dot( np.diag( cacc2 ), revMigrNormCorr2 )

    # Closure tests
    valuesd= getArray( histd )
    valuesh= getArray( histh )
    print "Rev. migr., cacc"
    printMatrix( revMigr )
    print cacc
    print "Closure test hadron vs corrected"
    print valuesh
    print CR.dot( valuesd )

    valuesd2= getArray( histd2 )
    valuesh2= getArray( histh2 )
    print "Closure test hadron vs corrected mom**2"
    print valuesh2
    print CR2.dot( valuesd2 )
    
    # Fin
    return cf, CR, cf2, CR2, eventAcc, histh, histd, histh2, histd2

# Normalise matrix columns
def normaliseColumns( matrix ):
    shape= matrix.shape
    result= np.array( np.zeros( shape ) )
    for icol in range( shape[1] ):
        sumcol= sum( matrix[:,icol] )
        if sumcol != 0.0:
            result[:,icol]= matrix[:,icol]/sumcol
    return result

# Normalise matrix rows
def normaliseRows( matrix ):
    shape= matrix.shape
    result= np.array( np.zeros( shape ) )
    for irow in range( shape[0] ):
        sumrow= sum( matrix[irow,:] )
        if sumrow != 0.0:
            result[irow,:]= matrix[irow,:]/sumrow
    return result

# Normalise histogram by number of entries and calculate
# covariance matrix
def normalise( hist, errorMatrix=None ):
    nevent= hist.GetEntries()
    values= getArray( hist )
    errors= getArrayErrors( hist )
    if errorMatrix is None:
        errorMatrix= np.diag( errors**2 )
    errorMatrixNorm= calcErrorMatrixJacobean( values, errorMatrix, nevent )
    for i in range( nbin ):
        binwi= hist.GetBinWidth( i+1 )
        for j in range( nbin ):
            binwj= hist.GetBinWidth( j+1 )
            errorMatrixNorm[i,j]/= (binwi*binwj)
    for i in range( 1, hist.GetNbinsX()+1 ):
        binw= hist.GetBinWidth( i )
        value= hist.GetBinContent( i )/( nevent*binw )
        error= hist.GetBinError( i )/( nevent*binw )
        hist.SetBinContent( i, value )
        hist.SetBinError( i, error )
    return errorMatrixNorm

# Transform error matrix with Jacobean from normalisation
def calcErrorMatrixJacobean( values, errorMatrix, nevent ):
    sumw= sum( values )
    J= np.zeros( (nbin,nbin) )
    for i in range( nbin ):
        for j in range( nbin ):
            J[j,i]= ( sumw*delta( i, j ) - values[i] )/(nevent*sumw)
    errorMatrixNormalised= np.dot( J.transpose(), np.dot( errorMatrix, J ) )   
    return errorMatrixNormalised

# Kronecker delta function
def delta( indx1, indx2 ):
    if indx1 == indx2:
        return 1.0
    else:
        return 0.0

# Apply resolution effect, values may be negative or too large
# and end up in under- or overflow
def resolution( value, binw ):
    delta= gRandom.Gaus( 0.0, 1.0 )
    deltaValue= value + delta*binw/2.0
    # if deltaValue <= 0.0:
    #     deltaValue= 0.0001
    # if deltaValue > 0.5:
    #     deltaValue= 0.4999
    return deltaValue

# Simulate limited acceptance
def acceptance( acc=0.8 ):
    if gRandom.Rndm() < acc:
        return True
    else:
        return False

    
