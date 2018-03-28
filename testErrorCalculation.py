
from ROOT import TH1D, TH2D, TF1, TMath, TCanvas, SetOwnership, TRandom3, gRandom, gStyle

TCanvas.__init__._creates= False

import numpy as np

from math import sqrt

binEdges= np.array( [ 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.09, 0.12, 0.15, 0.22, 0.30, 0.50 ] )
nbin= len(binEdges)-1

def inBinRange( value ):
    return binEdges[0] < value and value < binEdges[nbin]

def testWeighted( mnevent=1000, nsample=1000, resolutionFactor=0.02 ):

    # Prepare
    gRandom= TRandom3()
    fun= TF1( "fun", "1000*x*exp(-x*20)", 0.0 , 0.5 )
    histh= TH1D( "histh", "1000*x*exp(-x*20) h", nbin, binEdges )
    histhw1= TH1D( "histhw1", "1000*x*exp(-x*20) w1 h", nbin, binEdges )
    histd= TH1D( "histd", "1000*x*exp(-x*20) d", nbin, binEdges )
    histdw1= TH1D( "histdw1", "1000*x*exp(-x*20) w1 d", nbin, binEdges )
    histm= TH2D( "histm", "migration h vs d", nbin, binEdges, nbin, binEdges )
    histmw1h= TH2D( "histmw1h", "migration h vs d w1 h", nbin, binEdges, nbin, binEdges )
    histmw1d= TH2D( "histmw1d", "migration h vs d w1 d", nbin, binEdges, nbin, binEdges )
    SetOwnership( histh, False )
    SetOwnership( histhw1, False )
    SetOwnership( histd, False )
    SetOwnership( histdw1, False )
    SetOwnership( histm, False )
    SetOwnership( histmw1h, False )
    SetOwnership( histmw1d, False )

    # Generate
    for i in range( mnevent ):
        value= fun.GetRandom()
        histh.Fill( value )
        resValue= resolution( value, resolutionFactor )
        histd.Fill( resValue )
        histhw1.Fill( value, value )
        histdw1.Fill( resValue, resValue )
        histm.Fill( resValue, value )
        histmw1h.Fill( resValue, value, value )
        histmw1d.Fill( resValue, value, resValue )
        
    # Print raw migration matrices
    print "Migration"
    migrMatrix= getMatrix( histm )
    valuesh= getArray( histh )
    valuesd= getArray( histd )
    printMatrixVectors( migrMatrix, valuesd, valuesh )
    print "Migration w1 h"
    migrMatrixw1h= getMatrix( histmw1h )
    valueshw1= getArray( histhw1 )
    valuesdw1= getArray( histdw1 )
    printMatrixVectors( migrMatrixw1h, valuesdw1, valueshw1 )
    print "Migration w1 d"
    migrMatrixw1d= getMatrix( histmw1d )
    printMatrixVectors( migrMatrixw1d, valuesdw1, valueshw1 )
    
    # Normalised det-level migration
    migrMatrixw1dNorm= normaliseColumns( migrMatrixw1d )
    print "Migration w1 d normalised column-wise"
    printMatrix( migrMatrixw1dNorm )

    # Ratio of had- to det-level migrations
    hdratioMatrix= divideChecked( migrMatrixw1h, migrMatrixw1d )
    print "Migration w1 h/d ratio"
    printMatrix( hdratioMatrix )

    # Multiply normalised det-level migration by had- to det-level ratio
    migrMatrixw1Rec= hdratioMatrix*migrMatrixw1dNorm
    print "Migration w1 d->h corr. normalised column-wise"
    printMatrix( migrMatrixw1Rec )
    
    # Efficiency correction
    histmw1y= histmw1h.ProjectionY( "Migr w1 h y proj", 1, nbin )
    cacc= divideHistos( histhw1, histmw1y )
    print "Effciciency h", cacc

    # Final correction matrix and closure test
    CR= np.dot( np.diag( cacc ), migrMatrixw1Rec )
    valueshw1Rec= np.dot( CR, valuesdw1 )
    print "Closure"
    print valueshw1Rec
    print valueshw1

    # Put results with propagated errors into histo for normalisation
    histhrec= TH1D( "histhrec", "1000*x*exp(-x*20) h rec", nbin, binEdges )
    for i in range( nbin ):
        histhrec.SetBinContent( i+1, valueshw1Rec[i] )
    histhrec.SetEntries( mnevent )
    errorsdw1= getArrayErrors( histdw1 )
    diagErrorMatrixdw1= np.diag( errorsdw1**2 )
    errorMatrixhw1= np.dot( CR, np.dot( diagErrorMatrixdw1, CR.transpose() ) )
    errorMatrixhw1Norm= normalise( histhrec, errorMatrixhw1 )

    # Print results
    width, precision= 7, 4
    fmt="{:"+str(width)+"."+str(precision)+"f}"
    print "Result corrected normalised w1 distribution"
    for i in range( nbin ):
        valuehw1RecNorm= histhrec.GetBinContent( i+1 )
        errorhw1RecNorm= sqrt( errorMatrixhw1Norm[i,i] )
        histhrec.SetBinError( i+1, errorhw1RecNorm )
        print fmt.format( valuehw1RecNorm ), fmt.format( errorhw1RecNorm )        
    print "Correlation matrix"
    corr= cov2corr( errorMatrixhw1Norm )
    printMatrix( corr )

    # Integral of corrected weighted is 1st moment of unweighted
    # Error from diagonal bin errors only
    sumw1= calcIntegral( histhrec )
    errors= getArrayErrors( histhrec )
    binw= binEdges[1:]-binEdges[:nbin]
    errorsbw= errors*binw
    errSumw1= sqrt( sum(errorsbw**2) )

    # Error of integral given by Sum V_ij*binw_i*binw_j
    # but does not work in single closure sample
    # errSumw1= 0.0
    # for i in range( nbin ):
    #     for j in range( nbin ):
    #         term= corr[i,j]*errorsbw[i]*errorsbw[j]
    #         errSumw1+= term
    #         print corr[i,j], errorsbw[i], errorsbw[j], term, errSumw1
    # errSumw1= np.sum( errorMatrixhw1Norm )
    # errSumw1= sqrt( errSumw1 )
    width, precision= 7, 4
    fmt="{:"+str(width)+"."+str(precision)+"f}"
    print "Integrated corrected w1 distribution with diag. errors:",
    print fmt.format( sumw1 ), "+/-", fmt.format( errSumw1 )
    print "Mean of had-level histo from ROOT:",
    print fmt.format( histh.GetMean() ), "+/-", fmt.format( histh.GetMeanError() )

    # Plots
    canv= TCanvas( "canv", "Weighted correction", 800, 1000 )
    canv.Divide( 2, 3 )
    canv.cd( 1 )
    normalise( histhw1 )
    histhw1.Draw()
    canv.cd( 2 )
    normalise( histdw1 )
    histdw1.Draw()
    canv.cd( 3 )
    histm.Draw( "box" )
    canv.cd( 4 )
    histmw1h.Draw( "box" )
    canv.cd( 5 )
    histmw1d.Draw( "box" )
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
        raise "shapes not equal"
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
            fillweight= weight( value, opt )
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

    # Initialise ROOT random numbers
    gRandom= TRandom3()

    # Function to sample
    if "e" in opt:
        fun= TF1( "fun", "1000*x*exp(-x*20)", 0.0 , 0.5 )
    elif "g" in opt:
        fun= TF1( "fun", "TMath::Gaus(x,0.25,0.15)", 0.0 , 0.5 )
    fun.SetNormalized(True)
    
    # Correction factors, matrix and reference histo
    cf, CR, eventAcc, histref, histdet= generate( fun, 1000000, acc, resolutionFactor, opt )

    # Create and fill sample histos
    if "s" in opt:
        histds, hisths= fillSampleHistosSmear( fun, nsample, mnevent,
                                                            resolutionFactor, opt )
    elif "d" in opt:
        histds, hisths= fillSampleHistosDraw( histref, histdet, nsample,
                                                           mnevent )
        
    # Error matrices after filling
    errorMatrixhSample= sampleErrorMatrix( hisths )
    errorMatrixdSample= sampleErrorMatrix( histds )
    
    # Get results with migration matrix or bin-by-bin correction
    if "b" in opt:
        histcs, errorMatricesCorr= bbbCorrection( histds, cf, eventAcc,
                                                      histdet, opt )
    elif "m" in opt:
        histcs, errorMatricesCorr= matrixCorrection( histds, cf, CR, eventAcc,
                                                         histdet, opt )

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
    for isample in range( nsample ):
        histc= histcs[isample]
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
    mvaluesh/= float(nsample)
    merrorsh/= float(nsample)
    mvalues/= float(nsample)
    merrors/= float(nsample)
    merrors-= mvalues**2
    merrors= np.sqrt( merrors )
    merrorsh-= mvaluesh**2
    merrorsh= np.sqrt( merrorsh )
    merrors2/= float(nsample)

    errorMatrixNormSample= sampleErrorMatrix( histcs )
    errorMatrixNormAvg= averageMatrix( errorMatricesNorm )
    
    # Reference from histref
    normalise( histref )
        
    # Pulls
    pulls= np.array( nbin*[0.0] )
    pullsstd= np.array( nbin*[0.0] )    
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
    pulls/= float(nsample)
    pullsstd/= float(nsample)
    pullsstd-= pulls**2
    pullsstd= np.sqrt( pullsstd )
            
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
def fillSampleHistosDraw( histref, histdet, nsample, mnevent ):
    histds= dict()
    hisths= dict()
    for isample in range( nsample ):
        histh= TH1D( "histh"+str(isample), "1000*x*exp(-x*20) h", nbin, binEdges )
        histd= TH1D( "histd"+str(isample), "1000*x*exp(-x*20) d", nbin, binEdges )
        nevent= gRandom.Poisson( mnevent )
        for ievent in range( nevent ):
            histh.Fill( histref.GetRandom() )
            histd.Fill( histdet.GetRandom() )
        hisths[isample]= histh
        histds[isample]= histd
    return histds, hisths
        
# Sample by generating from "true" function and smearing
def fillSampleHistosSmear( fun, nsample, mnevent, resolutionFactor, opt ):
    histds= dict()
    hisths= dict()
    for isample in range( nsample ):
        histh= TH1D( "histh"+str(isample), "1000*x*exp(-x*20) h", nbin, binEdges )
        histd= TH1D( "histd"+str(isample), "1000*x*exp(-x*20) d", nbin, binEdges )
        nevent= gRandom.Poisson( mnevent )
        neventsFilled= 0.0
        for i in range( nevent ):
            value= fun.GetRandom()
            hWeight= weight( value, opt )
            histh.Fill( value, hWeight )
            resValue= resolution( value, resolutionFactor )
            dWeight= weight( resValue, opt )
            histd.Fill( resValue, dWeight )
            if inBinRange( resValue ):
                neventsFilled+= 1.0
        hisths[isample]= histh
        histd.SetEntries( neventsFilled )
        histds[isample]= histd
    return histds, hisths

# Weight or not
def weight( value, opt ):
    if "w1" in opt:
        return value
    elif "w2" in opt:
        return value**2
    elif "w3" in opt:
        return value**3
    elif "w4" in opt:
        return value**4
    elif "w5" in opt:
        return value**5
    else:
        return 1.0

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
        histc= TH1D( "histc"+str(isample), "1000*x*exp(-x*20) corr", nbin, binEdges )
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
        # errorMatrixh= np.zeros( (nbin,nbin) )
        # for i in range( nbin ):
        #     for j in range( nbin ):
        #         sumk= 0.0
        #         for k in range( nbin ):
        #             sumk+= errorsd[k]**2*CR[i,k]*CR[j,k]
        #         errorMatrixh[i,j]= sumk                
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
        histc= TH1D( "histc", "1000*x*exp(-x*20) corr", nbin, binEdges )
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
            corr[i,j]= m[i,j]/sqrt( m[i,i]* m[j,j])
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
    
    # Histos at hadron and detector level, migrations
    histh= TH1D( "histh", "1000*x*exp(-x*20) h", nbin, binEdges )
    histd= TH1D( "histd", "1000*x*exp(-x*20) d", nbin, binEdges )
    histm= TH2D( "histm", "migration h vs d d", nbin, binEdges, nbin, binEdges )
    histmh= TH2D( "histm", "migration h vs d h", nbin, binEdges, nbin, binEdges )
    neventsd= 0.0
    for i in range( n ):
        hValue= fun.GetRandom()
        hWeight= weight( hValue, opt )
        histh.Fill( hValue, hWeight )
        if acceptance( acc ):
            dValue= resolution( hValue, resolutionFactor )
            dWeight= weight( dValue, opt )
            histd.Fill( dValue, dWeight )
            histm.Fill( dValue, hValue, dWeight )
            histmh.Fill( dValue, hValue, hWeight )
            if inBinRange( dValue ):
                neventsd+= 1.0

    # Prepare for expected errors calculation, exclude
    # over- and underflows
    histd.SetEntries( neventsd )

    # Event acceptance
    eventAcc= histh.GetEntries()/neventsd
    print "Events: ", histh.GetEntries(), neventsd
    print "Event acceptance:", eventAcc
    
    # Bin-by-bin correction factors
    cf= divideHistos( histh, histd )
        
    # Migration matrix normalised row-wise as correction
    revMigr= getMatrix( histm )
    revMigrNorm= normaliseColumns( revMigr )

    # h->d correction of migration matrix if weighted
    if "w" in opt:
        revMigrh= getMatrix( histmh )
        hdRatioMatrix= divideChecked( revMigrh, revMigr )
        revMigrNormCorr= hdRatioMatrix*revMigrNorm
    else:
        revMigrNormCorr= revMigrNorm

    # Acceptance correction
    histmhy= histmh.ProjectionY( "revMigr y proj", 1, nbin )
    cacc= divideHistos( histh, histmhy )

    # Combine matrix and acceptance correction
    C= np.diag( cacc )
    CR= np.dot( C, revMigrNormCorr )

    # Closure test
    valuesd= getArray( histd )
    valuesh= getArray( histh )
    print "Rev. migr., cacc"
    printMatrix( revMigr )
    print cacc
    print "Closure"
    print valuesh
    print CR.dot( valuesd )
    
    # Fin
    return cf, CR, eventAcc, histh, histd

# Normalise matrix columns
def normaliseColumns( matrix ):
    shape= matrix.shape
    result= np.array( np.zeros( shape ) )
    for i in range( shape[1] ):
        result[:,i]= matrix[:,i]/sum( matrix[:,i] )
    return result

# Normalise histogram to unit area
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

    
