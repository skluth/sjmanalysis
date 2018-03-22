
from ROOT import TH1D, TH2D, TF1, TMath, TCanvas, SetOwnership, TRandom3, gRandom

TCanvas.__init__._creates= False

import numpy as np

from math import sqrt

binEdges= np.array( [ 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.09, 0.12, 0.15, 0.22, 0.30, 0.50 ] )
nbin= len(binEdges)-1

# Test sampling
def testSampling( mnevent=1000, nsample=1000, opt="p" ):
    if "p" in opt:
        print "Variable number of events from Poisson mu=", mnevent
    else:
        print "Fixed number of events n=", mnevent
    fun= TF1( "fun", "1000*x*exp(-x*20)", 0.0 , 0.5 )
    gRandom= TRandom3()
    hists= dict()
    for isample in range( nsample ):
        hist= TH1D( "histh"+str(isample), "1000*x*exp(-x*20) h", nbin, binEdges )
        if "p" in opt:
            nevent= gRandom.Poisson( mnevent )
        else:
            nevent= mnevent
        for i in range( int( nevent ) ):
            value= fun.GetRandom()
            hist.Fill( value )
        hists[isample]= hist
    errorMatrixSample= sampleErrorMatrix( hists )
    printMatrix( errorMatrixSample, 7, 3 )
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
    cf, CR, eventAcc, histref, histdet= generate( fun, 1000000, acc, resolutionFactor )

    # Create and fill sample histos
    if "s" in opt:
        histds, hisths, nevents= fillSampleHistosSmear( fun, nsample, mnevent,
                                                            resolutionFactor )
    elif "d" in opt:
        histds, hisths, nevents= fillSampleHistosDraw( histref, histdet, nsample, mnevent )
        
    # Error matrices after filling
    errorMatrixdSample= sampleErrorMatrix( histds )
    errorMatrixhSample= sampleErrorMatrix( hisths )
    
    # Get results with migration matrix or bin-by-bin correction
    if "b" in opt:
        histcs, errorMatricesCorr= bbbCorrection( histds, cf, eventAcc,
                                                      nevents, histdet, opt )
    elif "m" in opt:
        histcs, errorMatricesCorr= matrixCorrection( histds, cf, CR, eventAcc,
                                                         nevents, histdet, opt )

    # Error matrices after correction
    errorMatrixCorrAvg= averageMatrix( errorMatricesCorr )
    errorMatrixCorrSample= sampleErrorMatrix( histcs )

    print "Error matrix hadron"
    printMatrix( errorMatrixhSample, 7, 3 )
    print "Error matrix detector"
    printMatrix( errorMatrixdSample, 7, 3 )
    print "Error matrix after correction calculated avg."
    printMatrix( errorMatrixCorrAvg, 7, 3 )
    print "Error matrix after correction sampled"
    printMatrix( errorMatrixCorrSample, 7, 3 )

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
        errorMatrix= errorMatricesNorm[isample]
        # errorMatrix= errorMatrixNormAvg
        # errorMatrix= merrMatrixSample
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
    canv= TCanvas( "canv", "Chi^2", 800, 800 )
    canv.Divide( 2, 2 )
    canv.cd( 1 )
    hchi2diag.Draw()
    canv.cd( 2 )
    hchi2matrix.Draw()
    canv.cd( 3 )
    hchi2matrixs.Draw()
    canv.cd( 4 )
    hchi2matrixa.Draw()

    # Fin
    return

# Sample by generating from reference histos
def fillSampleHistosDraw( histref, histdet, nsample, mnevent ):
    histds= dict()
    hisths= dict()
    nevents= dict()
    for isample in range( nsample ):
        histh= TH1D( "histh"+str(isample), "1000*x*exp(-x*20) h", nbin, binEdges )
        histd= TH1D( "histd"+str(isample), "1000*x*exp(-x*20) d", nbin, binEdges )
        nevent= gRandom.Poisson( mnevent )
        for ievent in range( nevent ):
            histh.Fill( histref.GetRandom() )
            histd.Fill( histdet.GetRandom() )
        hisths[isample]= histh
        histds[isample]= histd
        nevents[isample]= nevent
    return histds, hisths, nevents
        
# Sample by generating from "true" function and smearing
def fillSampleHistosSmear( fun, nsample, mnevent, resolutionFactor ):
    histds= dict()
    hisths= dict()
    nevents= dict()
    for isample in range( nsample ):
        histh= TH1D( "histh"+str(isample), "1000*x*exp(-x*20) h", nbin, binEdges )
        histd= TH1D( "histd"+str(isample), "1000*x*exp(-x*20) d", nbin, binEdges )
        nevent= gRandom.Poisson( mnevent )
        for i in range( nevent ):
            value= fun.GetRandom()
            histh.Fill( value )
            resValue= resolution( value, resolutionFactor )
            histd.Fill( resValue )
        hisths[isample]= histh
        histds[isample]= histd
        nevents[isample]= nevent
    return histds, hisths, nevents

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
def matrixCorrection( histds, cf, CR, eventAcc, nevents, histdet, opt ):
    histcs= dict()
    errorMatrices= dict()
    nsample= len( histds )
    for isample in range( nsample ):
        histd= histds[isample]
        histc= TH1D( "histc"+str(isample), "1000*x*exp(-x*20) corr", nbin, binEdges )
        valuesd= np.array( nbin*[0.0] )
        errorsd= np.array( nbin*[0.0] )
        for i in range( nbin ):
            valuesd[i]= histd.GetBinContent( i+1 )
            errorsd[i]= histd.GetBinError( i+1 )
            if valuesd[i] < 10.0:
                print "bin is low", i, valuesd[i], isample
        valuesh= CR.dot( valuesd )
        # Use expected error
        if "u" in opt:
            scf= float(nevents[isample])/histdet.Integral()
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
        # nevent= nevents[isample]
        # neventCorrected= nevent*eventAcc
        # histc.SetEntries( neventCorrected )        
        histc.SetEntries( sum( valuesh ) )
        histcs[isample]= histc
        errorMatrices[isample]= errorMatrixh
    # Fin
    return histcs, errorMatrices

# Bin-by-bin correction
def bbbCorrection( histds, cf, eventAcc, nevents, histdet, opt ):
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
                scf= float(nevents[isample])/histdet.Integral()
                errorh= histdet.GetBinError( i+1 )*sqrt( scf )
            else:
                errorh= histd.GetBinError( i+1 )
            errorh*= cf[i]
            histc.SetBinContent( i+1, valueh )
            histc.SetBinError( i+1, errorh )
            valuesh[i]= valueh
            errorsh[i]= errorh
        # nevent= nevents[isample]
        # neventCorrected= nevent*eventAcc
        # histc.SetEntries( neventCorrected )
        histc.SetEntries( sum( valuesh ) )
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
def printMatrix( m, width=6, precision=3 ):
    fmt="{:"+str(width)+"."+str(precision)+"f}"
    ( imax, jmax )= m.shape
    for i in range( imax ):
        for j in range( jmax ):
            print fmt.format( m[j][i] ),
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
def generate( fun, n=1000, acc=0.8, resolutionFactor=0.01 ):
    
    # Histos at hadron and detector level, migrations
    histh= TH1D( "histh", "1000*x*exp(-x*20) h", nbin, binEdges )
    histd= TH1D( "histd", "1000*x*exp(-x*20) d", nbin, binEdges )
    histm= TH2D( "histm", "migration h vs d", nbin, binEdges, nbin, binEdges )
    for i in range( n ):
        hValue= fun.GetRandom()
        histh.Fill( hValue )
        if acceptance( acc ):
            dValue= resolution( hValue, resolutionFactor )
            histd.Fill( dValue )
            histm.Fill( dValue, hValue )

    # Event acceptance
    eventAcc= histh.GetEntries()/histd.Integral()
    print "Events: ", histh.GetEntries(), histd.Integral()
    print "Event acceptance:", eventAcc
    
    # Bin-by-bin correction factors
    cf= np.array( nbin*[0.0] )
    for i in range( nbin ):
        cf[i]= histh.GetBinContent( i+1 ) / histd.GetBinContent( i+1 )
        
    # Migration matrix normalised row-wise as correction
    revMigr= np.zeros( ( nbin, nbin ) )
    for i in range( nbin ):
        for j in range( nbin ):
            revMigr[j,i]= histm.GetBinContent( i+1, j+1 )            
    for i in range( nbin ):
        revMigr[:,i]/= sum( revMigr[:,i] )
        
    # Acceptance correction
    histmx= histm.ProjectionY( "revMigr y proj", 1, nbin )
    cacc= np.array( nbin*[0.0] )
    for i in range( nbin ):
        cacc[i]= histh.GetBinContent( i+1 ) / histmx.GetBinContent( i+1 )

    # Combine matrix and acceptance
    C= np.diag( cacc )
    CR= np.dot( C, revMigr )

    # Closure test
    valuesd= np.array( nbin*[0.0] )
    valuesh= np.array( nbin*[0.0] )
    for i in range( nbin ):
        valuesd[i]= histd.GetBinContent( i+1 )
        valuesh[i]= histh.GetBinContent( i+1 )
    print "Rev. migr., cacc"
    printMatrix( revMigr )
    print cacc
    print "Closure"
    print valuesh
    print CR.dot( valuesd )
    
    # Fin
    return cf, CR, eventAcc, histh, histd

# Normalise histogram to unit area
def normalise( hist, errorMatrix=None ):
    nevent= hist.GetEntries()
    values= np.array( nbin*[0.0] )
    errors= np.array( nbin*[0.0] )
    for i in range( nbin ):
        values[i]= hist.GetBinContent( i+1 )
        errors[i]= hist.GetBinError( i+1 )
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
    return deltaValue

# Simulate limited acceptance
def acceptance( acc=0.8 ):
    if gRandom.Rndm() < acc:
        return True
    else:
        return False
    
