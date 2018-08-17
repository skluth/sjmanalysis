CDECK  ID>, PXLTH4.
      SUBROUTINE PXLTH4 (NTRAK,ITKDM,PTRAK,THRVAL,THRVEC,IERR)
*.*********************************************************
*. ------
*. PXLTH4
*. ------
*. Routine to determine the Thrust Principal, Major and
*. Minor axes and values using the Jetset algorithm.
*. The implementation here is without a common block, however.
*. Thus this routine may be used regardless of whether the
*. Jetset6.3 or Jetset7.1 library might be linked.  It is
*. not necessary to link to Jetset, however.
*. Usage     :
*.
*.      INTEGER  ITKDM,MXTRAK
*.      PARAMETER  (ITKDM=3.or.more,MXTRAK=1.or.more)
*.      INTEGER NTRAK,IERR
*.      REAL  PTRAK (ITKDM,MXTRAK),
*.     +      THRVEC (3,3.or.more),
*.     +      THRVAL (3.or.more)
*.
*.      NTRAK = 1.to.MXTRAK
*.      CALL  PXLTH4 (NTRAK,ITKDM,PTRAK,THRVAL,THRVEC,IERR)
*.
*. The thrust vectors THRVEC are ordered according to the
*. corresponding thrust values THRVAL such that
*.     THRVAL (1) < THRVAL (2) < THRVAL (3)
*. Thus THRVEC (*,3) is the Thrust Principal axis;
*. Thus THRVEC (*,2) is the Thrust Major axis;
*. Thus THRVEC (*,1) is the Thrust Minor axis;
*.
*. INPUT     : NTRAK    Total number of particles
*. INPUT     : ITKDM    First dimension of PTRAK array
*. INPUT     : PTRAK    Particle momentum array: Px,Py,Pz,E
*. OUTPUT    : THRVAL   Thrust values
*. OUTPUT    : THRVEC   Associated Thrust vectors
*. OUTPUT    : IERR     = 0 if all is OK ;   = -1 otherwise
*.
*. CALLS     : PXLUT3
*. CALLED    : By User
*.
*. AUTHOR    :  J.W.Gary
*. CREATED   :  18-Jun-88
*. LAST MOD  :  04-Feb-89
*.
*. Modification Log.
*. 04-Feb-89  Integrate with PXLUT3  J.W.Gary
*.
*.*********************************************************
      IMPLICIT NONE
      INTEGER  NRLUDM,IOLUN
      PARAMETER (NRLUDM=1000,IOLUN=6)
      INTEGER  NTRAK,AXIS,IX1,IX2,ITKDM,IERR
      REAL  PTRAK (ITKDM,*),PLUND (NRLUDM,5),
     +      THRVEC (3,*),THRVAL (*)
      REAL  THRUST,OBLATE
      LOGICAL  LPRT
      DATA  LPRT / .FALSE. /

      IERR = 0
      IF (NTRAK.LE.1.OR.NTRAK.GT.NRLUDM) THEN
          IERR = -1
          WRITE (IOLUN,FMT='('' PXLTH4: Error, NTRAK ='',I6/
     +           ''  Max. allowed ='',I6)') NTRAK,NRLUDM
          GO TO 990
      END IF
*  Pack 4-momenta in Jetset format
*  ---- --------- -- ------ ------
      DO 110  IX1 = 1,NTRAK
          DO 100  AXIS = 1,4
              PLUND (IX1,AXIS) = PTRAK (AXIS,IX1)
 100      CONTINUE
 110  CONTINUE
*  Jetset algorithm for Thrust
*  ------ --------- --- ------
      CALL PXLUT3 (NTRAK,NRLUDM,PLUND,THRUST,OBLATE)
      IF (LPRT) WRITE (IOLUN,FMT='('' PXLTH4:  THRUST,'',
     +     ''OBLATE ='',2E12.4)') THRUST,OBLATE
      IF (THRUST.LT.0) THEN
          IERR = -1
          GO TO 990
      END IF
*  Buffer eigenvectors for output
*  ------ ------------ --- ------
      DO 210  IX1 = 1,3
          IX2 = NTRAK + (4 - IX1)
          THRVAL (IX1) = PLUND (IX2,4)
          DO 200  AXIS = 1,3
              THRVEC (AXIS,IX1) = PLUND (IX2,AXIS)
 200      CONTINUE
 210  CONTINUE
      IF (LPRT) THEN
          WRITE (IOLUN,FMT='('' PXLTH4: THRVAL,THRVEC ='',
     +          3(/5X,4E12.4))') (THRVAL (IX1),
     +          (THRVEC (IX2,IX1),IX2=1,3),IX1=1,3)
      END IF

 990  RETURN
      END
CDECK  ID>, PXLUT3.
      SUBROUTINE PXLUT3 (N,NRLUDM,P,THR,OBL)
*.*********************************************************
*. ------
*. PXLUT3
*. ------
*. An "in-house" version of the Jetset thrust finding algorithm
*. which works entirely through an argument list rather than
*. with e.g. the Jetset common blocks.  This routine calculates
*. the standard linear thrust vectors and values. Its operation
*. is entirely decoupled from any MST or MSTJ variables etc.
*. which might be set by a user using Jetset.
*. The main purpose of developing an in-house version of the
*. Jetset thrust algorithm is so as to have a version
*. which is compatible with both Jetset6.3 and Jetset7.1 etc.
*. (because of the change in the Jetset common blocks between
*. these two versions, the Jetset version of this thrust
*. algorithm LUTHRU is version specific).
*.
*. The Jetset thrust algorithm implements an "approximate" method
*. for thrust determination because not all particle combinations
*. are tested.  It is therefore logically possible that the thrust
*. axes and values determined by this routine will correspond
*. to local rather than to absolute maxima of the thrust function.
*. However in practice this is unlikely because several starting
*. points are used and the algorithm iterated to cross check one
*. convergence vs. another.  Thus this routine offers a very fast
*. and in practice quite accurate algorithm for thrust (much faster
*. than so-called "exact" algorithms).
*. Usage     :
*.
*.      INTEGER  NRLUDM
*.      PARAMETER (NRLUDM=1000.or.so)
*.      REAL PLUND (NRLUDM,5)
*.      INTEGER  NTRAK
*.      REAL  THRUST,OBLATE
*.
*.      (define NTRAK, fill PLUND)
*.      CALL PXLUT3 (NTRAK,NRLUDM,PLUND,THRUST,OBLATE)
*.
*. INPUT     : NTRAK    Number of tracks
*. INPUT     : NRLUDM   First dimension of PLUND
*. INPUT     : P        4-momenta in Jetset format
*. OUTPUT    : THRUST   Thrust value
*. OUTPUT    : OBLATE   Oblateness value
*.
*. CALLS     : PXANXY,PXPLU3,PXRMX3,PXROF3,PXROB3
*. CALLED    : PXLTH4
*.
*. AUTHOR    : Modified from LUTHRU (T.Sjostrand) by J.W.Gary
*. CREATED   : 31-Jan-89
*. LAST MOD  : 22-Sep-00
*.
*. Modification Log.
*. 04-Feb-89  In-house version for PX library  J.W.Gary
*. 12-Mar-89  Get rid of calls to RLU          J.W.Gary
*. 27-Nov-95 M.Schroder  Clear part of the array P above tracks
*. 05-May-97 D.Chrisman  Remove declaration of unused variable
*.                       SG and function RLU.
*. 22-Sep-00 R.Yaari     Add double prec. calc. for Track ordering and P
*.                       calculations
*.
*.*********************************************************
      IMPLICIT NONE
      INTEGER  N,NP,MSTU44,MSTU45,ILC,ILD,ILF,ILG,I,J,
     +         IAGR,NC,IPP,IERR,NRLUDM
      DOUBLE PRECISION DPP(3),PARU48,DTDS,DTDI(3),TPR(3),THP,THPS,DTEMP
      REAL  P (NRLUDM,*),TDI (3),PVEC (3),
     +      RVEC (3),RMX (3,3)
      REAL  PS,PARU42,TDS,SGN,OBL,THR,
     +      PHI,CP,SP
      DATA  PARU42 / 1. /, PARU48 / 0.0001 /,
     +      MSTU44 / 4  /, MSTU45 / 2 /

      IF(2*N+MSTU44+15.GE.NRLUDM-5) THEN
          WRITE (6,FMT='('' PXLUT3: Error, not enough buffer'',
     +           ''space for Thrust calculation'')')
          THR=-2.
          OBL=-2.
          GO TO 990
      ENDIF
C  M.Schroder (these elements are always used, but sometimes not set...)
      DO 50 I = N+1, 2*N+2
          P(I,1) = 0.
          P(I,2) = 0.
          P(I,3) = 0.
          P(I,4) = 0.
          P(I,5) = 0.
   50 CONTINUE
C...Take copy of particles that are to be considered in thrust analysis.
      NP = 0
      PS = 0.
      DO 100 I = 1,N
          NP = NP + 1
          P (N+NP,1) = P (I,1)
          P (N+NP,2) = P (I,2)
          P (N+NP,3) = P (I,3)
          DPP(1)     = P (I,1)
          DPP(2)     = P (I,2)
          DPP(3)     = P (I,3)
C improve selection by using double precision calculation
          P (N+NP,4) = SQRT (DPP(1)**2 +DPP(2)**2 + DPP(3)**2)
          P (N+NP,5) = 1.
          IF (ABS (PARU42-1.).GT.0.001)
     +        P (N+NP,5) = P (N+NP,4)**(PARU42-1.)
          PS = PS + P (N+NP,4) * P (N+NP,5)
  100 CONTINUE
C...Loop over thrust and major. T axis along z direction in latter case.
      DO 280 ILD=1,2
          IF (ILD.EQ.2) THEN
              CALL PXANXY (P (N+NP+1,1),P (N+NP+1,2),PHI,IERR)
              CALL PXPLU3 (N+NP+1,NRLUDM,P,PVEC,'U')
              CP = COS (PHI)
              SP = SIN (PHI)
              CALL PXRMX3 (PVEC,CP,SP,RMX)
              DO 105 IPP = N+1,N+NP+1
                  CALL PXPLU3 (IPP,NRLUDM,P,PVEC,'U')
                  CALL PXROF3 (RMX,PVEC,RVEC)
                  CALL PXPLU3 (IPP,NRLUDM,P,RVEC,'P')
  105         CONTINUE
          ENDIF
C...Find and order particles with highest p (pT for major).
          DO 110 ILF = N+NP+4,N+NP+MSTU44+4
              P (ILF,4) = 0.
  110     CONTINUE
          DO 150 I = N+1,N+NP
C improve ordering by using double precision calculation
              IF (ILD.EQ.2) THEN
                DPP(1) = P(I,1)
                DPP(2) = P(I,2)
                P(I,4) = SQRT (DPP(1)**2 + DPP(2)**2)
              ENDIF
              DO 120 ILF = N+NP+MSTU44+3,N+NP+4,-1
                  IF (P (I,4).LE.P (ILF,4)) GO TO 130
                  DO 115 J=1,5
                      P(ILF+1,J)=P(ILF,J)
  115             CONTINUE
  120         CONTINUE
              ILF = N + NP + 3
  130         DO 140 J=1,5
                  P (ILF+1,J) = P (I,J)
  140         CONTINUE
  150     CONTINUE
C...Find and order initial axes with highest thrust (major).
          DO 160 ILG=N+NP+MSTU44+5,N+NP+MSTU44+15
              P(ILG,4)=0.
  160     CONTINUE
          NC = 2**(MIN (MSTU44,NP) - 1)
          DO 220 ILC=1,NC
              DO 170 J=1,3
                   TDI(J)=0.
                  DTDI(J)=0.D0
  170         CONTINUE
              DO 180 ILF=1,MIN(MSTU44,NP)
                  SGN = P (N+NP+ILF+3,5)
                  IF (2**ILF*((ILC+2**(ILF-1)-1)/2**ILF).GE.ILC)
     +                SGN = -SGN
                  DO 175 J = 1,4-ILD
C SGN is really only 1 or -1. How I love these multiplexed arrays :(
                      DTDI (J) = DTDI (J) + SGN * P (N+NP+ILF+3,J)
                       TDI (J) = DTDI (J)
  175             CONTINUE
  180         CONTINUE
              DTDS = DTDI (1)**2 + DTDI (2)**2 + DTDI (3)**2
               TDS = DTDS
              DO 190 ILG = N+NP+MSTU44+MIN(ILC,10)+4,N+NP+MSTU44+5,-1
                  IF (TDS.LE.P (ILG,4)) GO TO 200
                  DO 185 J=1,4
                      P (ILG+1,J) = P (ILG,J)
  185             CONTINUE
  190         CONTINUE
              ILG=N + NP + MSTU44 + 4
  200         DO 210 J=1,3
                  P (ILG+1,J) = TDI (J)
  210         CONTINUE
              P (ILG+1,4) = TDS
  220     CONTINUE
C...Iterate direction of axis until stable maximum.
          P (N+NP+ILD,4) = 0.
          ILG = 0
  230     ILG = ILG + 1
          THP = 0.
  240     THPS = THP
          DO 250 J=1,3
              IF (THP.LE.1E-10) TDI (J) = P (N+NP+MSTU44+4+ILG,J)
              IF (THP.GT.1E-10) TDI (J) = TPR (J)
              TPR (J) = 0.
  250     CONTINUE
          DO 260 I = N+1,N+NP
              SGN = SIGN (P(I,5),
     +                TDI(1)*P(I,1)+TDI(2)*P(I,2)+TDI(3)*P(I,3))
              DO 255 J = 1,4-ILD
                  TPR (J) = TPR (J) + SGN * P (I,J)
  255         CONTINUE
  260     CONTINUE
          THP = SQRT (TPR (1)**2 + TPR (2)**2 + TPR (3)**2) / PS
          IF (THP.GE.THPS+PARU48) GO TO 240
C...Save good axis. Try new initial axis until a number of tries agree.
          IF (THP.LT.P(N+NP+ILD,4)-PARU48.AND.ILG.LT.MIN(10,NC))
     +          GO TO 230
          IF (THP.GT.P(N+NP+ILD,4)+PARU48) THEN
              IAGR = 0
**JWG              SGN = (-1.)**INT (RLU(0)+0.5)
              SGN = 1.
              DO 270 J=1,3
                  P (N+NP+ILD,J) = SGN * TPR (J) / (PS*THP)
  270         CONTINUE
              P(N+NP+ILD,4) = THP
              P(N+NP+ILD,5) = 0.
          ENDIF
          IAGR = IAGR + 1
          IF (IAGR.LT.MSTU45.AND.ILG.LT.MIN(10,NC)) GO TO 230
  280 CONTINUE
C...Find minor axis and value by orthogonality.
**JWG      SGN = (-1.)**INT (RLU(0)+0.5)
      SGN = 1.
      P (N+NP+3,1) = -SGN * P (N+NP+2,2)
      P (N+NP+3,2) = SGN * P (N+NP+2,1)
      P (N+NP+3,3) = 0.
      THP = 0.
      DO 290 I = N+1,N+NP
          THP = THP + P (I,5)
     +        * ABS (P (N+NP+3,1) * P (I,1) + P (N+NP+3,2) * P (I,2))
  290 CONTINUE
      P (N+NP+3,4) = THP / PS
      P (N+NP+3,5) = 0.
C...Fill axis information. Rotate back to original coordinate system.
      DO 300 ILD = 1,3
          DO 295 J = 1,5
              P (N+ILD,J) = P (N+NP+ILD,J)
  295     CONTINUE
  300 CONTINUE
      DO 305 IPP = N+1,N+3
          CALL PXPLU3 (IPP,NRLUDM,P,PVEC,'U')
          CALL PXROB3 (RMX,PVEC,RVEC)
          CALL PXPLU3 (IPP,NRLUDM,P,RVEC,'P')
  305 CONTINUE
C...Select storing option. Calculate thurst and oblateness.
      THR = P (N+1,4)
      OBL = P (N+2,4) - P (N+3,4)

  990 RETURN
      END
CDECK  ID>, PXPLU3.
      SUBROUTINE PXPLU3 (INDX,NRLUDM,PLUND,PVEC,CHAR)
*.*********************************************************
*. ------
*. PXPLU3
*. ------
*. A utility routine to repack a Jetset 3-momentum as a
*. standard ordered array or vice-versa.  This routine
*. is used to translate between the array conventions
*. of the Jetset thrust algorithm and of the other routines
*. in this library
*. Usage     :
*.
*.      CALL PXPLU3 (INDX,NRLUDM,PLUND,PVEC,CHAR)
*.
*. INPUT     : INDX     The Jetset vector number
*. IN/OUT    : NRLUDM   First argument of PLUND
*. IN/OUT    : PLUND    The Jetset 3-momentum array
*. IN/OUT    : PVEC     The input or output array
*. CONTROL   : CHAR     ='U' means unpack Jetset array
*.                      = anything else means pack Jetset array
*.
*. CALLS     : none
*. CALLED    : PXLUT3
*.
*. AUTHOR    :  J.W.Gary
*. CREATED   :  04-Feb-89
*. LAST MOD  :  04-Feb-89
*.
*. Modification Log.
*.
*.*********************************************************
      IMPLICIT NONE
      INTEGER  INDX,NRLUDM,IX
      REAL  PVEC (*),PLUND (NRLUDM,*)
      CHARACTER*1  CHAR
      DO 120  IX = 1,3
          IF (CHAR.EQ.'U') THEN
              PVEC (IX) = PLUND (INDX,IX)
          ELSE
              PLUND (INDX,IX) = PVEC (IX)
          END IF
 120  CONTINUE
      RETURN
      END
CDECK  ID>, PXPLU5.
      SUBROUTINE  PXPLU5 (IFST,ILST,NLUPDM,PLUND,NTRKDM,PTRK,COPT)
*.*********************************************************
*. ------
*. PXPLU5
*. ------
*. Routine to pack or unpack the Jetset "P array," which is packed as
*. P (itrak,ix),ix=1,5 = P [(itrak,(Px,Py,Pz,E,M)], into a standard
*. ordered 5-vector array PTRK, i.e. PTRK [(Px,Py,Pz,E,M),itrak].
*. This routine works with an argument list rather than directly
*. with the Jetset common block; thus it can be used with any version
*. of Jetset (JETSET63, JETSET71, etc.).
*.
*. INPUT     : IFST    First Jetset vector to copy
*. INPUT     : ILST    Last Jetset vector to copy
*. INPUT     : NLUPDM  1st dimension of Jetset "P" array
*. IN/OUTPUT : PLUND   The Jetset "P" array from /LUJETS/
*. INPUT     : NTRKDM  1st dimension of PTRK (= 5 or more)
*. IN/OUTPUT : PTRK    The standard ordered array
*. INPUT     : COPT    A switch between packing or unpacking
*.                     COPT = 'P' for packing PTRK into PLUND
*.                      = anything else for packing PLUND into PTRK
*.
*. CALLS:  none
*. CALLED: Various
*.
*. AUTHOR    :  J.W.Gary
*. CREATED   :  04-Feb-88
*. LAST MOD  :  04-Feb-88
*.
*. Modification Log.
*.
*.*********************************************************
      IMPLICIT NONE
      INTEGER  NLUPDM,NTRKDM,IFST,ILST,IP,IX,NVEC
      REAL  PLUND (NLUPDM,*),PTRK (NTRKDM,*)
      CHARACTER*1  COPT

      NVEC = ILST - IFST + 1
      IF (COPT.EQ.'P') THEN
          DO 120  IP = 1,NVEC
              DO 110  IX = 1,5
                  PLUND (IFST+IP-1,IX) = PTRK (IX,IP)
 110          CONTINUE
 120      CONTINUE
      ELSE
          DO 140  IP = 1,NVEC
              DO 130  IX = 1,5
                  PTRK (IX,IP) = PLUND (IFST+IP-1,IX)
 130          CONTINUE
 140      CONTINUE
      END IF

      RETURN
      END
CDECK  ID>, PXANXY.
      SUBROUTINE PXANXY (XX,YY,ANG,IERR)
*.*********************************************************
*. ------
*. PXANXY
*. ------
*. SOURCE: Jetset7.1 (T. Sjostrand)
*. Reconstruct the azimuthal angle of a vector,
*. given the X and Y components of a vector
*. Usage     :
*.
*.      INTEGER  IERR
*.      REAL  XX,YY,ANG
*.
*.      CALL PXANXY (XX,YY,ANG,IERR)
*.
*. INPUT     : XX      The X component of a vector
*. INPUT     : YY      The Y component of a vector
*. OUTPUT    : ANG     The azimuthal angle
*. OUTPUT    : IERR    = 0 if all is OK ;   = -1 otherwise
*.
*.*********************************************************
      IMPLICIT NONE
      REAL  PIII
      PARAMETER  (PIII=3.1415927)
      INTEGER  IERR
      REAL  XX,YY,ANG
      DOUBLE PRECISION  ULANGL,RRR,XXX,YYY

      IERR = 0
      XXX = XX
      YYY = YY
      RRR = DSQRT (XXX**2 + YYY**2)
      IF (RRR.LT.1E-20) GO TO 990
      IF ((DABS (XXX)/RRR).LT.0.8) THEN
          ULANGL = DSIGN (DACOS (XXX/RRR),YYY)
      ELSE
          ULANGL = DASIN (YYY/RRR)
          IF (XXX.LT.0..AND.ULANGL.GE.0.) THEN
              ULANGL = PIII - ULANGL
          ELSE IF (XXX.LT.0.) THEN
              ULANGL = - PIII - ULANGL
          END IF
      END IF
      ANG = ULANGL

      RETURN
 990  IERR = -1
      RETURN
      END
CDECK  ID>, PXRMX3.
      SUBROUTINE PXRMX3 (SVECT,CP,SP,RMX)
*.*********************************************************
*. ------
*. PXRMX3
*. ------
*. SOURCE: HERWIG (B.Webber,G.Marchesini)
*. Calculate the rotation matrix to get from vector VECT
*. to the Z axis, followed by a rotation PHI about the Z axis,
*. where CP, SP are the cosine and sine of PHI, respectively.
*. Usage     :
*.
*.      REAL  SVECT (3.or.more),
*.     +      RMX (3,3.or.more)
*.      REAL  CP,SP
*.
*.      CALL PXRMX3 (SVECT,CP,SP,RMX)
*.
*. INPUT     : VECT   The vector for which the rotation matrix
*.                    is to be calculated
*. INPUT     : CP     Cosine of phi
*. INPUT     : SP     Sine of phi
*. OUTPUT    : RMX    The rotation matrix
*.
*. LAST MOD :  9-Dec-00
*.
*. Modification Log:
*.  9-Dec-00 R.Yaari   Make PTCUT also double precision
*. 22-Sep-00 R.Yaari   Use double precision copy of array SVECT for calculations
*.*********************************************************
      IMPLICIT NONE
      REAL  SVECT (*),RMX (3,*)
      DOUBLE PRECISION  CT,ST,CF,SF,PP,PT,VECT(3),PTCUT
      REAL  CP,SP
      DATA  PTCUT / 1.D-10 /
      VECT(1) = SVECT(1)
      VECT(2) = SVECT(2)
      VECT(3) = SVECT(3)
      PT = VECT (1)**2 + VECT (2)**2
      IF (PT.LT.PTCUT) THEN
         CT = SIGN (1.0D0,VECT (3))
         ST = 0.
         CF = 1.
         SF = 0.
      ELSE
         PP = SQRT (VECT (3)**2 + PT)
         PT = SQRT (PT)
         CT = VECT (3) / PP
         ST = PT / PP
         CF = VECT (1) / PT
         SF = VECT (2) / PT
      END IF
      RMX (1,1) =  (CP * CF * CT) + (SP * SF)
      RMX (1,2) =  (CP * SF * CT) - (SP * CF)
      RMX (1,3) = -(CP * ST)
      RMX (2,1) = -(CP * SF) + (SP * CF * CT)
      RMX (2,2) =  (CP * CF) + (SP * SF * CT)
      RMX (2,3) = -(SP * ST)
      RMX (3,1) =  (CF * ST)
      RMX (3,2) =  (SF * ST)
      RMX (3,3) =  CT
      RETURN
      END
CDECK  ID>, PXROB3.
      SUBROUTINE PXROB3 (SRMX,SVECT,RVEC)
*.*********************************************************
*. ------
*. PXROB3
*. ------
*. SOURCE: HERWIG (B.Webber,G.Marchesini)
*. Rotate 3-vector VECT by inverse of rotation matrix RMX,
*.      RVEC = (RMX)-1 * VECT
*. Usage     :
*.
*.      REAL  VECT (3.or.more),
*.     +      RVEC (3.or.more),
*.     +      RMX  (3,3.or.more)
*.
*.      CALL PXROB3 (RMX,VECT,RVEC)
*.
*. INPUT     : RMX    The rotation matrix
*. INPUT     : VECT   The vector to be rotated
*. OUTPUT    : RVEC   The rotated vector
*.
*. LAST MOD : 22-Sep-00
*.
*. Modification Log:
*.
*. 22-Sep-00 R.Yaari   Use double precision copy of arrays for calculations
*.*********************************************************
      IMPLICIT NONE
      DOUBLE PRECISION  S1,S2,S3, RMX(3,3), VECT(3)
      REAL SRMX (3,*),SVECT (*),RVEC (*)
      INTEGER I,J
      DO 100 I = 1,3
        VECT(I) = SVECT(I)
        DO 150 J=1,3
          RMX(I,J) = SRMX(I,J)
  150   CONTINUE
  100 CONTINUE
      S1 = VECT (1) * RMX (1,1) + VECT (2) * RMX (2,1)
     +   + VECT (3) * RMX (3,1)
      S2 = VECT (1) * RMX (1,2) + VECT (2) * RMX (2,2)
     +   + VECT (3) * RMX (3,2)
      S3 = VECT (1) * RMX (1,3) + VECT (2) * RMX (2,3)
     +   + VECT (3) * RMX (3,3)
      RVEC (1) = S1
      RVEC (2) = S2
      RVEC (3) = S3
      RETURN
      END
CDECK  ID>, PXROF3.
      SUBROUTINE PXROF3 (SRMX,SVECT,RVEC)
*.*********************************************************
*. ------
*. PXROF3
*. ------
*. SOURCE: HERWIG (B.Webber,G.Marchesini)
*. Rotate 3-vector VECT by rotation matrix RMX,
*.      RVEC = RMX * VECT
*. Usage     :
*.
*.      REAL  VECT (3.or.more),
*.     +      RVEC (3.or.more),
*.     +      RMX (3,3.or.more)
*.
*.      CALL PXROF3 (RMX,VECT,RVEC)
*.
*. INPUT     : RMX    The rotation matrix
*. INPUT     : VECT   The vector to be rotated
*. OUTPUT    : RVEC   The rotated vector
*.
*. LAST MOD : 22-Sep-00
*.
*. Modification Log:
*.
*. 22-Sep-00 R.Yaari   Use double precision copy of arrays for calculations
*.*********************************************************
      IMPLICIT NONE
      DOUBLE PRECISION  S1,S2,S3, RMX(3,3), VECT(3)
      REAL SRMX (3,*),SVECT (*),RVEC (*)
      INTEGER I,J

      DO 100 I = 1,3
        VECT(I) = SVECT(I)
        DO 150 J=1,3
          RMX(I,J) = SRMX(I,J)
  150   CONTINUE
  100 CONTINUE
 
      S1 = RMX (1,1) * VECT (1) + RMX (1,2) * VECT (2)
     +   + RMX (1,3) * VECT (3)
      S2 = RMX (2,1) * VECT (1) + RMX (2,2) * VECT (2)
     +   + RMX (2,3) * VECT (3)
      S3 = RMX (3,1) * VECT (1) + RMX (3,2) * VECT (2)
     +   + RMX (3,3) * VECT (3)
      RVEC (1) = S1
      RVEC (2) = S2
      RVEC (3) = S3
      RETURN
      END
