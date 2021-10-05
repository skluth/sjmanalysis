C===================================================================
C
C  CODE FOR YCLUS JET FINDING ALGORITHM, FEATURING THE E0-, P-, P0-
C AND THE E-VARIANTS OF THE ORIGINAL JADE JET FINDER, AS WELL AS THE
C        NEW DURHAM (D) AND GENEVA (G) JET FINDING SCHEMES.
C
C REFERENCE:  SEE "JETS IN Z0 DECAYS", S. BETHKE, OPAL-CR084 AND
C             HD-PY 92/12 (AACHEN QCD WORKSHOP 1992), AND REFS QUOTED
C             THEREIN.
C
C THE YCLUS JET FINDING PACKAGE CONSISTS OF 7 SUBROUTINES.
C INTERFACE ROUTINES TO YCLUS, NAMELY PXJRC4, PXJRZ4, PXJTR4 AND PXYCR4,
C ARE PROVIDED WITH THIS PACKAGE, WHICH FUNCTION EXACTLY LIKE THE OLD
C PXLIB ROUTINES OF THE SAME NAME.
C YCLUS PROVIDES THE JADE JET FINDER VARIANTS E0, P, P0 AND E, AS WELL
C AS THE NEW DURHAM AND THE GENEVA SCHEMES. THESE ROUTINES ARE ABOUT
C FACTORS OF 2 TO 3 FASTER THAN THE OLD PXLIB ROUTINES, FOR AVERAGE
C PARTICLE MULTIPLICITIES AS FOUND IN Z0 EVENTS AT LEP.
C FOR HIGHER MULTIPLICITIES, THE GAIN IN SPEED IS EVEN LARGER!.
C IN ORDER TO MAKE FULL USE OF ALL THE NEW FEATURES OF THE YCLUS PACKAGE
C IT IS RECOMMENDED TO SHIFT TO THE LATTER ONE; THE YCLUS ROUTINES
C ARE EASY AND STRAIGHT FORWARD TO USE; THEY CONTAIN FEWER CALLING
C PARAMETERS AND PROVIDE MORE FEATURES THAN THE PXLIB ROUTINES.
C
C THE SUBROUTINES OF THE YCLUS PACKAGE (IT'S SELF-CONTAINED) ARE:
C
C  YKERN(IMODE,NT,ITKDM,PP,IERR)  RECONSTRUCTS AN EVENT DOWN TO 1-JET.
C                           MUST BE CALLED FIRST FOR EACH EVENT. USEFUL
C                           INFORMATION FOR FURTHER USE AND CALLS IS
C                           STORED IN COMMON /YCL/.
C                           IMODE: JET FINDER MODE (1...4 = JADE-E0, P,
C                                  P0 AND E, RESPECTIVELY; 5=D, 6=G).
C                           NT:    NUMBER OF INPUT PARTICLE 4-VECTORS.
C                           ITKDM: LENGTH OF VECTORS IN ARRAY PP.
C                           PP(ITKDM,*): ARRAY OF PARTICLE 4-MOMENTA.
C                           IERR:  RETURNED ERROR CODE (=0 IF O.K.)
C
C ONCE YKERN WAS CALLED AND IERR=0, THE FOLLOWING ROUTINES MAY BE USED
C IN ORDER TO EXTRACT FURTHER EVENT AND JET INFORMATION; THESE ROUTINES
C UTILIZE THE DATA WHICH WERE GENERATED AND STORED WITHIN YKERN.
C
C  YNJET(YCUT,NJET,IERR):   RETURNS THE NUMBER OF JETS FOR A GIVEN YCUT
C
C  YYJET(NJET,YL,YH,IERR):  RETURNS THE TWO VALUES OF YCUT BETWEEN WHICH
C                           THE EVENT IS CLASSIFIED AS NJET EVENT; YL<YH
C
C  YAXES(NJET,PNJ,IERR):    RETURNS THE JET AXES IF EVENT IS CLASSIFIED
C                           AS NJET EVENT (NJET=1 ... 10).
C                           PNJ(IL,IV): ARRAY OF 4-VECTORS OF AXES.
C
C  YASSO(NJET,PNJ,BL,IERR): RETURNS THE JET AXES FOR EVENT IF CLASSIFIED
C                           AS NJET EVENT (NJET=1 ... 10), AS WELL
C                           AS A POINTER WHICH, FOR EACH PARTICLE IN PP
C                           (SEE YKERN), POINTS TO THE JET TO WHICH IT
C                           IS ASSOCIATED WITH. JET AXES ARE ORDERED IN
C                           ENERGY:FIRST VECTOR IN PNJ IS MOST ENERGETIC
C                           BL(N): FOR N-TH PARTICLE IN ARRAY PP,
C                                  CONTAINS THE POSITION OF JET IN PNJ
C                                  WHICH PARTICLE N IS ASSOCIATED WITH.
C
C  YREAS(NJET,PNJ,YMIN,BL,IERR): RETURNS JET AXES FOR EVENT CLAS-
C                           SIFIED AS N-JET (SEE YASSO); HOWEVER AFTER
C                           REASSIGNING PARTICLES TO THEIR CLOSEST JET
C                           (BY ANGLE), STARTING WITH THE INITIAL AXES
C                           FROM LAST CALL TO YKERN, AND ITERATING NEW
C                           AXES UNTIL STABILITY IS REACHED. INPUTS (NT,
C                           PP) AS FOR YKERN AND NJET AS FOR YASSO; OUT-
C                           PUTS PNJ AND BL AS FOR YASSO. YMIN GIVES THE
C                           NEW VALUE OF YCUT WHERE EVENT, AFTER REASSI-
C                           GNMENT, FLIPS TO (NJET-1)-CONFIGURATION
C                           (COMPARE TO YREC(NJET-1) OF ORIGINAL CONF.).
C                           YREAS DOES NOT CHANGE COMMON /YCL/.
C
C
C  YTREE(LPRINT,PTR,IERR):  GENERATES COMPLETE EVENT HISTORY AND TREE
C                           STRUCTURE FROM JET RECONSTRUCTION DOWN TO
C                           1-JET, SIMILAR TO TREE INFORMATION AS KNOWN
C                           FROM MONTE CARLO GENERATORS. IF NT IS NUMBER
C                           OF PARTICLES WITH WHICH YKERN WAS CALLED BE-
C                           FORE, THE TREE IS STORED IN ARRAY PTR(10,K),
C                           WHERE K=1 TO NT-1 ARE ALL VECTORS FROM WITH-
C                           IN THE RECOMBINATION STAGE, AND VECTORS K=NT
C                           TO 2*NT-1 ARE THE ORIGINAL, FINAL PARTICLES;
C                           EACH VECTOR CARRIES ADDITIONAL INFORMATION
C                           ABOUT ITS PARENT AND DAUGHTERS, ITS MASS AND
C                           THE VALUE OF Y OF THE NEXT SPLITTING.
C
C COMMON /YCL/YREC(10),PJET(10,10,10),HISTOR(2,NYCLMX)
C        THIS COMMON IS FILLED BY YKERN AND CONTAINS ALL THE INFORMATION
C        WHICH IS REQUIRED FOR FURTHER ANALYSIS IN TYPICAL JET STUDIES.
C        YREC(N)     : VALUE OF YCUT FOR WHICH EVENT FLIPS FROM N+1 TO
C                      N-JET CONFIGURATION.
C        PJET(I,J,K) : 4-VECTORS (I=1..4) OF JET AXIS J WHEN EVENT IS
C                      CLASSIFIED TO HAVE K JETS (I=7 HOLDS # OF PARTI-
C                      CLES ASSOCIATED WITH THIS JET)
C        HISTOR(2,NYCLMX):CONTAINS CODED HISTORY OF RECOMBINATIONS DONE
C                      DOWN TO 1-JET CONFIGURATION; BETTER LEAVE IT
C                      TO YASSO, YREAS AND YTREE TO UNPACK THIS AND TO
C                      PROVIDE ACTUAL PARTICLE-JET ASSOCIATIONS.
C                      (NYCLMX = 250 AT THIS TIME; ENLARGE IF NEEDED).
C
C USE THE INFORMATION IN /YCL/ DIRECTLY OR CALL THE OTHER YXXXX ROUTINES
C FOR CONVENIENT UNPACKING OF INFORMATION.
C
C FOR FURTHER DETAILS, AND FOR THE PX-LIB INTERFACE ROUTINES, SEE THE
C HEADERS OF THE SUBROUTINES THEMSELVES.
C
C SUGGESTIONS/COMPLAINTS/PROBLEMS/DONATIONS TO SIGGI@CERNVM, PLEASE.
C
C HISTORY:
C ========
C 10 NOV 92  S. BETHKE  - 1ST SET-UP OF YCLUS PACKAGE
C 11 NOV 92  S. BETHKE  - REMOVE BUG IN CODE FOR JADE E0 SCHEME
C 12 NOV 92  S. BETHKE  - REMOVE INCONSISTENCIES FOR LOW MULTIPL. EVENTS
C                         (NT<=10) IN YKERN, YASSO
C                       - ADD UTILITY ROUTINE YREAS
C 25 NOV 92  S. BETHKE  - VARIABLE INPUT VECTOR LENGTH IN YKERN.
C                         (!! NOTE !! ONE MORE INPUT PARAMETER TO YKERN)
C                       - REDUCE SIZE OF AND RESTRUCTURE ARRAY "HISTOR"
C                       - ADD UTILITY ROUTINE YTREE.
C                       - STORE INPUT VECTOR ARRAY IN COMMON /YINT/ WHEN
C                         CALLING YKERN; ALL OTHER ROUTINES ARE THEN
C                         CALLED WITHOUT PASSING INPUT VECTORS AGAIN
C                         (NOTE: CALLING PARAMETERS FOR YREAS REDUCED).
C 09 MAR 93  S. BETHKE  - CODE CHANGED TO STANDARD FORTRAN 77
C 05 FEB 96 D. Chrisman - Add a WRITE statement in YYJET to print
C                         a message if (NJET.GT.NTO)
C 05 FEB 96 D. Chrisman - Initialize NJET=1 just before the
C                         "DO 5001 ..."-loop in subroutine YNJET.
C 12 May 97 D. Chrisman - remove declaration of unused variable NPMAX.
C 08 Oct 97 D. Chrisman - Increase NYCLMX from 250 to 500, everywhere.
C 09 Apr 98 D. Chrisman - Add PXCAMJ: CAMBRIDGE JET CLUSTERING ALGORITHM
C 21 JUN 98 D. Chrisman - For LEP2 studies we need to handle more
C                         than 10 jets. The max number of jets allowed
C                         was hard coded in
C                         many of the YKERN subroutines. The parameter
C                         NJETMX was introduced in YKERN, YNJET, YYJET,
C                         YAXES, YASSO, YREAS and YTREE in order to remove
C                         the harded coded limit of 10 jets.
C======================================================================
CDECK  ID>, YKERN.
      SUBROUTINE YKERN(IMODE,NT,ITKDM,PP,IERR)
C
C  JET FINDING ROUTINE A LA THE ORIGINAL JADE (E0) SCHEME, INCLUDING
C  THE E-,P- AND P0-VARIANTS. ALSO FEATURES THE NEW DURHAM (D)
C  AND GENEVA (G) SCHEMES.
C
C  INPUTS: IMODE (JET FINDING SCHEME); PP-ARRAY CONTAINING THE
C  FOUR-MOMENTA OF THE SELECTED PARTICLES, NT SPECIFYING HOW
C  MANY LOCATIONS IN PP ARE FILLED.
C
C  OUTPUTS:  YREC(I) CONTAINS VALUE OF Y WHERE EVENT FLIPS FROM
C  (I+1)-JET TO I-JET; PJET(K,I,J) CONTAINS THE JET AXES FOUR VECTOR
C  (K=1-4) OF JET NUMBER I WHEN THE EVENT IS CLASSIFIED AS J-JET.
C  LOCATION K=7 GIVES THE NUMBER OF PARTICLES ASSIGNED TO THIS JET.
C  NOTE THAT EACH EVENT IS ALWAYS FULLY RECONSTRUCTED DOWN TO
C  1-JET CONFIGURATION.
C  PROGRAM IS SELF-CONTAINED.
C
C  IMODE=1: E0 (=JADE) 2: P  3: P0  4: E  5: DURHAM (KT)  6: GENEVA
C
C  LAST MOD  :  21-Jun-98
C
C  Modification Log.
C  08-Oct-97 D. Chrisman, Increase NYCLMX from 250 to 500.
C  21-Jun-98 D. Chrisman, Want to be able to handle more than 10 jets.
C                         Introduce parameter NJETMX. Use this
C                         in declaration of PJET and YREC in order to remove
C                         hard coded numbers. Remove variable JCHECK and
C                         replace with NJETMX.
C
      IMPLICIT NONE
C     IMPLICIT NONE
      INTEGER NT,ITKDM,NYCLMX,NJETMX,NCALL,NPRINT,I,J,K,IERR,IMODEO
      INTEGER IMODE
      PARAMETER (NYCLMX = 500, NJETMX = 20)
      INTEGER NTO,NJETO,NJJ,IM,JM,NOLD,KK,II,HISTOR(2,NYCLMX)
      REAL PL(7,NYCLMX),Y(NYCLMX,NYCLMX),PINT(10,NYCLMX)
      REAL YREC(NJETMX),EVISO
      REAL PP(ITKDM,*),JADE,D,G,E,EVIS,PVIS,PJET(10,NJETMX,NJETMX),YMINI
      COMMON /YCL/YREC,PJET,HISTOR
      COMMON /YINT/ IMODEO,NTO,NJETO,PINT
      CHARACTER*7 CM
      SAVE NCALL,NPRINT,Y
      DATA NCALL,NPRINT /0,0/

Cf2py intent(in) IMODE
Cf2py integer intent(in), depend(pp) :: nt=shape(pp,1), itkdm=shape(pp,0)
Cf2py intent(in) pp
Cf2py intent(out) IERR
      
C
C  JET RESOLUTION FUNCTIONS
C
      JADE(I,J) = 2.*PL(4,I)*PL(4,J)*MAX(0.,(1.-
     + (PL(1,I)*PL(1,J)+PL(2,I)*PL(2,J)+PL(3,I)*PL(3,J))/
     + (PL(6,I)*PL(6,J))))
      D(I,J) = 2.*MIN(PL(4,I)*PL(4,I),PL(4,J)*PL(4,J))*MAX(0.,(1.-
     + (PL(1,I)*PL(1,J)+PL(2,I)*PL(2,J)+PL(3,I)*PL(3,J))/
     + (PL(6,I)*PL(6,J))))
      E(I,J) = MAX(0.,(PL(4,I)+PL(4,J))**2-(PL(1,I)+PL(1,J))**2-
     + (PL(2,I)+PL(2,J))**2-(PL(3,I)+PL(3,J))**2)
      G(I,J) = 8.*PL(4,I)*PL(4,J)*MAX(0.,(1.-
     + (PL(1,I)*PL(1,J)+PL(2,I)*PL(2,J)+PL(3,I)*PL(3,J))/
     + (PL(6,I)*PL(6,J))))/(9.*(PL(4,I)+PL(4,J))**2)
C
C INITIALIZE
C
      IERR = -1
      IF(NCALL.LE.0) THEN
        IMODEO = 0
        NTO = 0
        NJETO = 0
      ENDIF
      DO 5001 I=1,NJETMX
        YREC(I) = 0.
 5001 CONTINUE
      NCALL = NCALL + 1
C
C PRINT JET SCHEME
C
      IF(IMODE.NE.IMODEO .AND. NPRINT.LE.7) THEN
        IF(IMODE.EQ.1) THEN
          CM = 'JADE E0'
        ELSEIF(IMODE.EQ.2) THEN
          CM = 'JADE P '
        ELSEIF(IMODE.EQ.3) THEN
          CM = 'JADE P0'
        ELSEIF(IMODE.EQ.4) THEN
          CM = 'JADE E '
        ELSEIF(IMODE.EQ.5) THEN
          CM = 'DURHAM '
        ELSEIF(IMODE.EQ.6) THEN
          CM = 'GENEVA '
        ELSE
          WRITE(6,281) IMODE
 281      FORMAT(/,' ### YKERN: IMODE =',I3,' INVALID; SET TO 1 ###')
          IMODE = 1
          CM = 'JADE E0'
        ENDIF
        PRINT 789, CM
 789    FORMAT(/,8X,54('#'),/,8X,
     +  '# YCLUS JET FINDER WITH ',A7,' RECOMBINATION SCHEME #',/,
     +  8X,54('#'),/)
        NPRINT = NPRINT + 1
      ENDIF
C
C CHECK INPUT PARAMETERS
C
      IF(NT.LE.1) THEN
        PRINT 701, NT
 701    FORMAT(/,' ### YKERN: NUMBER OF INPUT PARTICLES ',I4,' < 2;',
     +  ' NO JET RECONSTRUCTION DONE. ###')
        RETURN
      ENDIF
      IF(NT.GT.NYCLMX) THEN
        PRINT 700, NT, NYCLMX,NYCLMX
  700   FORMAT(/,' #### YKERN: NUMBER OF INPUT PARTICLES ',I4,' > ',I4,
     +  '; SET TO ',I4,/,12X,' INCREASE NYCLMX IF THIS',
     +  ' WARNING OCCURS MORE OFTEN. #### ')
        NT = NYCLMX
      ENDIF
C
C COPY INPUT VECTORS INTO INTERNAL MOMENTUM ARRAY (PL)
C  (POSITION  7 OF EACH VECTOR IN PL WILL BE OVERWRITTEN BY THIS ROUTINE
C   IN ORDER TO RECORD NUMBER OF INITIAL PARTICLES BELONGING TO THIS
C   POSITION [= 1. INITIALLY])
C
      EVIS = 0.
      PVIS = 0.
      DO 5002 I=1,NT
        PL(6,I)=SQRT(PP(1,I)**2 + PP(2,I)**2 + PP(3,I)**2)
        PL(1,I)=PP(1,I)
        PL(2,I)=PP(2,I)
        PL(3,I)=PP(3,I)
        IF(IMODE.EQ.2 .OR. IMODE.EQ.3) THEN
          PL(4,I) = PL(6,I)
        ELSE
          PL(4,I) = PP(4,I)
        ENDIF
        PL(7,I) = 1.
        DO 5003 K=1,7
          PINT(K,I) = PL(K,I)
 5003   CONTINUE
        EVIS=EVIS+PL(4,I)
        PVIS=PVIS+PL(6,I)
 5002 CONTINUE
      NJJ = NT
C
      IF(EVIS.LE.0. .OR. PVIS.GT.1.001*EVIS) THEN
        WRITE(6,283) EVIS,PVIS
 283    FORMAT(' #### YKERN: INCOMPATIBLE SUMS OF ENERGIES',
     +  ' AND/OR MOMENTA:',
     +  /,12X,' SUM(E) = ',F11.4,' SUM(P) = ',F11.4,/,12X,' CHECK',
     +  ' INPUT VECTORS AND/OR ARRAY DIMENSIONS. ####')
        RETURN
      ENDIF
C
      IF(NJJ.LE.NJETMX) THEN
        DO 5004 J=1,NJJ
          DO 5005 I=1,7
            PJET(I,J,NJJ) =  PL(I,J)
 5005     CONTINUE
 5004   CONTINUE
      ENDIF
C
C CALCULATE AND STORE PAIR MASSES
C
       DO 5006 I=1,NJJ-1
         DO 5007 J=I+1,NJJ
           IF(IMODE.EQ.1) THEN
             Y(I,J) = JADE(I,J)
           ELSEIF(IMODE.EQ.5) THEN
             Y(I,J) = D(I,J)
           ELSEIF(IMODE.EQ.6) THEN
             Y(I,J) = G(I,J)
           ELSE
             Y(I,J) = E(I,J)
           ENDIF
 5007    CONTINUE
 5006  CONTINUE
C
C FIND LOCAL MINIMUM OF PAIR MASSES
C
      IM = 0
      JM = 0
 1000 CONTINUE
      YMINI = 2.*EVIS**2
      DO 5008 I=1,NJJ-1
        DO 5009 J=I+1,NJJ
          IF(Y(I,J).LT.YMINI) THEN
           YMINI = Y(I,J)
           IM = I
           JM = J
          ENDIF
 5009   CONTINUE
 5008 CONTINUE
C
      IF(IM.EQ.0 .OR. JM.EQ.0) THEN
        WRITE(6,284)
  284   FORMAT(' #### YKERN: ERROR; NO MINIMUM FOUND IN Y-ARRAY! ####')
        RETURN
      ENDIF
C
C RECOMBINE PARTICLES IM AND JM, STORE AT POSITION IM
C
      PL(1,IM) = PL(1,IM) + PL(1,JM)
      PL(2,IM) = PL(2,IM) + PL(2,JM)
      PL(3,IM) = PL(3,IM) + PL(3,JM)
      PL(6,IM) = SQRT(PL(1,IM)**2 + PL(2,IM)**2
     +         + PL(3,IM)**2)
      IF(IMODE.EQ.2) THEN
        PL(4,IM) = PL(6,IM)
      ELSEIF(IMODE.EQ.3) THEN
        EVISO = EVIS
        EVIS = EVIS-PL(4,IM)-PL(4,JM)+PL(6,IM)
        PL(4,IM) = PL(6,IM)
      ELSE
        PL(4,IM) = PL(4,IM) + PL(4,JM)
      ENDIF
C                  KEEP TRACK OF # OF PARTICLES ASSIGNED TO NEW CLUSTER:
      PL(7,IM) = PL(7,IM) + PL(7,JM)
C                  MOVE LAST PARTICLE IN LIST (NJJ) TO EMPTY SLOT (JM)
      NOLD = 0
      IF(JM.NE.NJJ) THEN
        DO 5010 KK=1,7
          PL(KK,JM) = PL(KK,NJJ)
 5010   CONTINUE
        NOLD = NJJ
      ENDIF
      HISTOR(1,NJJ) = JM
      HISTOR(2,NJJ) = IM
      NJJ = NJJ - 1
C                  REMEMBER JET AXES AND VALUE OF YIJ OF FLIP
      IF(NJJ.LE.NJETMX) THEN
        IF(IMODE.EQ.3) THEN
          YREC(NJJ) = YMINI/EVISO**2
        ELSEIF(IMODE.EQ.6) THEN
          YREC(NJJ) = YMINI
        ELSE
          YREC(NJJ) = YMINI/EVIS**2
        ENDIF
        DO 5011 I=1,NJJ
          DO 5012 K=1,7
            PJET(K,I,NJJ) = PL(K,I)
 5012     CONTINUE
 5011   CONTINUE
      ENDIF
C
C END IF 1-JET CASE REACHED
C
      IF(NJJ.LE.1) GOTO 9000
C
C NOW CALCULATE RELEVANT NEW MASS-COMBINATIONS
C
      DO 5013 II=1,NJJ
        IF(II.NE.IM) THEN
         I = MIN(II,IM)
         J = MAX(II,IM)
         IF(IMODE.EQ.1) THEN
           Y(I,J) = JADE(I,J)
         ELSEIF(IMODE.EQ.5) THEN
           Y(I,J) = D(I,J)
         ELSEIF(IMODE.EQ.6) THEN
           Y(I,J) = G(I,J)
         ELSE
           Y(I,J) = E(I,J)
         ENDIF
        ENDIF
        IF(NOLD.NE.0) THEN
          I = MIN(II,JM)
          J = MAX(II,JM)
          Y(I,J) = Y(II,NOLD)
        ENDIF
 5013 CONTINUE
C
C  BACK TO START OF LOOP
C
      GOTO 1000
C
 9000 CONTINUE
C
      IMODEO = IMODE
      NTO = NT
      IERR = 0
CCCCC
C     WRITE(6,825) (J,(HISTOR(I,J),I=1,2),J=1,NT)
C825  FORMAT(/,' HISTORY:',/,250(I3,4X,2I4,/))
CCCCC
      RETURN
      END
C
CDECK  ID>, YNJET.
      SUBROUTINE YNJET(YCUT,NJET,IERR)
C
C  ROUTINE TO RETURN THE NUMBER OF JETS FOR A GIVEN YCUT
C
C  LAST MOD  :  21-Jun-98
C
C  Modification Log.
C  08-Oct-97 D. Chrisman, Increase NYCLMX from 250 to 500.
C  21-Jun-98 D. Chrisman, Want to be able to handle more than 10 jets.
C                         Introduce parameter NJETMX. Use this
C                         in declaration of PJET and YREC in order to remove
C                         hard coded numbers. Also replace the hard
C                         coded number in the "D0 5001" loop with the
C                         parameter NJETMX.
C
      IMPLICIT NONE
C     IMPLICIT NONE
      INTEGER NJET,IERR,NYCLMX,NJETMX,IMODEO,NJETO,NTO,I
      PARAMETER (NYCLMX = 500, NJETMX = 20)
      REAL YCUT,YREC(NJETMX),PINT(10,NYCLMX),PJET(10,NJETMX,NJETMX)
      INTEGER HISTOR(2,NYCLMX)
      COMMON /YCL/YREC,PJET,HISTOR
      COMMON /YINT/ IMODEO,NTO,NJETO,PINT
      IERR = -1
C
CHECK IF CALL WAS MADE TO YKERN
C
      IF(IMODEO.LE.0 .OR. IMODEO.GT.7) THEN
        WRITE(6,111)
 111    FORMAT(' #### YNJET: YKERN MUST BE CALLED FIRST ! ####')
        NJET = -1
        RETURN
      ENDIF
C
CHECK IF INPUT MAKES SENSE
C
      IF(YCUT.LE.0. .OR. YCUT.GT.1.) THEN
        WRITE(6,1) YCUT
 1      FORMAT(' #### YNJET: INPUT YCUT=',E12.4,' INVALID ####')
        NJET = -1
        RETURN
      ENDIF
C
      NJET = 1
      DO 5001 I=1,NJETMX
        IF(YCUT.LT.YREC(I)) NJET = I+1
 5001 CONTINUE
      IERR = 0
      RETURN
      END
C
CDECK  ID>, YYJET.
      SUBROUTINE YYJET(NJET,YL,YH,IERR)
C
C  ROUTINE TO RETURN THE VALUES OF YCUT BETWEEN WHICH EVENT IS
C  CLASSIFIED AS N-JET (YL < YH)
C
C  LAST MOD  :  21-Jun-98
C
C  Modification Log.
C  08-Oct-97 D. Chrisman, Increase NYCLMX from 250 to 500.
C  21-Jun-98 D. Chrisman, Want to be able to handle more than 10 jets.
C                         Introduce parameter NJETMX. Use this
C                         in declaration of PJET and YREC in order to remove
C                         hard coded numbers.
C
      IMPLICIT NONE
C     IMPLICIT NONE
      INTEGER NYCLMX,NJETMX,NJET,IERR,IMODEO,NJETO,NTO
      PARAMETER (NYCLMX = 500, NJETMX = 20)
      INTEGER HISTOR(2,NYCLMX)
      REAL YREC(NJETMX),YL,YH,PINT(10,NYCLMX),PJET(10,NJETMX,NJETMX)
      COMMON /YCL/YREC,PJET,HISTOR
      COMMON /YINT/ IMODEO,NTO,NJETO,PINT
      IERR = -1
C
CHECK IF CALL WAS MADE TO YKERN
C
      IF(IMODEO.LE.0 .OR. IMODEO.GT.7) THEN
        WRITE(6,111)
 111    FORMAT(' #### YYJET: YKERN MUST BE CALLED FIRST ! ####')
        YL = -1.
        YH = -1.
        RETURN
      ENDIF
C
CHECK IF INPUT MAKES SENSE
C
      IF(NJET.LE.0 .OR. NJET.GT.NJETMX) THEN
        WRITE(6,1) NJET
 1      FORMAT(' #### YYJET: REQUEST FOR NJET=',I12,
     +  ' NOT SUPPORTED ####')
        YL = -1.
        YH = -1.
        RETURN
      ENDIF
      IF(NJET.GT.NTO) THEN
        WRITE(6,2) NJET,NTO
 2      FORMAT(' #### YYJET:',I3,' JETS OUT OF',I3,' PARTICLES NOT',
     +  ' POSSIBLE. ####')
        YL = -1.
        YH = -1.
        RETURN
      ENDIF
C
      IF(NJET.EQ.1) THEN
        YH = 1.
      ELSE
        YH = YREC(NJET-1)
      ENDIF
      YL = YREC(NJET)
C
      IERR = 0
      RETURN
      END
C
CDECK  ID>, YAXES.
      SUBROUTINE YAXES(NJET,PNJ,IERR)
C
C  ROUTINE TO RETURN THE JET AXES WHEN EVENT IS CLASSIFIED AS N-JET
C
C  LAST MOD  :  21-Jun-98
C
C  Modification Log.
C  08-Oct-97 D. Chrisman, Increase NYCLMX from 250 to 500.
C  21-Jun-98 D. Chrisman, Want to be able to handle more than 10 jets.
C                         Introduce parameter NJETMX. Use this
C                         in declaration of PJET, YREC and PNJ in order to
C                         remove hard coded numbers.
C
      IMPLICIT NONE
C     IMPLICIT NONE
      INTEGER NYCLMX,NJETMX,NJET,IERR,IMODEO,NJETO,NTO,ICHECK,I,J,K,N
      PARAMETER (NYCLMX = 500, NJETMX = 20)
      INTEGER HISTOR(2,NYCLMX)
      REAL YREC(NJETMX),PINT(10,NYCLMX),PJET(10,NJETMX,NJETMX)
      REAL PNJ(10,NJETMX)
      COMMON /YCL/YREC,PJET,HISTOR
      COMMON /YINT/ IMODEO,NTO,NJETO,PINT
      IERR = -1
C
CHECK IF CALL WAS MADE TO YKERN
C
      ICHECK = 0
      IF(IMODEO.LE.0 .OR. IMODEO.GT.7) THEN
        WRITE(6,111)
 111    FORMAT(' #### YAXES: YKERN MUST BE CALLED FIRST ! ####')
        ICHECK = -1
      ENDIF
C
CHECK IF INPUT MAKES SENSE
C
      IF(NJET.LE.0 .OR. NJET.GT.NJETMX) THEN
        WRITE(6,1) NJET
 1      FORMAT(' #### YAXES: REQUEST FOR NJET=',I4,
     +  ' NOT SUPPORTED ####')
        ICHECK = -1
      ENDIF
      IF(NJET.GT.NTO) THEN
        WRITE(6,2) NJET,NTO
 2      FORMAT(' #### YAXES:',I4,' JETS FOR',I4,' PARTICLES NOT',
     +  ' POSSIBLE. ####')
        ICHECK = -1
      ENDIF
      IF(ICHECK.NE.0) THEN
        DO 5001 I=1,10
          DO 5002 J=1,NJETMX
            PNJ(I,J) = -1.
 5002     CONTINUE
 5001   CONTINUE
        RETURN
      ENDIF
C
      DO 5003 N=1,NJET
        DO 5004 K=1,7
          PNJ(K,N) = PJET(K,N,NJET)
 5004   CONTINUE
 5003 CONTINUE
C
      IERR = 0
      RETURN
      END
C
CDECK  ID>, YASSO.
      SUBROUTINE YASSO(NJET,PNJ,BL,IERR)
C
C  ROUTINE TO RETURN THE ASSIGNMENT OF PARTICLES TO JET AXES, FOR
C  CLASSIFICATION AS N-JET (FOR PARTICLE K IN MOMENTUM ARRAY,
C  J=BL(K) POINTS TO JET NUMBER J IN ARRAY PNJ).
C  JETS ARE ORDERED ACCORDING TO THEIR ENERGIES: E1 >= E2 >= ...
C
C  LAST MOD  : 28-Jul-99
C
C  Modification Log.
C  28-Jul-99 D. Chrisman, ITAG and IREORD are now dimensioned with NJETMX 
C                         elements.
C  08-Oct-97 D. Chrisman, Increase NYCLMX from 250 to 500.
C  21-Jun-98 D. Chrisman, Want to be able to handle more than 10 jets.
C                         Introduce parameter NJETMX. Use this
C                         in declaration of PJET, YREC and PNJ in order to
C                         remove hard coded numbers.
C
      IMPLICIT NONE
C     IMPLICIT NONE
      INTEGER NYCLMX,NJETMX,NJET,IERR,IMODEO,NJETO,NTO,ICHECK,I,J,K,N      
      INTEGER I1,I2
      PARAMETER (NYCLMX = 500, NJETMX = 20)
      INTEGER HISTOR(2,NYCLMX),IMAX,BL(NYCLMX),ITAG(NJETMX),
     +         IREORD(NJETMX)
      REAL YREC(NJETMX),PINT(10,NYCLMX),PJET(10,NJETMX,NJETMX)    
      REAL PNJ(10,NJETMX),EMAX
      COMMON /YCL/YREC,PJET,HISTOR
      COMMON /YINT/ IMODEO,NTO,NJETO,PINT
      IERR = -1
C
CHECK IF CALL WAS MADE TO YKERN
C
      ICHECK = 0
      IF(IMODEO.LE.0 .OR. IMODEO.GT.7) THEN
        WRITE(6,111)
 111    FORMAT(' #### YASSO: YKERN MUST BE CALLED FIRST ! ####')
        ICHECK = -1
      ENDIF
C
CHECK IF INPUT MAKES SENSE
C
      IF(NJET.LE.0 .OR. NJET.GT.NJETMX) THEN
        WRITE(6,1) NJET
 1      FORMAT(' #### YASSO: REQUEST FOR NJET=',I4,
     +  ' NOT SUPPORTED ####')
        ICHECK = -1
      ENDIF
      IF(NJET.GT.NTO) THEN
        WRITE(6,2) NJET,NTO
 2      FORMAT(' #### YASSO:',I4,' JETS OUT OF',I4,' PARTICLES NOT',
     +  ' POSSIBLE. ####')
        ICHECK = -1
      ENDIF
      IF(ICHECK.NE.0) THEN
        DO 5001 I=1,10
          DO 5002 J=1,NJETMX
            PNJ(I,J) = -1.
 5002     CONTINUE
 5001   CONTINUE
        DO 5003 I=1,NTO
          BL(I) = -1
 5003   CONTINUE
        RETURN
      ENDIF
C
      DO 5004 I=1,NTO
        BL(I) = I
 5004 CONTINUE
C
      IF(NJET.NE.NTO) THEN
      DO 5005 I=NTO,NJET+1,-1
        I1 = HISTOR(1,I)
        I2 = HISTOR(2,I)
        DO 5006 N=1,NTO
          IF(BL(N).EQ.I1) BL(N) = I2
          IF(I1.NE.I) THEN
            IF(BL(N).EQ.I) BL(N) = I1
          ENDIF
 5006   CONTINUE
 5005 CONTINUE
      ENDIF
C
C ORDER JETS ACCORDING TO THEIR ENERGY (FIRST IS LARGEST)
C AND CHANGE POINTERS ACCORDINGLY
C
      DO 5007 I=1,NJET
        ITAG(I) = 1
 5007 CONTINUE
      DO 5008 I=1,NJET
C       IF(ITAG(I).NE.0) THEN
          EMAX = 0.
          IMAX = 0
          DO 5009 J=1,NJET
            IF(ITAG(J).NE.0 .AND. EMAX.LT.PJET(4,J,NJET)) THEN
               EMAX = PJET(4,J,NJET)
               IMAX = J
            ENDIF
 5009     CONTINUE
          IF(IMAX.LE.0) THEN
            WRITE(6,9)
 9          FORMAT(' #### YASSO: JET AXIS WITH ZERO OR NEGATIVE ',
     +      'ENERGY COMPONENT DETECTED; NO ORDERING DONE. ####')
            DO 5010 N=1,NJET
              DO 5011 K=1,7
                PNJ(K,N) = PJET(K,N,NJET)
 5011         CONTINUE
 5010       CONTINUE
            RETURN
          ENDIF
          ITAG(IMAX) = 0
          IREORD(IMAX) = I
C       ENDIF
 5008 CONTINUE
C
      DO 5012 N=1,NJET
        DO 5013 K=1,7
          PNJ(K,IREORD(N)) = PJET(K,N,NJET)
 5013   CONTINUE
 5012 CONTINUE
C
      DO 5014 I=1,NTO
        BL(I) = IREORD(BL(I))
 5014 CONTINUE
C
      IERR = 0
      RETURN
      END
C
CDECK  ID>, YREAS.
      SUBROUTINE YREAS(NJET,PNJ,YMIN,BL,IERR)
C
C  REASSIGNES PARTICLES TO CLOSEST JET (IN ANGLE);
C  STARTS FROM INITIAL EVENT CONFIGURATION ACCORDING TO LAST CALL
C  OF YKERN; RE-IERATES NEW JET AXES AND PARTICLE ASSIGNMENTS UNTIL
C  STABLE CONFIGURATION IS REACHED. IF AN AXES WITH ZERO PARTICLES
C  AND MOMENTUM EMERGES, ROUTINE RETAINS THE INITIAL N-JET CONFIGURA-
C  TION AND RETURNS IERR=-2. IERR=-1 MEANS INVALID INPUT CONDITIONS;
C  RESULTS MUST BE SKIPPED.
C  INPUT:   NJET IS NUMBER OF JETS REQUIRED
C  OUTPUTS: PNJ(K,NJ) CONTAINS NEW JET AXES; POS. K=7 HOLDS NUMBER OF
C                     PARTICLES ASSIGNED TO AXES NJ.
C           YMIN IS VALUE OF JET RESOLUTION YCUT WHERE EVENT FLIPS
C                TO (NJET-1) JET CONFIGURATION [COMPARE WITH
C                YREC(NJET-1) OF CONFIGURATION BEFORE REASSIGNMENT].
C           BL(N) CONTAINS POINTER TO JET AXES TO WHICH PARTICLE N
C                 IS ASSIGNED.
C
C  LAST MOD  :  21-Jun-98
C
C  Modification Log.
C  12-May-97 D. Chrisman, remove declaration of unused variable EMAX.
C  08-Oct-97 D. Chrisman, Increase NYCLMX from 250 to 500.
C  21-Jun-98 D. Chrisman, Want to be able to handle more than 10 jets.
C                         Introduce parameter NJETMX. Use this
C                         in declaration of PJET, YREC, PNJ and PNN in order to
C                         remove hard coded numbers.
C
      IMPLICIT NONE
C     IMPLICIT NONE
      INTEGER NYCLMX,NJETMX,NJET,IERR,IMODEO,NJETO,NTO,ICHECK,I,J    
      INTEGER NITMAX,IER
      PARAMETER (NYCLMX = 500, NJETMX = 20)
      PARAMETER (NITMAX=10)
      INTEGER HISTOR(2,NYCLMX),NIT,ICFLAG,JMAX,K,NZERO,BL(NYCLMX)
      REAL YREC(NJETMX),PINT(10,NYCLMX),PJET(10,NJETMX,NJETMX)
      REAL PNJ(10,NJETMX)
      REAL PNN(10,NJETMX),COSMAX,COSIJ,EVIS,YIJ,YMIN
      COMMON /YCL/YREC,PJET,HISTOR
      COMMON /YINT/ IMODEO,NTO,NJETO,PINT
      IERR = -1
C
CHECK IF CALL WAS MADE TO YKERN
C
      ICHECK = 0
      IF(IMODEO.LE.0 .OR. IMODEO.GT.7) THEN
        WRITE(6,111)
 111    FORMAT(' #### YREAS: YKERN MUST BE CALLED FIRST ! ####')
        ICHECK = -1
      ENDIF
C
CHECK IF INPUT MAKES SENSE
C
      IF(NJET.LE.0 .OR. NJET.GT.NJETMX) THEN
        WRITE(6,1) NJET
 1      FORMAT(' #### YREAS: REQUEST FOR NJET=',I4,
     +  ' NOT SUPPORTED ####')
        ICHECK = -1
      ENDIF
      IF(NJET.GT.NTO) THEN
        WRITE(6,2) NJET,NTO
 2      FORMAT(' #### YREAS:',I4,' JETS OUT OF',I4,' PARTICLES NOT',
     +  ' POSSIBLE. ###')
        ICHECK = -1
      ENDIF
      IF(ICHECK.NE.0) THEN
        DO 5001 I=1,10
          DO 5002 J=1,NJETMX
            PNJ(I,J) = -1.
 5002     CONTINUE
 5001   CONTINUE
        DO 5003 I=1,250
          BL(I) = -1
 5003   CONTINUE
        RETURN
      ENDIF
C
C GET INITIAL JET AXES FROM OUTPUT OF JET FINDER
C
      CALL YASSO(NJET,PNJ,BL,IER)
C
      NIT = 0
C
  100 CONTINUE
      NIT = NIT + 1
      DO 5004 J=1,NJET
        DO 5005 I=1,7
          PNN(I,J) = 0.
 5005   CONTINUE
 5004 CONTINUE
      ICFLAG = 0
      DO 5006 I=1,NTO
        COSMAX = -2.
        DO 5007 J=1,NJET
          COSIJ = (PINT(1,I)*PNJ(1,J)+PINT(2,I)*PNJ(2,J)+
     +    PINT(3,I)*PNJ(3,J)) / (PINT(6,I)*PNJ(6,J))
          IF(COSIJ.GT.COSMAX) THEN
            COSMAX = COSIJ
            JMAX = J
          ENDIF
 5007   CONTINUE
C ADD PARTICLE I TO NEW JET JMAX
        DO 5008 K=1,4
          PNN(K,JMAX) = PNN(K,JMAX) + PINT(K,I)
 5008   CONTINUE
        PNN(7,JMAX) = PNN(7,JMAX) + 1.
CHECK WHETHER ASSIGNMENT HAS CHANGED
        IF(BL(I). NE. JMAX) THEN
          ICFLAG = ICFLAG + 1
          BL(I) = JMAX
        ENDIF
 5006 CONTINUE
C
CHECK IF JET WITH ZERO PARTICLES OCCURED
C
      NZERO = 0
      DO 5009 J=1,NJET
        IF(PNN(7,J).LT.1.) NZERO = NZERO + 1
 5009 CONTINUE
C
      IF(NZERO.NE.0) THEN
        WRITE(6,3)
 3      FORMAT(' #### YREAS: JET WITH ZERO PARTICLES OCCURED;',
     +  ' RETAIN ORIGINAL AXES. ####')
        CALL YASSO(NJET,PNJ,BL,IER)
        IERR = -2
        ICFLAG = 0
        YMIN = YREC(NJET-1)
      ELSE
COPY NEW JET MOMENTA TO OUTPUT ARRAY
        EVIS = 0.
        DO 5010 J=1,NJET
          PNN(6,J)=SQRT(PNN(1,J)**2 + PNN(2,J)**2 + PNN(3,J)**2)
          PNJ(1,J) = PNN(1,J)
          PNJ(2,J) = PNN(2,J)
          PNJ(3,J) = PNN(3,J)
          PNJ(6,J) = PNN(6,J)
          PNJ(7,J) = PNN(7,J)
          IF(IMODEO.EQ.2 .OR. IMODEO.EQ.3) THEN
            PNJ(4,J) = PNN(6,J)
          ELSE
            PNJ(4,J) = PNN(4,J)
          ENDIF
          EVIS = EVIS + PNJ(4,J)
 5010   CONTINUE
      ENDIF
C  RE-ITERATE IF A CHANGE OCCURED
      IF(ICFLAG.NE.0 .AND. NIT.LT.NITMAX) GOTO 100
C  END OF ITERATION
      IF(NIT.GE.NITMAX .AND. ICFLAG.NE.0) WRITE(6,200) NIT,ICFLAG
 200  FORMAT(' #### YREAS: JET ASSIGNMENT NOT STABLE AFTER',I3,
     +' ITERATIONS (STILL',I3,' CHANGES) ####')
C
CALCULATE FINAL VALUE OF JET RESOLUTION FOR TRANSITION TO NJET-1 JETS
C
      IF(IERR.NE.-2) THEN
      YMIN = EVIS**2
      DO 5011 I=1,NJET-1
        DO 5012 J=I+1,NJET
          IF(IMODEO.EQ.1) THEN
            YIJ = 2.*PNJ(4,I)*PNJ(4,J)*MAX(0.,(1.-
     +      (PNJ(1,I)*PNJ(1,J)+PNJ(2,I)*PNJ(2,J)+PNJ(3,I)*PNJ(3,J))/
     +      (PNJ(6,I)*PNJ(6,J))))
          ELSEIF(IMODEO.EQ.5) THEN
            YIJ = 2.*MIN(PNJ(4,I)*PNJ(4,I),PNJ(4,J)*PNJ(4,J))*MAX(0.,
     +      (1.-(PNJ(1,I)*PNJ(1,J)+PNJ(2,I)*PNJ(2,J)+PNJ(3,I)*PNJ(3,J))/
     +      (PNJ(6,I)*PNJ(6,J))))
          ELSEIF(IMODEO.EQ.6) THEN
            YIJ = 8.*PNJ(4,I)*PNJ(4,J)*MAX(0.,(1.-
     +      (PNJ(1,I)*PNJ(1,J)+PNJ(2,I)*PNJ(2,J)+PNJ(3,I)*PNJ(3,J))/
     +      (PNJ(6,I)*PNJ(6,J))))/(9.*(PNJ(4,I)+PNJ(4,J))**2)
          ELSE
            YIJ = MAX(0.,(PNJ(4,I)+PNJ(4,J))**2-(PNJ(1,I)+PNJ(1,J))**2-
     +      (PNJ(2,I)+PNJ(2,J))**2-(PNJ(3,I)+PNJ(3,J))**2)
          ENDIF
          IF(YIJ.LT.YMIN) YMIN = YIJ
 5012   CONTINUE
 5011 CONTINUE
      IF(IMODEO.NE.6) YMIN = YMIN / EVIS**2
      IERR = 0
      ENDIF
C
      RETURN
      END
C
CDECK  ID>, YTREE.
      SUBROUTINE YTREE(LPRINT,PTR,IERR)
C
C  ROUTINE TO GENERATE, RETURN AND (IF REQUESTED) PRINT OUT THE
C  ENTIRE EVENT HISTORY OF PARTICLE/JET RECOMBINATIONS, FOR THE
C  EVENT AND PARTICLES FROM LAST CALL TO YKERN.
C  INPUT:  LPRINT = .TRUE.: EVENT HISTORY (I.E. PTR) SHALL BE PRINTED
C  OUTPUT: PTR(10,*) CONTAINS PARTICLE/JET 4-MOMENTA, MASSES, MOMENTA
C                    AND HISTORY INFORMATION.
C          PTR(1-6,K):  P_X, P_Y, P_Z, E, MASS, |P| OF VECTOR K
C          PTR(7,*)  :  VALUE OF Y_IJ WHEN SPLIT INTO DAUGHTERS (IF ANY)
C          PTR(8,*)  :  POINTER TO PARENT JET/CLUSTER (0 FOR 1ST VECTOR)
C          PTR(9-10,*): POINTERS TO DAUGHTERS (0 FOR FINAL STATE PART.).
C          VECTORS K= 1 TO NT-1 (NT IS THE NUMBER OF PARTICLES FOR WHICH
C          YKERN WAS CALLED) CORRESPOND TO THE CLUSTERS OR JETS FROM THE
C          ENTIRE RECOMBINATION CHAIN, DOWN TO 1-JET; THE 1ST VECTOR
C          THUS HOLDS THE SUM OF 4-VECTORS OF ALL FINAL STATE PARTICLES.
C          VECTORS NT TO 2*NT-1 ARE AN EXACT COPY OF THE INPUT PARTICLES
C          TO PREVIOUS CALL OF YKERN.
C          NOTE: E AND Y_IJ CONFORM WITH THE JET SCHEME CHOSEN WHEN
C          YKERN WAS CALLED; E.G. E=|P| IN P- AND P0-SCHEME.
C
C  LAST MOD  :  21-Jun-98
C
C  Modification Log.
C  12-May-97 D. Chrisman, remove declaration of unused variable NJET.
C  08-Oct-97 D. Chrisman, Increase NYCLMX from 250 to 500.
C  21-Jun-98 D. Chrisman, Want to be able to handle more than 10 jets.
C                         Introduce parameter NJETMX. Use this
C                         in declaration of PJET and YREC in order to
C                         remove hard coded numbers.
C
      IMPLICIT NONE
C     IMPLICIT NONE
      INTEGER NYCLMX,NJETMX,IERR,IMODEO,NJETO,NTO,ICHECK,I,J,K,I1,I2,II
      PARAMETER (NYCLMX = 500, NJETMX = 20)
      INTEGER HISTOR(2,NYCLMX),BL(NYCLMX)
      REAL YREC(NJETMX),PINT(10,NYCLMX),PJET(10,NJETMX,NJETMX),EVIS
      COMMON /YCL/YREC,PJET,HISTOR
      COMMON /YINT/ IMODEO,NTO,NJETO,PINT
      REAL PTR(10,NYCLMX)
      LOGICAL LPRINT
      IERR = -1

Cf2py intent(out) PTR
Cf2py intent(out) IERR
      
C
CHECK IF CALL WAS MADE TO YKERN
C
      ICHECK = 0
      IF(IMODEO.LE.0 .OR. IMODEO.GT.7) THEN
        WRITE(6,111)
 111    FORMAT(' #### YTREE:  YKERN MUST BE CALLED FIRST ! ####')
        ICHECK = -1
      ENDIF
C
CHECK IF INPUT MAKES SENSE
C
      IF(NYCLMX.LT.2*NTO-1) THEN
        WRITE(6,1) NYCLMX,2*NTO
 1      FORMAT(' #### YTREE: NOT ENOUGH SPACE FOR EVENT TREE; ',
     +  'INCREASE ARRAY DIMENSIONS FROM',I4,' TO',I4,' ####')
        ICHECK = -1
      ENDIF
      IF(ICHECK.NE.0) THEN
        DO 5001 I=1,10
          DO 5002 J=1,NYCLMX
            PTR(I,J) = -1.
 5002     CONTINUE
 5001   CONTINUE
        RETURN
      ENDIF
C
      DO 5003 I=1,NTO
        BL(I) = NTO+I
 5003 CONTINUE
      EVIS = 0.
C
C  COPY INPUT MOMENTA AND INITIALIZE OUTPUT ARRAY
C
      DO 5004 I=1,NTO
        K = NTO + I -1
        PTR(1,K) = PINT(1,I)
        PTR(2,K) = PINT(2,I)
        PTR(3,K) = PINT(3,I)
        PTR(6,K) = SQRT(PINT(1,I)**2+PINT(2,I)**2+PINT(3,I)**2)
        IF(IMODEO.EQ.2 .OR. IMODEO.EQ.3) THEN
          PTR(4,K) = PTR(6,K)
          PTR(5,K) = 0.
        ELSE
          PTR(4,K) = PINT(4,I)
          PTR(5,K) = SQRT(MAX(0.,PTR(4,K)**2-PTR(6,K)**2))
        ENDIF
        PTR(7,K) = 0.
        PTR(8,K) = 0.
        PTR(9,K) = 0.
        PTR(10,K) = 0.
        EVIS = EVIS + PTR(4,K)
 5004 CONTINUE
C
C RECONSTRUCT AND STORE EVENT TREE
C
      DO 5005 I=NTO,2,-1
        I1 =  BL(HISTOR(1,I))-1
        I2 =  BL(HISTOR(2,I))-1
        II = I-1
        PTR(1,II) = PTR(1,I1) + PTR(1,I2)
        PTR(2,II) = PTR(2,I1) + PTR(2,I2)
        PTR(3,II) = PTR(3,I1) + PTR(3,I2)
        PTR(6,II) = SQRT(PTR(1,II)**2+PTR(2,II)**2+PTR(3,II)**2)
        IF(IMODEO.EQ.2 .OR. IMODEO.EQ.3) THEN
          PTR(4,II) = PTR(6,II)
          PTR(5,II) = 0.
        ELSE
          PTR(4,II) = PTR(4,I1) + PTR(4,I2)
          PTR(5,II) = SQRT(MAX(0.,PTR(4,II)**2-PTR(6,II)**2))
        ENDIF
C FLAG THIS VECTOR (I) AS PARENT OF DAUGHTER VECTORS (BIGGEST FIRST)
        IF(PTR(4,I1).GT.PTR(4,I2)) THEN
          PTR(9,II) = FLOAT(I1)
          PTR(10,II)= FLOAT(I2)
        ELSE
          PTR(10,II)= FLOAT(I1)
          PTR(9,II) = FLOAT(I2)
        ENDIF
C FLAG DAUGHTERS OF THIS VECTOR (I)
        PTR(8,I1)= FLOAT(II)
        PTR(8,I2)= FLOAT(II)
C UPDATE RELATIVE POSITION OF VECTORS
        BL(HISTOR(2,I)) = I
        BL(HISTOR(1,I)) = BL(I)
C CALC AND STORE RESOLUTION PARAMETER Y_IJ FOR PARENT -> DAUGHTERS I,J
        IF(IMODEO.EQ.1) THEN
         PTR(7,II) = 2.*PTR(4,I1)*PTR(4,I2)*MAX(0.,(1.-
     +   (PTR(1,I1)*PTR(1,I2)+PTR(2,I1)*PTR(2,I2)+PTR(3,I1)*PTR(3,I2))/
     +   (PTR(6,I1)*PTR(6,I2)))) / EVIS**2
        ELSEIF(IMODEO.EQ.5) THEN
         PTR(7,II) = 2.*MIN(PTR(4,I1)*PTR(4,I1),PTR(4,I2)*PTR(4,I2))*
     +   MAX(0.,(1.-(PTR(1,I1)*PTR(1,I2)+PTR(2,I1)*PTR(2,I2)+PTR(3,I1)*
     +   PTR(3,I2))/(PTR(6,I1)*PTR(6,I2)))) / EVIS**2
        ELSEIF(IMODEO.EQ.6) THEN
         PTR(7,II) = 8.*PTR(4,I1)*PTR(4,I2)*MAX(0.,(1.-
     +   (PTR(1,I1)*PTR(1,I2)+PTR(2,I1)*PTR(2,I2)+PTR(3,I1)*PTR(3,I2))/
     +   (PTR(6,I1)*PTR(6,I2))))/(9.*(PTR(4,I1)+PTR(4,I2))**2)
        ELSE
         PTR(7,II) = MAX(0.,(PTR(4,I1)+PTR(4,I2))**2-(PTR(1,I1)+
     +   PTR(1,I2))**2-(PTR(2,I1)+PTR(2,I2))**2-
     +   (PTR(3,I1)+PTR(3,I2))**2) / EVIS**2
        ENDIF
C READJUST EVIS IF IN MODE 3 (P0)
        IF(IMODEO.EQ.3) THEN
          EVIS = EVIS - PTR(4,I1) - PTR(4,I2) + PTR(4,II)
        ENDIF
 5005 CONTINUE
C
      IF(LPRINT) THEN
C       WRITE(6,201) (I,(HISTOR(J,I),J=1,2),I=1,NTO)
C201    FORMAT(/,' HISTOR:',/,250(I4,2X,2I5/))
        WRITE(6,199)
 199    FORMAT(/,' YTREE: EVENT TREE INFORMATION',/,
     +  ' POS      P_X      P_Y      P_Z       E        M       |P| ',
     +  ' Y_D1D2 PAR  D1  D2')
        DO 5010 I=1,NTO-1
         WRITE(6,200) I,(PTR(J,I),J=1,7),(INT(PTR(J,I)),J=8,10)
 200     FORMAT(I4,6(1X,F8.3),1X,F7.4,3(1X,I3))
 5010   CONTINUE
        WRITE(6,202)
 202    FORMAT('        STABLE PARTICLES:')
        DO 5011 I=NTO,2*NTO-1
          WRITE(6,203) I,(PTR(J,I),J=1,7),INT(PTR(8,I))
 203      FORMAT(I4,6(1X,F8.3),1X,F7.4,1X,I3)
 5011   CONTINUE
      ENDIF
C
      IERR = 0
      RETURN
      END
