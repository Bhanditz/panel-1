C	Copyright 1993 by Robert Gentleman and Jane Lindsey
C		All Rights Reserved
        SUBROUTINE CMPGAMMA(TIME,STAGE,GARRAY,AMAT,AINV,AEVS,
     C    BEVS,QUV,LEN,SCORE,INFOTEST,NPAR,NPARG,NSTAGE,C1,C2)
        INTEGER LEN,NPARG,STAGE(LEN),CT1,IND1,IND2,CTR,NSTAGE
        INTEGER NPAR
        REAL*8 TIME(LEN),AMAT(NSTAGE,NSTAGE),AINV(NSTAGE,NSTAGE)
        REAL*8 AEVS(NSTAGE),SCORE(NPARG),T1
        REAL*8 PPARR(50,50,50),PMAT(50,50)
	REAL*8 GARRAY(NPAR,NSTAGE,NSTAGE),AEEVS(50)
	REAL*8 C1(2,NSTAGE,NSTAGE),C2(2,NSTAGE,NSTAGE),TMP2(50,50)
        REAL*8 E(2,50,50),TMP1(2,50,50),TMPA(50,50)
	REAL*8 TM1,TM1SQ,TM2,TM2SQ,BEVS(NSTAGE)
	REAL*8 INFOTEST(2,(NPAR+1)),BEEVS(50),TMPL,TMPG
        REAL*8 QUV,COMPTOL
C
C        CALCULATE SCORE VECTOR AND INFORMATION MATRIX FOR THETA
C
         COMPTOL=1e-14
         DO 13  CTR = 1, (LEN-1)
           TM1 = TIME(CTR)
           TM1SQ = TM1*TM1
           TM2= TIME(CTR+1)
           TM2SQ = TM2*TM2
           T1 = TM2-TM1
           IND1 = STAGE(CTR)
           IND2 = STAGE(CTR+1)
           DO 131 I=1,NSTAGE
             AEEVS(I)=DEXP(T1*AEVS(I))
 131         BEEVS(I)=DEXP(T1*BEVS(I))
C
C        CALCULATE PPARR(1:NPAR,NSTAGE,NSTAGE)
C
           DO 402 I=1,NSTAGE
             DO 403 J=1,NSTAGE
               IF (DABS(AEVS(I)-AEVS(J)) .GT. COMPTOL) THEN
                 TMPA(I,J)= ((AEEVS(I)-AEEVS(J))/(AEVS(I)-AEVS(J)))
               ELSE
                 TMPA(I,J)= T1*AEEVS(I)
               END IF
 403           CONTINUE
 402         CONTINUE
C
           DO 404 I=1,NPAR
             DO 405 I1=1,NSTAGE
               DO 406 I2=1,NSTAGE
                 TMP2(I1,I2)=0.0
                 DO 407 CT1=1,NSTAGE
                   TMP2(I1,I2)=TMP2(I1,I2)
     C               +AMAT(I1,CT1)*TMPA(CT1,I2)*GARRAY(I,CT1,I2)
 407               CONTINUE
 406             CONTINUE
 405           CONTINUE
C
             DO 408 I1=1,NSTAGE
               DO 409 I2=1,NSTAGE
                 PPARR(I,I1,I2)=0.0
                 DO 410 CT1=1,NSTAGE
                   PPARR(I,I1,I2)=PPARR(I,I1,I2)+TMP2(I1,CT1)
     C                   *AINV(CT1,I2)
 410               CONTINUE
 409             CONTINUE
 408           CONTINUE
 404         CONTINUE
C
C        CALCULATE PMAT
C
           DO 15 J1=1,NSTAGE
             PMAT(IND1,J1)=0.0
             DO 16 CT1=1,NSTAGE
               PMAT(IND1,J1)=PMAT(IND1,J1)+AEEVS(CT1)
     C               *AMAT(IND1,CT1)*AINV(CT1,J1)
 16            CONTINUE
 15          CONTINUE
C
C        CALCULATE COMPONENT OF SCORE VECTOR AND INFORMATION MATRIX
C        FOR TWO ALTERNATE HYPOTHESES
C
C        C1 C2 AND E MATRICES
C
           DO 201 I=1,NSTAGE
	     DO 202 J=1,NSTAGE
               IF (DABS(BEVS(I)-BEVS(J)) .LT. COMPTOL) THEN
                 E(2,I,J) = BEEVS(I)/2 * (TM2SQ - TM1SQ)
               ELSE
                 E(2,I,J) = (((BEVS(I)-BEVS(J))*TM2 - 1) 
     C                        * BEEVS(I)
     C                   - ((BEVS(I)-BEVS(J))*TM1-1) 
     C                          * BEEVS(J))
     C                      /(BEVS(I)-BEVS(J))**2
                 END IF
		IF (DABS(AEVS(I)-AEVS(J)) .LT. COMPTOL) THEN
                  E(1,I,J) = AEEVS(I)/2 * (TM2SQ - TM1SQ)
		ELSE
                  E(1,I,J) = (((AEVS(I)-AEVS(J))*TM2 - 1) * AEEVS(I)
     C                   - ((AEVS(I)-AEVS(J))*TM1-1)  * AEEVS(J))
     C                      /(AEVS(I)-AEVS(J))**2
                 ENDIF
 202           CONTINUE
 201       CONTINUE
C
C        CALCULATE PPARR(NPAR+1:NPAR+2,,) = C1*E*C2
C
           DO 203 I=1,2
             DO 204 I1=1,NSTAGE
               DO 205 I2=1,NSTAGE
                  TMP1(I,I1,I2)=0.0
                  DO 206 CT1=1,NSTAGE
                    TMP1(I,I1,I2) =TMP1(I,I1,I2)+
     C                  C1(I,I1,CT1)*E(I,CT1,I2)
 206                    CONTINUE
 205                CONTINUE
 204           CONTINUE
C
           DO 207 I1=1,NSTAGE
             DO 208 I2=1,NSTAGE
               PPARR((NPAR+I),I1,I2)=0.0
               DO 209 CT1=1,NSTAGE
                 PPARR((NPAR+I),I1,I2) = PPARR((NPAR+I),I1,I2) + 
     C                     TMP1(I,I1,CT1) * C2(I,CT1,I2)
 209             CONTINUE
                 IF (I .GT. 1) THEN
                   PPARR((NPAR+I),I1,I2) = PPARR((NPAR+I),I1,I2)*QUV
                   END IF
 208           CONTINUE
 207         CONTINUE
 203       CONTINUE
C
C         CALCULATE SCORE(NPAR+1:NPAR+2)
C
           DO 210 I=1,2
             SCORE(NPAR+I) = SCORE(NPAR+I) + PPARR((NPAR+I),IND1,IND2)
     C                        /PMAT(IND1,IND2)
 210         CONTINUE
C
C         CALCULATE INFORMATION TEST MATRIX
C
           DO 211 I=1,NPAR
             TMPL=0.0
             TMPG=0.0
             DO 212 I1=1,NSTAGE
               IF (PMAT(IND1,I1) .GT. 0.0) THEN
                 TMPL = TMPL + (PPARR(I,IND1,I1)*PPARR(NPAR+1,IND1,I1))/
     C                         PMAT(IND1,I1)
                 TMPG = TMPG + (PPARR(I,IND1,I1)*PPARR(NPAR+2,IND1,I1))/
     C                         PMAT(IND1,I1)
                END IF
 212          CONTINUE
            INFOTEST(1,I) = INFOTEST(1,I) + TMPL
            INFOTEST(2,I) = INFOTEST(2,I) + TMPG
 211        CONTINUE
C
C         CALCULATE INFORMATION(NPARG,NPARG)
C
           TMPL=0.0
           TMPG=0.0
           DO 213 I1=1,NSTAGE
             IF (PMAT(IND1,I1) .GT. 0.0) THEN
               TMPL = TMPL + (PPARR(NPAR+1,IND1,I1)*
     C            PPARR(NPAR+1,IND1,I1))/ PMAT(IND1,I1)
               TMPG = TMPG + (PPARR(NPAR+2,IND1,I1)*
     C            PPARR(NPAR+2,IND1,I1))/PMAT(IND1,I1)
               END IF
 213         CONTINUE
          INFOTEST(1,NPAR+1) = INFOTEST(1,NPAR+1) + TMPL
          INFOTEST(2,NPAR+1) = INFOTEST(2,NPAR+1) + TMPG
C         CALL DBLEPR("F SCORE",7,SCORE,7)
C         CALL DBLEPR("F INFO",6,INFO,49)
 13       CONTINUE
          RETURN
          END
