       SUBROUTINE CMPSCORE(TIME,STAGE,COV,GARRAY,AMAT,AINV,EVS,LEN,SCORE
     c       ,INFO,THTA,NCOV,NPAR,NSTAGE)
        INTEGER LEN,NPAR,STAGE(LEN),CT1,IND1,IND2,CTR,NSTAGE
        INTEGER NCOV,COV(LEN),CI
        REAL*8 AMAT(NCOV,NSTAGE,NSTAGE),AINV(NCOV,NSTAGE,NSTAGE)
        REAL*8 INFO(NPAR,NPAR),EVS(NCOV,NSTAGE),TMP,T1,EEVS(50)
        REAL*8 PPARR(50,50,50),GARRAY(NCOV,NPAR,NSTAGE,NSTAGE)
        REAL*8 THTA(NPAR),SCORE(NPAR),PMAT(50,50),TIME(LEN)
C Copyright 1994 Robert Gentleman
C to compute the score we loop through the observations on an individual
C and compute it step by step, at each step we need to know which prob.
C is appropriate (CI determines this) and then compute the time step
C 
C        CALL INTPR("LEN",3,LEN,1)
C        CALL INTPR("NCOV",4,NCOV,1)
C        CALL INTPR("NPAR",4,NPAR,1)
C        CALL INTPR("NSTAGE",5,NSTAGE,1)


         DO 13  CTR = 1, (LEN-1)
           CI = COV(CTR)
           T1 = TIME(CTR+1)-TIME(CTR)
           IND1 = STAGE(CTR)
           IND2 = STAGE(CTR+1)
           DO 131 I=1,NSTAGE
 131        EEVS(I)=DEXP(T1*EVS(CI,I))
            CALL COMPPP(PPARR,GARRAY,AMAT,AINV,EVS,EEVS,NSTAGE,
     c          NPAR,T1,CI,NCOV)
            DO 15 J1=1,NSTAGE
             PMAT(IND1,J1)=0
             DO 16 CT1=1,NSTAGE
              PMAT(IND1,J1)=PMAT(IND1,J1)+
     c         EEVS(CT1)*AMAT(CI,IND1,CT1)*AINV(CI,CT1,J1)
 16          CONTINUE
 15         CONTINUE
           DO 17 I=1,NPAR
            SCORE(I)=SCORE(I)+PPARR(I,IND1,IND2)/PMAT(IND1,IND2)
            DO 18 J=1,NPAR
             TMP=0
             DO 19 I1=1,NSTAGE
              IF (PMAT(IND1,I1) .GT. 0) THEN
               TMP=TMP+(PPARR(I,IND1,I1)*PPARR(J,IND1,I1))/PMAT(IND1,I1)   
              END IF
 19          CONTINUE
             INFO(I,J)=INFO(I,J)+TMP
 18         CONTINUE
 17        CONTINUE
c        call dblepr('score each time',-1,score,6)
c        call dblepr('info each time',-1,info,36)
 13       CONTINUE
c         what the heck is this doing?
         DO 98 I=1,NPAR
          THTA(I)=PMAT(1,I)
 98      CONTINUE
c  I think the intervening lines need to be chopped!
          RETURN
          END



       SUBROUTINE COMPPP(PPARR,GARRAY,AMAT,AINV,EVS,EEVS,
     c     NSTAGE,NPAR,T1,CI,NCOV)
         INTEGER NSTAGE,NPAR,CT1,CI
         REAL*8 GARRAY(NCOV,NPAR,NSTAGE,NSTAGE),AMAT(NCOV,NSTAGE,NSTAGE)
         REAL*8 AINV(NCOV,NSTAGE,NSTAGE),EVS(NCOV,NSTAGE),T1
         REAL*8 TMPA(50,50),EEVS(50),TMP1(50,50),PPARR(50,50,50)
         REAL*8 COMPTOL

c        CALL INTPR("CI",2,CI,1)
c        CALL INTPR("NPAR",4,NPAR,1)

         COMPTOL=1e-14
         DO 102 I=1,NSTAGE
          DO 103 J=1,NSTAGE
           IF (DABS(EVS(CI,I)-EVS(CI,J)) .GT. COMPTOL) THEN
            TMPA(I,J)= ((EEVS(I)-EEVS(J))/(EVS(CI,I)-EVS(CI,J)))
           ELSE
            TMPA(I,J)= T1*EEVS(I)
           END IF
 103      CONTINUE
 102     CONTINUE
         DO 104 I=1,NPAR
          DO 105 I1=1,NSTAGE
           DO 106 I2=1,NSTAGE
            TMP1(I1,I2)=0
            DO 107 CT1=1,NSTAGE
              TMP1(I1,I2)=TMP1(I1,I2)
     c         +AMAT(CI,I1,CT1)*TMPA(CT1,I2)*GARRAY(CI,I,CT1,I2)
 107        CONTINUE
 106       CONTINUE
 105      CONTINUE
          DO 108 I1=1,NSTAGE
           DO 109 I2=1,NSTAGE
            PPARR(I,I1,I2)=0
            DO 110 CT1=1,NSTAGE
             PPARR(I,I1,I2)=PPARR(I,I1,I2)+TMP1(I1,CT1)*AINV(CI,CT1,I2)
 110        CONTINUE
 109       CONTINUE
 108      CONTINUE
 104     CONTINUE
         RETURN
         END
