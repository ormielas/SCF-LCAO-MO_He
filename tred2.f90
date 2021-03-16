SUBROUTINE TRED2(A,N,NP,D,E)

  IMPLICIT NONE

  INTEGER                         :: N, NP
  REAL(KIND=8), DIMENSION(NP, NP) :: A
  REAL(KIND=8), DIMENSION(NP)     :: D, E
  
  INTEGER      :: I, J, K, L
  REAL(KIND=8) :: F, G, H, HH, SCALE

  IF(N.GT.1)THEN
     DO 18 I=N,2,-1  
        L=I-1
        H=0.
        SCALE=0.
        IF(L.GT.1)THEN
           DO 11 K=1,L
              SCALE=SCALE+ABS(A(I,K))
11         END DO
           IF(SCALE.EQ.0.)THEN
              E(I)=A(I,L)
           ELSE
              DO 12 K=1,L
                 A(I,K)=A(I,K)/SCALE
                 H=H+A(I,K)**2
12            END DO
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.
              DO 15 J=1,L
                 A(J,I)=A(I,J)/H
                 G=0.
                 DO 13 K=1,J
                    G=G+A(J,K)*A(I,K)
13               END DO
                 IF(L.GT.J)THEN
                    DO 14 K=J+1,L
                       G=G+A(K,J)*A(I,K)
14                  END DO
                 ENDIF
                 E(J)=G/H
                 F=F+E(J)*A(I,J)
15            END DO
              HH=F/(H+H)
              DO 17 J=1,L
                 F=A(I,J)
                 G=E(J)-HH*F
                 E(J)=G
                 DO 16 K=1,J
                    A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16               END DO
17            END DO
           ENDIF
        ELSE
           E(I)=A(I,L)
        ENDIF
        D(I)=H
18   END DO
  ENDIF
  D(1)=0.
  E(1)=0.
  DO 23 I=1,N
     L=I-1
     IF(D(I).NE.0.)THEN
        DO 21 J=1,L
           G=0.
           DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19         END DO
           DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20         END DO
21      END DO
     ENDIF
     D(I)=A(I,I)
     A(I,I)=1.
     IF(L.GE.1)THEN
        DO 22 J=1,L
           A(I,J)=0.
           A(J,I)=0.
22      END DO
     ENDIF
23 END DO

  RETURN
END SUBROUTINE TRED2
