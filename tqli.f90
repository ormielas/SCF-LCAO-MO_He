SUBROUTINE TQLI(D,E,N,NP,Z)

  IMPLICIT NONE
  
  INTEGER                             :: N, NP
  DOUBLE PRECISION, DIMENSION(NP, NP) :: Z
  DOUBLE PRECISION, DIMENSION(NP)     :: D, E
  
  INTEGER          :: I, ITER, K, L, M
  DOUBLE PRECISION :: B, C, DD, F, G, P, R, S
  
  IF ( N .GT. 1 ) THEN
     DO 11 I = 2, N
        E(I-1) = E(I)
11   END DO
     E(N) = 0.0D0
     DO 15 L = 1, N
        ITER = 0
1       DO 12 M = L, N-1
           DD = ABS(D(M))+ABS(D(M+1))
           IF ( ABS(E(M))+DD .EQ. DD ) GO TO 2
12      END DO
        M=N
2       IF ( M .NE. L) THEN
           IF ( ITER .EQ. 30 ) EXIT
           !IF(ITER.EQ.30) EXIT 'too many iterations'
           ITER = ITER + 1
           G = (D(L+1)-D(L))/(2.0D0*E(L))
           R = SQRT(G**2+1.0D0)
           G = D(M)-D(L)+E(L)/(G+SIGN(R,G))
           S = 1.0D0
           C = 1.0D0
           P = 0.0D0
           DO 14 I = M-1, L, -1
              F = S*E(I)
              B = C*E(I)
              IF ( ABS(F) .GE. ABS(G) ) THEN
                 C = G/F
                 R = SQRT(C**2 + 1.0D0)
                 E(I+1) = F*R
                 S = 1.0D0/R
                 C = C*S
              ELSE
                 S = F/G
                 R = SQRT(S**2 + 1.0D0)
                 E(I+1) = G*R
                 C = 1.0D0/R  
                 S = S*C
              ENDIF
              G     = D(I+1) - P
              R     = (D(I)-G)*S + 2.0D0*C*B
              P     = S*R
              D(I+1)= G + P
              G     = C*R - B
              DO 13 K = 1, N
                 F        = Z(K,I+1)
                 Z(K,I+1) = S*Z(K,I) + C*F
                 Z(K,I  ) = C*Z(K,I) - S*F
13            END DO
14         END DO
           D(L)=D(L)-P
           E(L)=G
           E(M)=0.0D0
           GO TO 1
        ENDIF
15   END DO
  ENDIF

  RETURN
END SUBROUTINE TQLI
