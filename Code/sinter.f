      SUBROUTINE SINTER(X,A,Y,B,M,N)
C     VERSION(03/02/90)
C
C----------------------------------------------------------------------|
C       A SPECIAL CASE OF INTERP ....NO EXTRAPOLATION BELOW DATA       |
C       THIS ROUTINE LINEARLY INTERPOLATES AND EXTRAPOLATES AN ARRAY B |
C       X(M) MUST BE DESCENDING                                        |
C       A(X) GIVEN FUNCTION                                            |
C       B(Y) FOUND BY LINEAR INTERPOLATION AND EXTRAPOLATION           |
C       Y(N) THE DESIRED DEPTHS                                        |
C       M    THE NUMBER OF POINTS IN X AND A                           |
C       N    THE NUMBER OF POINTS IN Y AND B                           |
C----------------------------------------------------------------------|
C
      DIMENSION X(M),A(M),Y(N),B(N)
C
C-------- EXTRAPOLATION CASES ------------------------------------------
      DO 30 I=1,N
      IF(Y(I).GT.X(1)) B(I)=A(1)+((A(1)-A(2))/(X(1)-X(2)))*(Y(I)-X(1))
      IF(Y(I).LT.X(M)) B(I)=A(M)
 30   CONTINUE
C
C-------- INTERPOLATION CASES ------------------------------------------
      NM=M-1
      DO 10 I=1,N
      DO 20 J=1,NM
      IF(Y(I).LE.X(J).AND.Y(I).GE.X(J+1))
     +  B(I)=A(J) -(A(J)-A(J+1))*(X(J)-Y(I))/(X(J)-X(J+1))
 20   CONTINUE
 10   CONTINUE
C
      RETURN
      END
