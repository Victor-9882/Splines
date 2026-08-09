PROGRAM SPLINE_2D_ESCALAR
  IMPLICIT DOUBLE PRECISION (a-h, o-z)

  PARAMETER (NMG=10, NMZ=10, NMA=NMG*NMZ, LWORK=10*NMA)

  DOUBLE PRECISION :: XMATRIX(NMA,NMA), ZMATRIX(NMA,NMA)
  DOUBLE PRECISION :: VR(NMA,NMA), VL(1,1)
  DOUBLE PRECISION :: ALPHAR(NMA), ALPHAI(NMA), BETA(NMA)
  DOUBLE PRECISION :: WR(NMA), WI(NMA), WORK(LWORK)
  DOUBLE PRECISION :: zv(NMZ+1), gv(NMG+1)
  DOUBLE PRECISION :: c(NMG,NMZ), XG(NMA), YG(NMA)
  DOUBLE PRECISION :: Xq(1000), DXq(1000), Yp(1000), DYp(1000)
  DOUBLE PRECISION :: Wv(1000), DWv(1000)

  INTEGER :: ii, i, j, NPARAM, INFO, iw
  INTEGER :: N_intervalG, N_intervalZ, NCOL
  INTEGER :: Nnz(100), Nng(100), Nnv(100)
  INTEGER :: Nz, Ng, Nv
  DOUBLE PRECISION :: e, m, Mtot, mu, kappa, PI, gam0, Alfa
  DOUBLE PRECISION, EXTERNAL :: KERNEL_LOWER, KERNEL_UPPER

  open(unit=10, file="autovalores.dat",  STATUS="UNKNOWN")
  open(unit=11, file="autovetores.dat",  STATUS="UNKNOWN")
  open(unit=12, file="alfa.dat",         STATUS="UNKNOWN")
  open(unit=16, file="coeficientes.dat", STATUS="UNKNOWN")
  open(unit=20, file="inputs.dat",       STATUS="UNKNOWN")

  e  = 0.0001d0
  PI = DACOS(-1.d0)

  READ(20,*) NPARAM
  DO i = 1, NPARAM
    READ(20,*) Nnz(i), Nng(i), Nnv(i)
  END DO
  CLOSE(20)

  Mtot  = 1.d0
  m     = 1.d0
  mu    = 0.5d0
  kappa = SQRT(m**2 - 0.25d0*Mtot**2)
  gam0  = 10.d0

  WRITE(10,'(A,I0,A,I0,A,I0)') "NMA: ",NMA," NMG: ",NMG," NMZ: ",NMZ
  WRITE(10,'(A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10)') &
    "Mtot: ",Mtot," m: ",m," mu: ",mu," kappa: ",kappa," e: ",e," gam0: ",gam0
  WRITE(10,*) "--------------------------------------------------"

  WRITE(12,'(A,I0,A,I0,A,I0)') "NMA: ",NMA," NMG: ",NMG," NMZ: ",NMZ
  WRITE(12,'(A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10,A,F20.10)') &
    "Mtot: ",Mtot," m: ",m," mu: ",mu," kappa: ",kappa," e: ",e," gam0: ",gam0
  WRITE(12,*) "--------------------------------------------------"

  iw          = 14
  N_intervalZ = (NMZ-1)/2
  N_intervalG = (NMG-1)/2
  NCOL        = 2

  DO ii = 1, NPARAM
    print*, ii

    Nz = Nnz(ii)
    Ng = Nng(ii)
    Nv = Nnv(ii)

    ! z-mesh via collocação
    call G1D(iw, -1.d0, N_intervalZ, 1.d0, 1.d0, Xq)
    call COLLOC(iw, NCOL, N_intervalZ, Xq, XG)
    do i = 1, 2*N_intervalZ
      zv(i+1) = XG(i)
    end do
    zv(1)   = -0.999999999d0
    zv(NMZ) =  0.99999999d0

    ! gamma-mesh via collocação
    call G1D(iw, 0.d0, N_intervalG, 1.d0, 3.d0, Yp)
    call COLLOC(iw, NCOL, N_intervalG, Yp, YG)
    do i = 1, 2*N_intervalG
      gv(i+1) = YG(i)
    end do
    gv(1)   = 0.d0
    gv(NMG) = 3.d0

    call SPLGR1(zv, NMZ)
    call SPLGR2(gv, NMG)

    ! Nós e pesos de Gauss-Legendre
    CALL legauss(0.d0, 1.d0, Nz, Xq,  DXq,  1.d-15)
    CALL legauss(0.d0, 3.d0, Ng, Yp,  DYp,  1.d-15)
    CALL legauss(0.d0, 1.d0, Nv, Wv,  DWv,  1.d-15)

    CALL BUILD_ZMATRIX(zv, gv, NMZ, NMG, NMA, m, mu, kappa, PI, &
                       Xq, DXq, Yp, DYp, Wv, DWv, Nz, Ng, Nv, &
                       ZMATRIX, XMATRIX)

    do i = 1, NMA
      XMATRIX(i,i) = XMATRIX(i,i) + e
    end do

    CALL DGGEV('N', 'V', NMA, ZMATRIX, NMA, XMATRIX, NMA, &
               ALPHAR, ALPHAI, BETA, VL, 1, VR, NMA, &
               WORK, LWORK, INFO)

    IF (INFO /= 0) THEN
      PRINT *, 'ERRO NO DGGEV: INFO = ', INFO
      STOP
    END IF

    DO i = 1, NMA
      if (ABS(BETA(i)) > 1.d-16) then
        WR(i) = ALPHAR(i) / BETA(i)
        WI(i) = ALPHAI(i) / BETA(i)
      else
        WR(i) =  1.d+16
        WI(i) =  0.d0
      end if
    END DO

    WRITE(10,'(A,I0,A,I0,A,I0,A)') "autovalores_Nz",Nz,"_Ng",Ng,"_Nv",Nv,".dat"
    WRITE(12,'(A,I0,A,I0,A,I0,A)') &
      "Numericamente.Collocation.autovalores_Nz",Nz,"_Ng",Ng,"_Nv",Nv,".dat"

    do i = 1, NMA
      WRITE(10,'(I4,2X,F20.12,2X,F20.12)') i, WR(i), WI(i)
    end do

    WRITE(10,*) ""
    Alfa = 1.0d0 / (WR(1) * 16.0d0 * PI)
    WRITE(10,'(A,F20.10)') "Valor de Alfa: ", Alfa
    WRITE(12,'(A,F20.10)') "Valor de Alfa: ", Alfa

    WRITE(10,'(9999ES16.8)') (VR(j,1), j=1, NMA)

    do j = 1, NMZ
      do i = 1, NMG
        c(i,j) = DABS(VR(i + (j-1)*NMG, 1))
      end do
    end do

    DO i = 1, NMG
      WRITE(16,'(9999ES16.8)') (c(i,j), j=1, NMZ)
    END DO

  END DO

  CLOSE(10)
  CLOSE(11)
  CLOSE(12)
  CLOSE(16)

END PROGRAM SPLINE_2D_ESCALAR


! -----------------------------------------------------------------------
! BUILD_ZMATRIX
! Monta ZMATRIX (lado direito) e XMATRIX (lado esquerdo).
! Para trocar o kernel físico, modifique KERNEL_UPPER e KERNEL_LOWER.
! -----------------------------------------------------------------------
SUBROUTINE BUILD_ZMATRIX(zv, gv, Nmz, Nmg, NMA, m, mu, kappa, PI, &
                          Xq, DXq, Yp, DYp, Wv, DWv, Nz, Ng, Nv, &
                          ZMATRIX, XMATRIX)
  IMPLICIT DOUBLE PRECISION (a-h, o-z)

  INTEGER, INTENT(IN) :: Nmz, Nmg, NMA, Nz, Ng, Nv
  DOUBLE PRECISION, INTENT(IN)  :: zv(*), gv(*)
  DOUBLE PRECISION, INTENT(IN)  :: m, mu, kappa, PI
  DOUBLE PRECISION, INTENT(IN)  :: Xq(Nz), DXq(Nz)
  DOUBLE PRECISION, INTENT(IN)  :: Yp(Ng), DYp(Ng)
  DOUBLE PRECISION, INTENT(IN)  :: Wv(Nv), DWv(Nv)
  DOUBLE PRECISION, INTENT(OUT) :: ZMATRIX(NMA,NMA), XMATRIX(NMA,NMA)

  DOUBLE PRECISION :: splz_z(Nmz), splz_q(Nmz), splg(Nmg)
  DOUBLE PRECISION :: z, g, gp, zq, v, dzq, dgp, dv
  DOUBLE PRECISION, EXTERNAL :: KERNEL_UPPER, KERNEL_LOWER
  INTEGER :: i, j, k, l, p, q, r, index1, index2

  ZMATRIX = 0.d0
  XMATRIX = 0.d0

  do i = 1, Nmg
    g = gv(i)
    print*, i
    do j = 1, Nmz
      z      = zv(j)
      index1 = (j-1)*Nmg + i
      call SPLMD1(zv, Nmz, z, splz_z)

      do k = 1, Nmg
        do l = 1, Nmz
          index2 = (l-1)*Nmg + k

          do p = 1, Ng
            gp  = Yp(p)
            dgp = DYp(p)
            call SPLMD2(gv, Nmg, gp, splg)

            ! Lado esquerdo: integra apenas em gp
            XMATRIX(index1,index2) = XMATRIX(index1,index2) + &
              splg(k)*splz_z(l)*dgp / &
              (g + gp + (1.d0 - z**2)*kappa**2 + m**2*z**2)**2

            do q = 1, Nz

              ! Região superior: z' em [z, 1]
              dzq = (1.d0 - z)*DXq(q)
              zq  = (1.d0 - z)*Xq(q) + z
              call SPLMD1(zv, Nmz, zq, splz_q)
              do r = 1, Nv
                v  = Wv(r)
                dv = DWv(r)
                ZMATRIX(index1,index2) = ZMATRIX(index1,index2) + &
                  KERNEL_UPPER(z, zq, g, gp, v, m, mu, kappa, PI) * &
                  splg(k)*splz_q(l)*dzq*dgp*dv
              end do

              ! Região inferior: z' em [-1, z]
              dzq = (z + 1.d0)*DXq(q)
              zq  = (z + 1.d0)*Xq(q) - 1.d0
              call SPLMD1(zv, Nmz, zq, splz_q)
              do r = 1, Nv
                v  = Wv(r)
                dv = DWv(r)
                ZMATRIX(index1,index2) = ZMATRIX(index1,index2) + &
                  KERNEL_LOWER(z, zq, g, gp, v, m, mu, kappa, PI) * &
                  splg(k)*splz_q(l)*dzq*dgp*dv
              end do

            end do  ! q
          end do  ! p
        end do  ! l
      end do  ! k
    end do  ! j
  end do  ! i

END SUBROUTINE BUILD_ZMATRIX


! -----------------------------------------------------------------------
! KERNEL_UPPER  —  integrando do kernel para z' em [z, 1]
! -----------------------------------------------------------------------
DOUBLE PRECISION FUNCTION KERNEL_UPPER(z, zp, g, gp, v, m, mu, kappa, PI)
  IMPLICIT DOUBLE PRECISION (a-h, o-z)
  DOUBLE PRECISION, INTENT(IN) :: m, mu, kappa

  D = v*(1.d0-v)*(1.d0+zp)*g &
    + v*(1.d0+z)*gp &
    + v*(1.d0+z)*(1.d0+zp)*(1.d0 - z*(1.d0-v) - v*zp)*kappa**2 &
    + v*((1.d0-v)*(1.d0+zp)*z**2 + v*zp**2*(1.d0+z))*m**2 &
    + (1.d0-v)*(1.d0+z)*mu**2

  KERNEL_UPPER = (1.d0+z)**2 &
    / (32.d0*PI**2 * (g + z**2*m**2 + (1.d0-z**2)*kappa**2)) &
    * v**2 / D**2

END FUNCTION KERNEL_UPPER


! -----------------------------------------------------------------------
! KERNEL_LOWER  —  integrando do kernel para z' em [-1, z]
! -----------------------------------------------------------------------
DOUBLE PRECISION FUNCTION KERNEL_LOWER(z, zp, g, gp, v, m, mu, kappa, PI)
  IMPLICIT DOUBLE PRECISION (a-h, o-z)
  DOUBLE PRECISION, INTENT(IN) :: m, mu, kappa

  D = v*(1.d0-v)*(1.d0-zp)*g &
    + v*(1.d0-z)*gp &
    + v*(1.d0-z)*(1.d0-zp)*(1.d0 + z*(1.d0-v) + v*zp)*kappa**2 &
    + v*((1.d0-v)*(1.d0-zp)*z**2 + v*zp**2*(1.d0-z))*m**2 &
    + (1.d0-v)*(1.d0-z)*mu**2

  KERNEL_LOWER = (1.d0-z)**2 &
    / (32.d0*PI**2 * (g + z**2*m**2 + (1.d0-z**2)*kappa**2)) &
    * v**2 / D**2

END FUNCTION KERNEL_LOWER


! -----------------------------------------------------------------------
! Rotinas de Spline
! -----------------------------------------------------------------------
      SUBROUTINE SPLGR1 (X,N)
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER (NP1=500)
      DIMENSION X(N+1),HI(NP1),U(NP1),Q(NP1,NP1),C(NP1,NP1)
      COMMON /FActz/ FAK1(NP1,NP1),FAK2(NP1,NP1),FAK3(NP1,NP1)
      U(1)=0.D0
      HI(2)=X(2)-X(1)
      DO 5 I=1,N
    5  Q(1,I)=0.D0
      DO 10 I=2,N-1
       AX=X(I+1)-X(I)
       HI(I+1)=AX
       BX=X(I+1)-X(I-1)
       CX=X(I)-X(I-1)
       AL=AX/BX
       AM=1.D0-AL
       PI=1.D0/(2.D0-AM*U(I-1))
       U(I)=AL*PI
       DO 15 J=1,N
   15   Q(I,J)=-PI*AM*Q(I-1,J)
       Q(I,I-1)=Q(I,I-1)+PI/(CX*BX)
       Q(I,I)=Q(I,I)-PI/(CX*AX)
   10  Q(I,I+1)=Q(I,I+1)+PI/(AX*BX)
      DO 20 J=1,N
       C(N,J)=0.D0
       FAK1(N,J)=0.D0
       FAK2(N,J)=0.D0
   20  FAK3(N,J)=0.D0
      DO 25 I=N-1,1,-1
       H1=1.D0/HI(I+1)
       DO 30 J=1,N
        C(I,J)=Q(I,J)-C(I+1,J)*U(I)
   30   FAK1(I,J)=-HI(I+1)*(2.D0*C(I,J)+C(I+1,J))
       FAK1(I,I)=FAK1(I,I)-H1
       FAK1(I,I+1)=FAK1(I,I+1)+H1
       DO 25 J=1,N
        FAK2(I,J)=3*C(I,J)
   25   FAK3(I,J)=(C(I+1,J)-C(I,J))*H1
      return
      END

      SUBROUTINE SPLMD1 (X,N,XA,SPL)
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER (NP1=500)
      DIMENSION X(N+1),SPL(N)
      COMMON /FActz/ FAK1(NP1,NP1),FAK2(NP1,NP1),FAK3(NP1,NP1)
      I=-1
   10  I=I+1
       IF (XA .GE. X(I+1) .AND. I .LT. N) GOTO 10
      IF (I .EQ. 0) I=1
      DX=XA-X(I)
      DO 20 J=1,N
   20  SPL(J)=((FAK3(I,J)*DX+FAK2(I,J))*DX+FAK1(I,J))*DX
      SPL(I)=SPL(I)+1.D0
      return
      END

      SUBROUTINE SPLGR2 (X,N)
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER (NP1=500)
      DIMENSION X(N+1),HI(NP1),U(NP1),Q(NP1,NP1),C(NP1,NP1)
      COMMON /FActg/ FAK4(NP1,NP1),FAK5(NP1,NP1),FAK6(NP1,NP1)
      U(1)=0.D0
      HI(2)=X(2)-X(1)
      DO 5 I=1,N
    5  Q(1,I)=0.D0
      DO 10 I=2,N-1
       AX=X(I+1)-X(I)
       HI(I+1)=AX
       BX=X(I+1)-X(I-1)
       CX=X(I)-X(I-1)
       AL=AX/BX
       AM=1.D0-AL
       PI=1.D0/(2.D0-AM*U(I-1))
       U(I)=AL*PI
       DO 15 J=1,N
   15   Q(I,J)=-PI*AM*Q(I-1,J)
       Q(I,I-1)=Q(I,I-1)+PI/(CX*BX)
       Q(I,I)=Q(I,I)-PI/(CX*AX)
   10  Q(I,I+1)=Q(I,I+1)+PI/(AX*BX)
      DO 20 J=1,N
       C(N,J)=0.D0
       FAK4(N,J)=0.D0
       FAK5(N,J)=0.D0
   20  FAK6(N,J)=0.D0
      DO 25 I=N-1,1,-1
       H1=1.D0/HI(I+1)
       DO 30 J=1,N
        C(I,J)=Q(I,J)-C(I+1,J)*U(I)
   30   FAK4(I,J)=-HI(I+1)*(2.D0*C(I,J)+C(I+1,J))
       FAK4(I,I)=FAK4(I,I)-H1
       FAK4(I,I+1)=FAK4(I,I+1)+H1
       DO 25 J=1,N
        FAK5(I,J)=3*C(I,J)
   25   FAK6(I,J)=(C(I+1,J)-C(I,J))*H1
      return
      END

      SUBROUTINE SPLMD2 (X,N,XA,SPL)
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER (NP1=500)
      DIMENSION X(N+1),SPL(N)
      COMMON /FActg/ FAK4(NP1,NP1),FAK5(NP1,NP1),FAK6(NP1,NP1)
      I=-1
   10  I=I+1
       IF (XA .GE. X(I+1) .AND. I .LT. N) GOTO 10
      IF (I .EQ. 0) I=1
      DX=XA-X(I)
      DO 20 J=1,N
   20  SPL(J)=((FAK6(I,J)*DX+FAK5(I,J))*DX+FAK4(I,J))*DX
      SPL(I)=SPL(I)+1.D0
      return
      END


! -----------------------------------------------------------------------
! Quadratura de Gauss-Legendre
! -----------------------------------------------------------------------
      SUBROUTINE legauss(XS,XL,N,X,DX,ZZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension X(N),DX(N)
      IF(N)10,10,20
 10   WRITE(5,600) N
 600  FORMAT(1H ,I10,' REJEITADOS PTOS.LEG-GAUSS')
      RETURN
 20   IF(N-2) 30,40,40
 30   X(1)=0.D0
      DX(1)=.5D0
      GO TO 140
 40   I=1
      G=-1.D0
      IC=(N+1)/2
 50   S=G
      T=1.D0
      U=1.D0
      V=0.D0
      DO 60 K=2,N
      A=K
      FACT1=(2.D0*A-1.D0)/A
      FACT2=(A-1.D0)/A
      P=FACT1*G*S-FACT2*T
      DP=FACT1*(S+G*U)-FACT2*V
      T=S
      S=P
      V=U
 60   U=DP
      SUM=0.D0
      IF(I-1)90,90,70
 70   IM1=I-1
      DO 80 K=1,IM1
 80   SUM=SUM+1.D0/(G-X(K))
 90   TEST=G
      G=G-P/(DP-P*SUM)
      R=DABS(TEST-G)
      IF(R.LT.ZZ)GOTO 100
      GOTO 50
 100  R=N
      X(I)=G
      DX(I)=2.D0/R/T/DP
      IF(IC-I)120,120,110
 110  FIM1=IM1
      G=G-(DP-P*SUM)/((2.D0*G*DP-A*(A+1.D0)*P)/(1.D0-G*G)-2.D0*DP*SUM-P*SUM**2+FIM1*P)
      I=I+1
      GOTO 50
 120  K0=2*IC-N+2*(N/2)+1
      IC=IC+1
      DO 130 I=IC,N
      K=K0-I
      X(I)=-X(K)
 130  DX(I)=DX(K)
 140  FACT1=(XL-XS)/2.D0
      FACT2=(XL+XS)/2.D0
      DO 150 I=1,N
      DX(I)=DX(I)*FACT1
 150  X(I)=X(I)*FACT1+FACT2
      RETURN
      END


! -----------------------------------------------------------------------
! Geração de malha 1D e pontos de colocação
! -----------------------------------------------------------------------
      SUBROUTINE G1D(IW,X0,N,A,XN,X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(0:N)
      X(0)=X0
      DX=(XN-X0)/DFLOAT(N)
      IF(A.NE.1.D0)DX=(XN-X0)*(A-1.D0)/(A**N-1.D0)
      DO I=1,N
      X(I)=X(I-1)+DX
      DX=DX*A
      enddo
      X(N)=XN
      HMIN=DMIN1(X(1)-X(0),X(N)-X(N-1))
      HMAX=DMAX1(X(1)-X(0),X(N)-X(N-1))
      IF(IW.EQ.0)RETURN
      WRITE(IW,100)
  100 FORMAT(/,2X,'ONE-DOMAIN GRID (G1D) CHARACTERISED BY'&
     //,2X,4X,'Xmin',7X,'Xmax ',6X,'N ',6X,'A',6X,'Hmin',4X,'Hmax',/)
      WRITE(IW,101)X0,XN,N,A,HMIN,HMAX
  101 FORMAT(2X,D10.5,2X,D10.5,2X,I4,4X,F6.4,2X,F10.4,2X,F10.4)
      WRITE(IW,102) (X(I),I=0,N)
  102 FORMAT(/,6(2X,D15.8))
      RETURN
      END

      SUBROUTINE COLLOC(IW,NCOL,N,X,XG)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(0:N),XG(NCOL*N)
      IF(NCOL.EQ.2)THEN
      IG=1
      DO 2 I=1,N
      A=X(I-1)
      B=X(I)
      BPA=B+A
      BMA=B-A
      U=-0.577350269189626D0
      XG(IG)=0.5D0*(BMA*U+BPA)
      IG=IG+1
      U=+0.577350269189626D0
      XG(IG)=0.5D0*(BMA*U+BPA)
      IG=IG+1
    2 CONTINUE
      ENDIF
      IF(NCOL.EQ.3)THEN
      IG=1
      DO 3 I=1,N
      A=X(I-1)
      B=X(I)
      BPA=B+A
      BMA=B-A
      U=-0.774596669241483D0
      XG(IG)=0.5D0*(BMA*U+BPA)
      IG=IG+1
      U=+0.0D0
      XG(IG)=0.5D0*(BMA*U+BPA)
      IG=IG+1
      U=+0.774596669241483D0
      XG(IG)=0.5D0*(BMA*U+BPA)
      IG=IG+1
    3 CONTINUE
      ENDIF
      IF(IW.EQ.0)RETURN
      WRITE(IW,100)
  100 FORMAT(/,2X,'COLLOCATION GRID',/)
      WRITE(IW,102) (XG(I),I=1,NCOL*N)
  102 FORMAT(6(2X,D15.8))
      RETURN
      END
