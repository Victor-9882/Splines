IMPLICIT DOUBLE PRECISION(a-h,o-z)

    INTEGER :: NMG, NMZ, NMA, LWORK
    DOUBLE PRECISION, ALLOCATABLE :: XMATRIX(:,:), ZMATRIX(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: VR(:,:), VL(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: WR(:), WI(:), WORK(:)
    DOUBLE PRECISION, ALLOCATABLE :: zv(:), splz(:), splg(:), gv(:)
    DOUBLE PRECISION, ALLOCATABLE :: dzv(:), dgv(:), ALPHAR(:), ALPHAI(:)
    DOUBLE PRECISION, ALLOCATABLE :: BETA(:), c(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: XG(:), YG(:)

    INTEGER, ALLOCATABLE :: ipvt(:)
    COMMON/PARAM/PI
    COMMON/ALPHAINT/X(1000),DX(1000), Y(1000),DY(1000) &
    , W(1000),DW(1000), Nz, Ng, Nv, det (1000)
    INTEGER :: ii, jj, i, j, k, l & 
    ,index1, index2, p, q, r, NPARAM, N_intervalZ, NCOL, N_intervalG
    INTEGER :: INFO, LDVL, IPRINTEIGEN, N_PLOT, k_max, k_min
    DOUBLE PRECISION :: LAMBDAR,  LAMBDAI, e, m, Mtot, mu, kappa, gam0, soma &
        , IntegralV, num2, max_gamma_visualizacao

    INTEGER :: Nnz (100), Nng (100), Nnv (100) 

    DOUBLE PRECISION, EXTERNAL :: f_map, Jacobian_map, inverse_map


    open (unit = 10, file = "autovalores.dat",STATUS="UNKNOWN")
    open (unit = 11, file = "autovetoresG.dat",STATUS="UNKNOWN")
    open (unit = 15, file = "autovetoresZ.dat",STATUS="UNKNOWN")
    open (unit = 16, file = "coeficientes.dat",STATUS="UNKNOWN")
    open (unit = 12, file = "alfa.dat",STATUS="UNKNOWN")
    open (unit = 13, file = "plotalfa.dat",STATUS="UNKNOWN")
    open (unit = 14, file = 'erros.dat', status='unknown')
    open (UNIT = 20, FILE = "inputs.dat", STATUS="UNKNOWN")

    e = 0.0001d0
    PI = DACOS(-1.D0)       !3.14159265358979323846264338

        !Parâmetros
        !Massas
        Mtot = 1.0d0
        m = 1.0d0
        mu = 0.50d0
        kappa = sqrt(m**2 - 0.25*Mtot**2)

        gam0 = 10.0d0

        Nz = 60
        Ng = 60
        Nv = 60

        READ(20,*) NPARAM
        DO I = 1, NPARAM
            READ(20,*) Nng(i), Nnz(i)
        END DO
        CLOSE(20)


    do ii = 1, NPARAM

        NMG = Nng(ii)
        NMZ = Nnz(ii)

        NMA = NMG * NMZ
        LWORK = 10 * NMA

        !ALOCAR VARIÁVEIS
        ALLOCATE( XMATRIX(NMA, NMA), ZMATRIX(NMA, NMA))
        ALLOCATE( VR(NMA, NMA), VL(NMA, 2*NMA) )
        ALLOCATE( WR(NMA), WI(NMA), WORK(LWORK) )
        ALLOCATE( zv(NMZ + 1), gv(NMG+1) )
        ALLOCATE( splz(NMZ), splg(NMG))
        ALLOCATE( ALPHAR(NMA), ALPHAI(NMA), BETA(NMA) )
        ALLOCATE( c(NMA, NMA) )
        ALLOCATE( XG(NMA), YG(NMA) )
        
        ALLOCATE( ipvt(NMA) )

    ! Ler matriz de coeficientes
      DO I = 1, NMG
         READ(16, *) (c(I,J), J=1, NMZ)
      END DO

      ! 3. Feche o arquivo
      CLOSE(16)

        iw = 14
        N_intervalZ = (NMZ-1)/2
        N_intervalG = (NMG-1)/2
        NCOL = 2    

        !Contrução das malhas
            
            call G1D(IW,-1.d0, N_intervalZ, 1.d0, 1.d0, X)
            call COLLOC(IW,2,N_intervalZ,X,XG)  
            
            do i = 1, 2*N_intervalZ
            zv(i+1) = XG(i)
            end do

            zv(1)=-0.999999d0
            zv(nmz)= 0.999999d0


            call G1D(IW,0.d0, N_intervalG, 1.0d0, 3.d0, Y)
            call COLLOC(IW,2,N_intervalG,Y,YG)  

            do i=1, 2*N_intervalG
            gv(i+1)=YG(i)
            enddo

            gv(1) = 0.0000001d0
            gv(nmg) = 3.d0
        

            call SPLGR1 (zv,Nmz)
            call SPLGR2 (gv,Nmg)


      !Plot de função x gamma para z fixo
      z_fixo = 0.4d0
      max_gamma_visualizacao = 3.d0
      N_PLOT = 1000
      call SPLMD1(zv, Nmz, z_fixo, splz)

          do p = 0, N_PLOT
              ! Cria uma distribuição de pontos
              gamma_plot = (dble(p) / dble(N_PLOT)) * max_gamma_visualizacao
              !t_plot = inverse_map(gamma_plot)
              
              ! Avalia as bases de Spline no ponto gamma_plot atual
              call SPLMD2(gv, Nmg, gamma_plot, splg)
              !call SPLMD2(gv, Nmg, t_plot, splg)

              soma = 0.d0
              do j = 1, Nmz  
                  ! Como z_fixo é constante, splz(j) já foi calculado fora do loop de p
                  do i = 1, Nmg
                      soma = soma + c(i,j) * splg(i) * splz(j)
                  end do
              end do
              if (soma < 0.0d0) then
                    ! Adiciona um espaço fixo ' ' antes de imprimir os números
                    write(11, '(ES25.17E3, 1X, ES25.17E3)') gamma_plot, soma
                else
                    ! Imprime normalmente (o descritor ES já deixa um espaço natural para positivos)
                    write(11, '(2ES25.17E3)') gamma_plot, soma
                end if
          end do

      ! Plot da função variando z para gamma fixo
      gamma_fixo = 1.d0
      N_PLOT = 1500
      call SPLMD2(gv, Nmg, gamma_fixo, splg)

      do p = 0, N_PLOT
          ! Cria uma distribuição linear para z_plot no domínio de -1.0 a 1.0
          z_plot = -cos( (dble(p) / dble(N_PLOT)) * PI )
          
          ! Avalia as bases de Spline no ponto z_plot atual
          call SPLMD1(zv, Nmz, z_plot, splz)

          soma = 0.d0
          do j = 1, Nmz  
              do i = 1, Nmg
                  ! splg(i) é constante nesta etapa, splz(j) está variando
                  soma = soma + c(i,j) * splg(i) * splz(j)
              end do
          end do
          
          if (soma < 0.0d0) then
              ! Adiciona um espaço fixo ' ' antes de imprimir os números
              write(15, '(ES25.17E3, 1X, ES25.17E3)') z_plot, soma
          else
              ! Imprime normalmente
              write(15, '(2ES25.17E3)') z_plot, soma
          end if
      end do

      ! -------------------------------------------------------------
      ! PLOT DA FUNÇÃO PSI NORMALIZADA: Psi(gamma, 0.2) / Psi(0, 0.5)
      ! -------------------------------------------------------------
      open(unit = 16, file = "plot_psi.dat", STATUS="UNKNOWN")

      ! 1. CÁLCULO DA NORMALIZAÇÃO (Denominador): Psi(0.0, 0.5)
      gamma_den = 0.0d0
      z_den     = 0.5d0
      psi_den   = 0.d0
      
      !Xi = (1-z)/2
      z_den = 1 - 2*z_den
      ! Avalia a spline em z_den (Fixo fora do loop da integral)
      call SPLMD1(zv, Nmz, z_den, splz)
      CALL legauss(0.d0,3.d0,Ng,Y,dY,1.d-15)
      ! Integrando em gama' usando os pontos de Gauss Y(p) em [-1, 1]
      do p = 1, Ng
          tp = Y(p)         
          wtp = dY(p)       
          
          gammap = tp
          dgp = wtp
          
          ! Avalia a spline em gama' (auxiliar tp)
          call SPLMD2(gv, Nmg, tp, splg)
          
          ! Constrói g(gammap, z_den)
          g_val = 0.d0
          do j = 1, Nmz
              do i = 1, Nmg
                  g_val = g_val + c(i,j) * splg(i) * splz(j)
              end do
          end do
          
          ! Denominador do propagador: [gamma + gamma' + m^2*z^2 + (1-z^2)*kappa^2]
          D_den = gamma_den + gammap + (m**2)*(z_den**2) + (1.d0 - z_den**2)*(kappa**2)
          
          ! Acumula a integral de dgamma' * g / [Propagador]^2
          psi_den = psi_den + (g_val * dgp) / (D_den**2)
      end do


      ! 2. CÁLCULO DO NUMERADOR E PLOT: Psi(gamma_plot, z0)
      z_num = 0.2d0
      max_gamma_visualizacao = 3.d0
      N_PLOT = 1000
      
      z_num = 1 - 2*z_num
      ! Avalia a spline no novo z_num (Fixo para todo o plot)
      call SPLMD1(zv, Nmz, z_num, splz)
      
      do k_plot = 0, N_PLOT
          gamma_plot = (dble(k_plot) / dble(N_PLOT)) * max_gamma_visualizacao
          psi_num = 0.d0
          
          ! Integrando em gama' para o gamma_plot atual
          do p = 1, Ng
              tp = Y(p)
              wtp = dY(p)
              
              gammap = tp
              dgp = tp * wtp
              
              call SPLMD2(gv, Nmg, tp, splg)
              
              g_val = 0.d0
              do j = 1, Nmz
                  do i = 1, Nmg
                      g_val = g_val + c(i,j) * splg(i) * splz(j)
                  end do
              end do
              
              D_num = gamma_plot + gammap + (m**2)*(z_num**2) + (1.d0 - z_num**2)*(kappa**2)
              
              psi_num = psi_num + (g_val * dgp) / (D_num**2)
          end do
          
          ! Faz a normalização final cancelando as constantes
          psi_norm = psi_num / psi_den
          
          ! Salva no arquivo (gamma_plot no Eixo X, psi_norm no Eixo Y)
          write(16, '(2ES25.17E3)') gamma_plot, psi_norm
      end do
      
      close(16)

      ! -------------------------------------------------------------
      ! PLOT DA FUNÇÃO PSI NORMALIZADA: Psi(gamma0, z) / Psi(0, 0.5)
      ! -------------------------------------------------------------
      open(unit = 17, file = "plot_psi_z.dat", STATUS="UNKNOWN")

      ! Numerador tem gamma fixo em gamma0
      gamma_num = 0.3d0
      N_PLOT = 1500
      
      do k_plot = 0, N_PLOT
          ! Variando z_plot entre -1 e 1
          ! Usando a distribuição cosseno (que aglomera pontos nas bordas, 
          ! ideal para capturar bem o comportamento das splines)
          z_plot = -cos( (dble(k_plot) / dble(N_PLOT)) * PI )
          
          psi_num = 0.d0
          
          ! Avalia a spline no z_plot atual (que muda a cada passo do loop)
          call SPLMD1(zv, Nmz, z_plot, splz)
          
          ! Integrando em gama' usando os pontos de Gauss Y(p) em [0, inf]
          do p = 1, Ng
              tp = Y(p)
              wtp = dY(p)
              
              gammap = tp
              dgp = wtp
              
              ! Avalia a spline em gamma'
              call SPLMD2(gv, Nmg, tp, splg)
              
              ! Constrói g(gammap, z_plot) combinando as splines e os coeficientes
              g_val = 0.d0
              do j = 1, Nmz
                  do i = 1, Nmg
                      g_val = g_val + c(i,j) * splg(i) * splz(j)
                  end do
              end do
              
              ! Denominador do propagador usando gamma_num fixo e z_plot variável
              D_num = gamma_num + gammap + (m**2)*(z_plot**2) + (1.d0 - z_plot**2)*(kappa**2)
              
              ! Acumula a integral
              psi_num = psi_num + (g_val * dgp) / (D_num**2)
          end do
          
          ! Normalização final
          psi_norm = psi_num / psi_den
          
          ! Salva no arquivo (z_plot no Eixo X, psi_norm no Eixo Y)
          write(17, '(2ES25.17E3)') z_plot, psi_norm
      end do
      
      close(17)


      DEALLOCATE(XMATRIX, ZMATRIX, WORK, VR, VL, WR, WI, splz)
      DEALLOCATE(splg, ALPHAR, ALPHAI, BETA, c, XG, YG, gv, zv)
      DEALLOCATE(IPVT)
      

      end do

      close (15)      
      CLOSE(10)
      CLOSE(12)
      close(13)
10     FORMAT(11E12.4)
18     format(5e15.6)
20     FORMAT(A70)

       
      close (14)
       Close(2)
    END



      !Rotinas
        !Spline 1
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
    
        !Spline 2
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
    
    
    
    
    
    
    
    
    

      SUBROUTINE legauss(XS,XL,N,X,DX,ZZ)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension X(N),DX(N)
      IF(N)10,10,20
 10   WRITE(5,600) N
      WRITE(2,600) N
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


        !=====================================================================          
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
!=====================================================================          
      SUBROUTINE COLLOC(IW,NCOL,N,X,XG)                                         
!                                                                               
!     RETURN ABCISSES OF NCOL=2,3 GAUSS COLLOCATION POINTS              
!     ON EACH OF THE N INTERVALS OF THE GRID X(0),X(1),...,X(N)                          
!                                                                               
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











        