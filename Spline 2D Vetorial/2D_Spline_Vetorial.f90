 IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(NMG=10, NMZ=10, NMA = NMG*NMZ,LWORK=10*NMA)   
	  COMPLEX*16 IMAG     
      	DIMENSION XMATRIX(NMA,NMA), ZMATRIX(NMA,NMA),AUNIT(NMA,NMA) &
         ,VR(NMA,NMA),VL(NMA,2*NMA), BETAA(NMA),VAUX(NMA),WR(NMA)&
         ,WI(NMA),WORK(LWORK),zv(nmz+1), gv(300), splz(nmz), splg(nmg) &
         ,splz_at(nmz), dzv(nmz),dgv(nmg), ALPHAR (NMA), ALPHAI(NMA) &
         ,BETA (NMA), a(nma, nma), b(nma, nma), c(nma, nma), T(300), DT(300) &
         , XG(NMA), YG(NMA)
       COMMON/PARAM/PI
       COMMON/ALPHAINT/X(1000),DX(1000), Y(1000),DY(1000) &
       , W(1000),DW(1000), Nz, Ng, Nv, det (1000)
       INTEGER :: ii, jj, i, j, k, l, iw & 
       ,index1, index2, p, q, r, NPARAM
        INTEGER ipvt(NMA), N_intervalG, N_intervalZ, NCOL
       INTEGER :: INFO, LDVL, IPRINTEIGEN, N_PLOT
       double precision LAMBDAR,  LAMBDAI, e, m, Mtot, mu, kappa, gam0, soma &
            , IntegralV, num2, max_gamma_visualizacao, Numerador
       CHARACTER(LEN=200) :: filename
       INTEGER :: Nnz (100), Nng (100), Nnv (100)

        open (unit = 10, file = "autovalores.dat",STATUS="UNKNOWN")
        open (unit = 12, file = "alfa.dat",STATUS="UNKNOWN")
        open (unit = 11, file = "autovetores.dat",STATUS="UNKNOWN")
        open (UNIT = 20, FILE = "inputs.dat", STATUS="UNKNOWN")
    
		    IMAG=DCMPLX(0.D0,1.D0)
        e = 0.000001d0
        PI = DACOS(-1.D0)       !3.14159265358979323846264338

      READ(20,*) NPARAM
          DO I = 1, NPARAM
           READ(20,*) nnz(I), nng(I), nnv(I)
          END DO
    CLOSE(20)
        
     !Parâmetros
          !Massas
          Mtot = 1.99d0
          m = 1.0d0
          m1 = 1.0d0
          m2 = 1.0d0
          mbar = (m1+m2)/2
          delta = (m2-m1)/2
          mu = 0.15d0
          kappa = sqrt(m**2 - 0.25*Mtot**2)

          gam0 = 10.0d0
        
          WRITE(10, '(A, I0, A, I0, A, I0)') "NMA: ", nma, " NMG: ", nmg, " NMZ: ", nmz
          WRITE(10, '(A, F20.10, A, F20.10, A, F20.10, A, F20.10, A, F20.10, A, F20.10)') &
         & "Mtot: ", Mtot, " mu: ", mu, " kappa: ", kappa, " e", e, " gam0", gam0, &
          "m: ", m
          WRITE(10, *) "--------------------------------------------------"

          WRITE(12, '(A, I0, A, I0, A, I0)') "NMA: ", nma, " NMG: ", nmg, " NMZ: ", nmz
          WRITE(12, '(A, F20.10, A, F20.10, A, F20.10, A, F20.10, A, F20.10, A, F20.10)') &
         & "Mtot: ", Mtot, " mu: ", mu, " kappa: ", kappa, " e", e, " gam0", gam0, &
          "m: ", m
          WRITE(12, *) "--------------------------------------------------"
          !print*, Mtot, m, mu, kappa


            iw = 14
            N_intervalZ = (NMZ-1)/2
            N_intervalG = (NMG-1)/2
            NCOL = 2
        do ii = 1, NPARAM
            print*, ii
      
          !Número de pontos de Gauss para integração em cada variável
          Nz = Nnz(ii)
          Ng = Nng (ii)
          Nv = Nnv (ii)
          
        
      !Contrução das malhas
            
         !Malha do z
          !CALL legauss(-1.d0,0.97d0,Nmz/2,zv,dzv,1.d-15)
        	!CALL legauss(0.97d0,1.d0,Nmz/2,X,dX,1.d-15)

            !CALL legauss(-1.d0,1d0,Nmz,zv,dzv,1.d-15)
            
        !do i=nmz/2+1,nmz
           !zv(i)=x(i-nmz/2)
        !end do
    
        !zv(1)=-1.d0
        !zv(nmz)=1.d0


        call G1D(IW,-1.d0, N_intervalZ, 1.2d0, 1.d0, X)
        call COLLOC(IW,2,N_intervalZ,X,XG)  
        
        do i = 1, 2*N_intervalZ
          zv(i+1) = XG(i)
        end do

        zv(1)=-0.999999999d0
        zv(nmz)=0.99999999d0
        
        !Malha da Gamma
        !CALL legauss(0.d0,3.d0,Nmg,gv,dY,1.d-15)
           

        call G1D(IW,0.d0, N_intervalG, 1.0d0, 3.d0, Y)
        call COLLOC(IW,2,N_intervalG,Y,YG)  

        do i=1, 2*N_intervalG
           gv(i+1)=YG(i)
        enddo

        gv(1) = 0.d0
        gv(nmg) = 3.d0


        !call setgaulag(0.d0,Nmg,gv,dY)

        
        !do i=1,nmg
          !gv(i)=2.d0/gam0*gv(i)
        !end do

        !do i=1, nmg
           !gv(i)=Y(i)
        !enddo

            !CALL legauss(0.d0, 1.d0, Nmg, T, dT, 1.d-15)

              !do p = 1, Nmg
                  !gv  = T(p) / (1.d0 - T(p))          ! γ′
                  !dgp = dT(p) / (1.d0 - T(p))**2      ! jacobiano
              !end do

        !print*, gv
    !Preparação das Splines
        call SPLGR1 (zv,Nmz)
        call SPLGR2 (gv,Nmg)
        
        
    !Montagem da Matriz

        !Pesos e absissas de Gauss-Legendre para cada variável
	      CALL legauss(0.d0,1.d0,Nz,X,dX,1.d-15)
        CALL legauss(0.d0,3.d0,Ng,Y,dY,1.d-15)
        !call setgaulag(0.d0,Ng,Y,dY)

        !CALL legauss(-1.d0,1.d0,Ng,Y,dY,1.d-15)

        !gam0 = 12.0d0
        !do i=1,Ng
          !Y(i)=2.d0/gam0*Y(i)
          !dY(i)=2.d0/gam0*dY(i)
        !end do
        !print*, Y
      

        CALL legauss(0.d0,1.d0,Nv,W,dW,1.d-15)
           
        
        zmatrix  = 0.d0
        xmatrix  = 0.d0
        !Loop para cada elemento da Matriz
        do i=1,Nmg
              g=gv(i)
              print*, i
           do j=1,Nmz
              z=zv(j)
              index1 = (j-1)*Nmg + i           !Juntei cada ponto (gi, zj) num vetor coluna de dimensão Nmg*Nmz
              !print*, j
              do k=1,Nmg
                do l=1, Nmz
                    index2 = (l-1)*Nmg + k    !Juntei cada Iteração das Splines Sg e Sz da integração em um vetor coluna de dimenção Nmg*Nmz
                    !print*, index2
                    
                    do p=1, Ng
                        do q=1, Nz
                            do r=1, Nv
    !Lado Direito
                
        ! v variando de 0 a 1
                    v = W(r)
                    dV = DW (r)
        ! gamma variando de 0 a infinito
                    gp = Y(p)
                    dgp = dY(p)

                    !gp = 3*(1.d0+Y(p))/(1.d0-Y(p))
                    !dgp = 3*(2.d0/((1.d0-Y(p))**2))*dY(p)

        ! theta z’ variando de z até 1
                     dzq=(1.d0-z)*dX(q)
                     zq=(1.d0-z)*X(q)+z
                     
                     call SPLMD1 (zv,Nmz,zq,SPLz) 
                     call SPLMD2 (gv,Nmg,gp,SPLg)
                     
        !Termos do Kernel
                        D = v*(1-v)*(1+zq)*g + v*(1+z)*gp + v*(1+z)*(1+zq)*&
                        (1-z*(1-v)-v*zq)*(kappa**2) + v*((1-v)*(1+zq)*(z**2) + &
                        v*(zq**2)*(1+z))*(m**2) + (1-v)*(1+z)*(mu**2)

                       f1_ku = - (m**2 * (v - 1.0) * v * (z - zq)) / (z + 1.0) &
                                + 0.25 * Mtot**2 * v * (zq - (zq + 1.0) * (-v*z + v*zq + z)) &
                                - v * (gp + kappa**2) &
                                + mu**2 * (v - 1.0) &
                                + (g * (v - 1.0) * v * (zq + 1.0)) / (z + 1.0)

                        Co_ku = ((v - 2.0) * (m**2 * (-v*z + v*zq + 2.0*z) + g * (v*zq + v - 2.0))) / (z + 1.0) &
                                + 0.25 * Mtot**2 * ((v*zq + v - 2.0) * ((v - 2.0)*z - v*zq) + 4.0)

                        !Numerador = Co_ku + 2*f1_ku
                        Numerador = 1.d0
            
                        zmatrix (index1, index2) = zmatrix (index1, index2)+ (1+z)**2 /&
                         (32*PI**2*(g + z**2*m**2 + (1 - z**2) * kappa**2)) * (v**2 / (D**2)) * &
                         Numerador * splg(k)*splz(l)*dzq*dgp*dv  

        ! theta z’ variando de -1 até z
                     dzq=(z+1.d0)*dx(q)
                     zq=(z+1.d0)*x(q)-1.d0
                     call SPLMD1 (zv,Nmz,zq,SPLz)

        !Termos do Kernel      
                        
                        D = v*(1-v)*(1-zq)*g + v*(1-z)*gp + v*(1-z)*(1-zq)*&
                        (1+z*(1-v)+v*zq)*(kappa**2) + v*((1-v)*(1-zq)*(z**2) +&
                         v*(zq**2)*(1-z))*(m**2) + (1-v)*(1-z)*(mu**2)
                      
                        f1_kd = - (m**2 * (v - 1.0) * v * (z - zq)) / (z - 1.0) &
                                  + 0.25 * Mtot**2 * v * ((v - 1.0) * z * (zq - 1.0) + zq * (v * (-zq) + v - 1.0)) &
                                  - v * (gp + kappa**2) &
                                  + mu**2 * (v - 1.0) &
                                  + (g * (v - 1.0) * v * (zq - 1.0)) / (z - 1.0)

                        Co_kd = ((v - 2.0) * (m**2 * (-v*z + v*zq + 2.0*z) + g * (v * (zq - 1.0) + 2.0))) / (z - 1.0) &
                                  + 0.25 * Mtot**2 * ((v * (zq - 1.0) + 2.0) * ((v - 2.0) * z - v * zq) + 4.0)

                         !Numerador = Co_kd + 2*f1_kd
                         Numerador = 1.d0

                        zmatrix (index1, index2) = zmatrix (index1, index2)+ &
                        (1-z)**2 / (32*PI**2*(g + z**2*m**2 + (1 - z**2) * kappa**2)) * (v**2 / (D**2)) * &
                        Numerador * splg(k)*splz(l)*dzq*dgp*dv
                    

                             end do
                        end do
        !Lado Esquerdo
                        call SPLMD1 (zv,Nmz,z,SPLz) 
                        xmatrix (index1, index2) = xmatrix (index1, index2) + &
                        1.0d0 / ((g +gp + ((1-z**2)*kappa**2) + m**2*z**2)**2)*splg(k)*splz(l)*dgp
                        end do 
                        
               enddo
           enddo
           enddo
    end do
           
          !Condicionar a matriz
            do i = 1, nma
                xmatrix (i,i) = xmatrix (i,i) + e
            end do

            
! --- Copiar XMATRIX -> a, ZMATRIX -> b ---
!do i=1,NMA
   !do j=1,NMA
      !b(i,j) = XMATRIX(i,j)
      !a(i,j) = ZMATRIX(i,j)
   !end do
!end do

!print*, a, b

! --- Inverter b com LAPACK ---
!call dgetrf(NMA, NMA, b, NMA, ipvt, info)
!if (info /= 0) then
   !write(6,*) 'ERRO em dgetrf, info = ', info
   !stop
!end if

!call dgetri(NMA, b, NMA, ipvt, work, LWORK, info)
!if (info /= 0) then
   !write(6,*) 'ERRO em dgetri, info = ', info
   !stop
!end if

! --- Formar c = b**(-1) * a ---
!do i=1,nma
   !do j=1,nma
      !sum = 0.d0
      !do k=1,nma
         !sum = sum + b(i,k) * a(k,j)
      !end do
      !c(i,j) = sum
   !end do
!end do

! --- Diagonalizar c ---
!call dgeev('N','V',nma,c,NMA,wr,wi,vl,1,vr,NMA,work,lwork,info)
!write(6,*) 'dgeev info = ',info

!write(6,*) nma,' eigenvalues '
!do i=1,nma
   !write(6,*) i, wr(i), wi(i)
!end do
        
        CALL DGGEV('N', 'V', NMA, ZMATRIX, NMA, XMATRIX, NMA, &
                   ALPHAR, ALPHAI, BETA, VL, NMA, VR, NMA, &
                   WORK, LWORK, INFO)

        ! Verificação de erro
        IF (INFO .NE. 0) THEN
            PRINT *, 'ERRO NO DGGEV: INFO = ', INFO
            STOP
        ENDIF

        ! O DGGEV retorna (ALPHAR + i*ALPHAI) e BETA.
        ! O autovalor real é lambda = alpha / beta.
        DO I = 1, NMA
            IF (ABS(BETA(I)) .GT. 1.D-16) THEN
                WR(I) = ALPHAR(I) / BETA(I)
                WI(I) = ALPHAI(I) / BETA(I)
            ELSE
                ! Evita divisão por zero (autovalor infinito)
                WR(I) = 1.D+16 
                WI(I) = 0.D0
            ENDIF
        END DO
      
        WRITE(10, '(A,I0,A,I0,A,I0,A)') "autovalores_Nz",nz,"_Ng",ng,"_Nv",nv,".dat"
        WRITE(12, '(A,I0,A,I0,A,I0,A)') "Numericamente.autovalores_Nz",nz,"_Ng",ng,"_Nv",nv,".dat"
          
	do I = 1, NMA
            WRITE(10,'(I4,2X,F20.12,2X,F20.12)') i, wr(i), wi(i)
   end do
      
      WRITE(10, *) ""
      Alfa = 1.0d0 / (wr(1) * 16.0d0 * PI)
      WRITE(10, '(A, F20.10)') "Valor de Alfa: ", Alfa
      WRITE(12, '(A, F20.10)') "Valor de Alfa: ", Alfa

      WRITE(10, '(9999ES16.8)') (vr(J,1), J=1, NMA)

      !Autovetores
      !Contrução dos termos cij para dps fazer o sum cij * Spline
      do j=1,Nmz
        do i = 1, Nmg
          c(i,j) = dabs(vr(i + (j-1)*Nmg, 1))
        enddo
      enddo
      
      !Printar Matriz
      DO I = 1, NMG
         WRITE(10, '(9999ES16.8)') (c(I,J), J=1, NMZ)
      END DO
     
      !Escolhendo Nplot pontos
      z_fixo = 0.5d0
      max_gamma_visualizacao = 10.d0
      N_PLOT = 300
      call SPLMD1(zv, Nmz, z_fixo, splz) ! Avalia os pesos das splines em z=0 uma única vez



          do p = 0, N_PLOT
              ! Cria uma distribuição de pontos (aqui linear, mas pode ser logarítmica)
              gamma_plot = (dble(p) / dble(N_PLOT)) * max_gamma_visualizacao
              
              ! Avalia as bases de Spline no ponto gamma_plot atual
              call SPLMD2(gv, Nmg, gamma_plot, splg)

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




      end do

      
      CLOSE(10)
      CLOSE(12)
10     FORMAT(11E12.4)
18     format(5e15.6)
20     FORMAT(A70)

       

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

    SUBROUTINE gauleg(x1,x2,x,w,n) !x1, x2 Limites de integração, x armazena as abssissas boas, w armazena os pesos, n é o número de testes

!     calculating gauss-legendre weights and abscissas - numerical recipes

      IMPLICIT REAL*8 (A-H,O-Z)
      integer :: m
      integer :: n
      integer :: i, j
      dimension x(600),w(600)
      
      eps= 3.d-15
      pi= DACOS(-1.D0)
      
      if(n.lt.1) then
      Write(6,*) 'n not a positive integer in gauleg.f'
      stop
      endif

      m=(n+1)/2 ! roots are symmetric in interval so find only half of them
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=dcos(pi*(i-.25d0)/(n+.5d0))      ! approximate the ith root
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n  ! recurrence relation for Legendre polynomial in z
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j 
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)       ! derivative of Legendre polynomial
          z1=z
          z=z1-p1/pp                      ! Newton's method to refine root
        if(dabs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z                      ! scale root to desired interval
        x(n+1-i)=xm+xl*z                  ! its symmetric counterpart
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)   ! compute weight
        w(n+1-i)=w(i)                     ! and symmetric counterpart
12    continue
      return
    END
    
    
    
    
          subroutine setgaulag(aa,n,xx,wei)
      implicit real*8(a-h,o-z) 

      dimension xx(300),wei(300),tdvr(300,300)
      dimension alf(300),bet(300),tri(300,300),aux(33*300)

      naux=33*300
      if(n.gt.300) then
      write(6,*) 'errore in setgauleg'
      stop
      endif

      do j=1,n
        alf(j)=2*j-1+aa
        bet(j)=-dsqrt(dfloat(j)*(dfloat(j)+aa))
      end do

      do i=1,n
        do j=1,n
          tri(i,j)=0.d0
        end do
        if(i.eq.1) then
          tri(1,1)=alf(1)
          tri(1,2)=bet(1)
        else if(i.eq.n) then
          tri(n,n-1)=bet(n-1)
          tri(n,n)=alf(n)
        else
          tri(i,i-1)=bet(i-1)
          tri(i,i)=alf(i)
          tri(i,i+1)=bet(i)
        end if
      end do

      call dsyev('V','U',n,tri,300,xx,aux,naux,info)    !lapack

      do i=1,n
      do j=1,n
        tdvr(i,j)=tri(j,i)
      end do
      end do

      return
    end
    
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