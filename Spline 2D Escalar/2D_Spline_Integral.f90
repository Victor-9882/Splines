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

      !Leitura Inputs
      READ(20,*) NPARAM
          DO I = 1, NPARAM
           READ(20,*) Nng(i), Nnz(i)
          END DO
    CLOSE(20)
        
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
          !print*,"a"

        WRITE(12, *) "------------------Analiticamente--------------------"
        WRITE(12, '(A, F20.10, A, F20.10, A, F20.10, A, F20.10, A, F20.10, A, F20.10)') &
         & "Mtot: ", Mtot, " m: ", m, " mu: ", mu, " kappa: ", kappa, " e", e, " gam0", gam0
        WRITE(12, '(A, I0, A, I0)') "nz: ", nz, " ng: ", ng

        
        WRITE(12, 110)
110     FORMAT(/, 'NMZ           NMG            ALFA')

        WRITE(13, '(A)') '#NMG NMZ ALFA'


      do ii = 1, NPARAM

      NMG = Nng(ii)
      NMZ = Nnz(ii)

      NMA = NMG * NMZ
      LWORK = 10 * NMA

      WRITE(10, '(A, I0, A, I0, A, I0)') "NMA: ", nma, " NMG: ", nmg, " NMZ: ", nmz
      WRITE(10, '(A, F20.10, A, F20.10, A, F20.10, A, F20.10, A, F20.10, A, F20.10)') &
         & "Mtot: ", Mtot, " m: ", m, " mu: ", mu, " kappa: ", kappa, " e", e, " gam0", gam0
          WRITE(10, *) "--------------------------------------------------"

     
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
          


          
          !print*, Mtot, m, mu, kappa

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


        !call G1D(IW,0.d0, N_intervalG, 0.999d0, 1.d0, Y)
        !call COLLOC(IW,2,N_intervalG,Y,YG)  

        call G1D(IW,0.d0, N_intervalG, 1.0d0, 3.d0, Y)
        call COLLOC(IW,2,N_intervalG,Y,YG)  

        do i=1, 2*N_intervalG
           gv(i+1)=YG(i)
        enddo

        gv(1) = 0.0000001d0
        gv(nmg) = 3.d0

        !gv(1)   =  1.d-3
        !gv(NMG) = 1.d0 - 1.d-3

       
    !Preparação das Splines

        call SPLGR1 (zv,Nmz)
        call SPLGR2 (gv,Nmg)
        
        
    !Montagem da Matriz

        !Pesos e absissas de Gauss-Legendre para cada variável
	    CALL legauss(0.d0,1.d0,Nz,X,dX,1.d-15)

        CALL legauss(0.d0,3.d0,Ng,Y,dY,1.d-15)
        !CALL legauss(0.d0,1.d0,Ng,Y,dY,1.d-15)
                
        teste = -huge(1.d0)
        zmatrix  = 0.d0
        xmatrix  = 0.d0
        !Loop para cada elemento da Matriz
        do i=1,Nmg
              g=gv(i)
              !ti = gv(i)
              !g  = f_map(ti)         ! γ físico do ponto de colocação
              !g=gv(i)

              print*, i
           do j=1,Nmz
              z=zv(j)
              index1 = (j-1)*Nmg + i           !Juntei cada ponto (gi, zj) num vetor coluna de dimensão Nmg*Nmz
              !call SPLMD1(zv, Nmz, z, splz_at)
              !print*, j
              do k=1,Nmg

           ! k_min = MAX(1, k-1)
           ! k_max = MIN(Nmg, k+1)
           ! DSK = (gv(k_max) + gv(k_min)) * 0.5d0
           ! DDK = (gv(k_max) - gv(k_min)) * 0.5d0


                do l=1, Nmz
                    index2 = (l-1)*Nmg + k    !Juntei cada Iteração das Splines Sg e Sz da integração em um vetor coluna de dimenção Nmg*Nmz
                    !print*, index2
                    
                    do p=1, Ng
                        do q=1, Nz               
    !Lado Direito
              
        ! gamma variando de 0 a infinito
                   ! tp = DDK * Y(p) + DSK
                    !tp  = Y(p)
                    !wtp = dY(p)

                    !gp = f_map(tp)
                   ! dgp = DDK*Jacobian_map(tp)*wtp
                    !dgp = Jacobian_map(tp)*wtp

                  gp = Y(p)
                  dgp = dY(p)
                    
                    !write(14,*) gp
                    
        ! theta z’ variando de z até 1
                     dzq=(1.d0-z)*dX(q)
                     zq=(1.d0-z)*X(q)+z
                     
                     call SPLMD1 (zv,Nmz,zq,SPLz) 
                     call SPLMD2 (gv,Nmg,gp,SPLg)
                     
        !Termos do Kernel
                        
                        b0 = (1.d0 + z) * (mu**2)
                        b1 = g + gp - (1.d0+ z)*(mu**2) + gp*z + g*zq + &
                        (1.d0 + zq) * ( (z**2)*(m**2) + (1.d0 - z**2)*(kappa**2) )

                        b2 = -g * (1.d0 + zq) + &
                        (z - zq) * ( (1.d0 + z)*(1.d0 + zq)*(kappa**2) - &
                        (z + zq + z*zq)*(m**2) )
                        
                                                
                        delta_sqrt =-(b1**2 - 4.d0 * b0 * b2)

                        if (delta_sqrt<0) then   
                        
                        !Karmov
                        bp = -1.d0 / (2.d0 * b2) * (b1 + sqrt(-delta_sqrt))
                        bm = -1.d0 / (2.d0 * b2) * (b1 - sqrt(-delta_sqrt))  
                            
                        IntegralV = (1.d0 / (b2**2 * (bp - bm)**3)) * &
                        ( (((bp - bm)*(2.d0*bp*bm - bp - bm)) / ((1.d0 - bp)*(1.d0 - bm))) + &
                        2.d0*bp*bm * log((bp*(1.d0 - bm)) / (bm*(1.d0 - bp))) )
                        
                        !Minha expressão
                        !termo_log1 = log( (sqrt(-delta_sqrt)- b1 - 2.d0*b2) / (sqrt(-delta_sqrt) + b1 + 2.d0*b2) )
                        !termo_log2 = log( (sqrt(-delta_sqrt) - b1) / (sqrt(-delta_sqrt) + b1) )
                        
                        !fracao1 = (2.d0 * b0 * (termo_log1 - termo_log2)) / (sqrt(-delta_sqrt)* delta_sqrt)
                        !fracao2 = (2.d0 * b0 + b1) / ((b0 + b1 + b2) * delta_sqrt)
                        
                        !IntegralV = fracao1 - fracao2
                        
                        TesteLog =  IntegralV

                
                            if (isnan(TesteLog)) then
                                write(14,*) 'AVISO: NaN! i=', i, 'j=', j, 'k=', k, 'l=', l, 'p=', p, 'q=', q
                                write(14,*) '  g =', g, '  gp =', gp, '  z =', z, '  zq =', zq
                                print*,"a"
                            end if


                            if(teste/=TesteLog) then
                              !WRITE(12, '(A, F20.10)')"Mtot: ", TesteLog
                              teste = TesteLog
                            endif

                        elseif (delta_sqrt == 0) then   
                        
                         IntegralV = (16.d0* b2**2) / (3.0d0 * b1 * (b1 + 2.d0 * b2)**3)
                          !print*, z, zq, g, gp                       

                        elseif (delta_sqrt>0) then
                            !print*, z, zq, g, gp  
                          termo1 = - (2.d0 * b0 + b1) / ((b0 + b1 + b2) * delta_sqrt)
                          num2 = 4.d0* b0 * ( atan(b1 / sqrt(delta_sqrt)) - atan((b1 + 2.d0 * b2) / sqrt(delta_sqrt)) )
                          den2 = delta_sqrt * sqrt(delta_sqrt)
                          termo2 = num2 / den2

                          IntegralV = termo1 - termo2

                          

                        endif

                        zmatrix(index1, index2) = zmatrix(index1, index2) + &
                            (1.d0/ (32.d0 * PI**2)) * &
                            (((1.d0 + z)**2) / (g + (z**2)*(m**2) + (1.d0 - z**2)*(kappa**2)))* &
                            IntegralV *splg(k) * splz(l) * dzq * dgp



                            
        ! theta z’ variando de -1 até z
                     dzq=(z+1.d0)*dx(q)
                     zq=(z+1.d0)*x(q)-1.d0
                     call SPLMD1 (zv,Nmz,zq,SPLz)

        !Termos do Kernel      
                        
                        b0 = (1.d0 - z) * (mu**2)

                        b1 = g + gp - (1.d0 - z)*(mu**2) - gp*z - g*zq + &
                            (1.d0 - zq) * ( (z**2)*(m**2) + (1.d0 - z**2)*(kappa**2) )

                        b2 = -g * (1.d0 - zq) - &
                            (z - zq) * ( (1.d0 - z)*(1.d0 - zq)*(kappa**2) + &
                            (z + zq - z*zq)*(m**2) )
                        
                        
                       

                        delta_sqrt =-(b1**2 - 4.d0 * b0 * b2)
                        
                        if (delta_sqrt<0) then   
                          
                        bp = -1.d0 / (2.d0 * b2) * (b1 + sqrt(-delta_sqrt))
                        bm = -1.d0 / (2.d0 * b2) * (b1 - sqrt(-delta_sqrt))  
                            
                        IntegralV = (1.d0 / (b2**2 * (bp - bm)**3)) * &
                        ( (((bp - bm)*(2.d0*bp*bm - bp - bm)) / ((1.d0 - bp)*(1.d0 - bm))) + &
                        2.d0*bp*bm * log((bp*(1.d0 - bm)) / (bm*(1.d0 - bp))) )

                        elseif (delta_sqrt == 0) then   
                        
                         IntegralV = (16.d0* b2**2) / (3.0d0 * b1 * (b1 + 2.d0 * b2)**3)

                        elseif (delta_sqrt>0) then

                          termo1 = - (2.d0 * b0 + b1) / ((b0 + b1 + b2) * delta_sqrt)
                          num2 = 4.d0* b0 * ( atan(b1 / sqrt(delta_sqrt)) - atan((b1 + 2.d0 * b2) / sqrt(delta_sqrt)) )
                          den2 = delta_sqrt * sqrt(delta_sqrt)
                          termo2 = num2 / den2

                          IntegralV = termo1 - termo2

                        
                        endif


                        zmatrix(index1, index2) = zmatrix(index1, index2) + &
                            (1.d0/ (32.d0 * PI**2)) * &
                            (((1.d0 - z)**2) / (g + (z**2)*(m**2) + (1.d0 - z**2)*(kappa**2)))* &
                            IntegralV * splg(k) * splz(l) * dzq * dgp
                           
                        enddo
        !Lado Esquerdo
                        call SPLMD1 (zv,Nmz,z,SPLz) 
                        xmatrix (index1, index2) = xmatrix (index1, index2) + &
                        1.0d0 / ((g +gp + ((1.d0 -z**2)*kappa**2) + m**2*z**2)**2)*splg(k)*splz(l)*dgp
                        end do 
                        
               enddo
           enddo
           enddo
    end do
           
          !print*, zmatrix
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
          
	do I = 1, NMA
            WRITE(10,'(I4,2X,F20.12,2X,F20.12)') i, wr(i), wi(i)
   end do
      
      WRITE(10, *) ""
      Alfa = 1.0d0 / (wr(1) * 16.0d0 * PI)
      WRITE(10, '(A, F20.10)') "Valor de Alfa: ", Alfa
      WRITE(12, 120) NMZ, NMG, Alfa
120        FORMAT(I5, 9X, I5, 11X, F10.8)
      WRITE(13, 130) NMG, NMZ, Alfa
130    FORMAT(I4, 1X, I4, 1X, F22.15)


      !WRITE(10, '(9999ES16.8)') (vr(J,1), J=1, NMA)

      !Autovetores
      !Contrução dos termos cij para dps fazer o sum cij * Spline
      do j=1,Nmz
        do i = 1, Nmg
          c(i,j) = (vr(i + (j-1)*Nmg, 1))
        enddo
      enddo
      
      !Printar Matriz
      DO I = 1, NMG
         WRITE(16, '(9999ES16.8)') (c(I,J), J=1, NMZ)
      END DO
     

!Plot de função x gamma para z fixo
      z_fixo = 0.8d0
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


      ! 2. CÁLCULO DO NUMERADOR E PLOT: Psi(gamma_plot, 0.2)
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
      ! PLOT DA FUNÇÃO PSI NORMALIZADA: Psi(0.54, z) / Psi(0, 0.5)
      ! -------------------------------------------------------------
      open(unit = 17, file = "plot_psi_z.dat", STATUS="UNKNOWN")

      ! Numerador tem gamma fixo em 0.54
      gamma_num = 0.54d0
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
    
    


    !Mapeamento

DOUBLE PRECISION FUNCTION inverse_map(gamma)
    DOUBLE PRECISION :: gamma
    
    ! O inverso exato da função f_map = 0.5d0*dlog((1.d0+t)/(1.d0-t))
    inverse_map = dtanh(gamma)
    
END FUNCTION

DOUBLE PRECISION FUNCTION f_map(t)
    DOUBLE PRECISION :: t, gam0
    gam0 = 10.0d0
    ! Mapeamento atual do seu código
   ! f_map = gam0 * (1.d0+t) / (1.d0 - t + 1.d-12)
    f_map=0.5d0*dlog((1.d0+t)/(1.d0-t))
END FUNCTION

DOUBLE PRECISION FUNCTION Jacobian_map(t)
    DOUBLE PRECISION :: t, gam0
    gam0 = 10.0d0
    ! Derivada dq/dt do mapeamento escolhido
    !Jacobian_map = 2* gam0 / ((1.d0 - t + 1.d-12)**2)
    Jacobian_map= 1/(1.d0-t*t)   
END FUNCTION
    
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
