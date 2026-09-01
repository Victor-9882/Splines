      PROGRAM SPLINE_VETORIAL
      USE OMP_LIB
      IMPLICIT NONE

!===================== Dimensoes fixas do problema =====================
      INTEGER, PARAMETER :: NMG   = 20            !Numero de splines em gamma
      INTEGER, PARAMETER :: NMZ   = 20             !Numero de splines em z
      INTEGER, PARAMETER :: NMA   = NMG*NMZ         !Dimensao do problema de autovalores
      INTEGER, PARAMETER :: LWORK = 10*NMA          !Tamanho do buffer de trabalho do LAPACK

!===================== Malhas e bases de splines =======================
      DOUBLE PRECISION :: zv(NMZ+1), gv(300)        !Nos das malhas em z e em gamma
      DOUBLE PRECISION :: splz(NMZ), splg(NMG)      !Splines avaliadas no ponto corrente
      DOUBLE PRECISION :: XG(NMA), YG(NMA)          !Pontos de colocacao gerados por COLLOC
      INTEGER :: N_intervalZ, N_intervalG           !Numero de intervalos das malhas
      INTEGER :: NCOL                               !Pontos de colocacao por intervalo
      INTEGER :: IW                                 !Unidade de saida dos diagnosticos de malha

!===================== Matrizes do problema generalizado ===============
      DOUBLE PRECISION, ALLOCATABLE :: XMATRIX(:,:)  !Lado esquerdo (norma)
      DOUBLE PRECISION, ALLOCATABLE :: ZMATRIX(:,:)  !Lado direito (kernel)
      DOUBLE PRECISION, ALLOCATABLE :: c(:,:)        !Coeficientes cij do autovetor fundamental

!===================== Saida do DGGEV ==================================
      DOUBLE PRECISION, ALLOCATABLE :: ALPHAR(:), ALPHAI(:) !Numeradores do autovalor generalizado
      DOUBLE PRECISION, ALLOCATABLE :: BETA(:)      !Denominador do autovalor generalizado
      DOUBLE PRECISION, ALLOCATABLE :: WR(:), WI(:) !Autovalores lambda = alpha/beta
      DOUBLE PRECISION, ALLOCATABLE :: VR(:,:)      !Autovetores a direita
      DOUBLE PRECISION, ALLOCATABLE :: VL(:,:)      !Autovetores a esquerda (nao calculados)
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)      !Buffer de trabalho do LAPACK
      INTEGER :: INFO                               !Codigo de erro do LAPACK

!===================== Pontos e pesos de Gauss (COMMON) ================
      DOUBLE PRECISION :: PI
      DOUBLE PRECISION :: X(1000), DX(1000)         !Abscissas e pesos em z'
      DOUBLE PRECISION :: Y(1000), DY(1000)         !Abscissas e pesos em gamma'
      DOUBLE PRECISION :: W(1000), DW(1000)         !Abscissas e pesos em v
      DOUBLE PRECISION :: det(1000)                 !Nao utilizado; mantido pelo COMMON
      INTEGER :: Nz, Ng, Nv                         !Numero de pontos de Gauss em cada variavel
      COMMON /PARAM/ PI
      COMMON /ALPHAINT/ X, DX, Y, DY, W, DW, Nz, Ng, Nv, det

!===================== Parametros fisicos ==============================
      DOUBLE PRECISION :: Mtot                      !Massa total do estado ligado
      DOUBLE PRECISION :: m1, m2, m                 !Massas dos constituintes e media
      DOUBLE PRECISION :: mu                        !Massa da particula de troca
      DOUBLE PRECISION :: kappa                     !sqrt(m**2 - Mtot**2/4)
      DOUBLE PRECISION :: gam0                      !Escala de gamma
      DOUBLE PRECISION :: e                         !Regularizacao da diagonal de XMATRIX
      DOUBLE PRECISION :: Alfa                      !Constante de acoplamento extraida de wr(1)

!===================== Conjuntos de parametros lidos de inputs.dat =====
      INTEGER :: NPARAM                             !Numero de conjuntos a processar
      INTEGER :: Nnz(100), Nng(100), Nnv(100)       !Pontos de Gauss de cada conjunto

!===================== Variaveis do loop de montagem ===================
      INTEGER :: ii                                 !Indice do conjunto de parametros
      INTEGER :: i, j                               !Indices do ponto de colocacao (g_i, z_j)
      INTEGER :: k, l                               !Indices das splines em gamma e z
      INTEGER :: p, q, r                            !Indices de Gauss em gamma', z' e v
      INTEGER :: index1, index2                     !Linha e coluna da matriz (achatamento 2D->1D)
      DOUBLE PRECISION :: z, g                      !Ponto de colocacao corrente
      DOUBLE PRECISION :: zq, dzq                   !Abscissa e peso em z'
      DOUBLE PRECISION :: gp, dgp                   !Abscissa e peso em gamma'
      DOUBLE PRECISION :: v, dv                     !Abscissa e peso em v
      DOUBLE PRECISION :: D0                        !Denominador comum do kernel
      DOUBLE PRECISION :: Du, Dd                    !Denominadores dos ramos z'>z e z'<z
      DOUBLE PRECISION :: f1_ku, Co_ku              !Termos do numerador no ramo superior
      DOUBLE PRECISION :: Co_kd                     !Termo do numerador no ramo inferior
      DOUBLE PRECISION :: contrib_escu, contrib_escd  !Contribuicoes escalares de cada ramo
      DOUBLE PRECISION :: contrib_C0_kd, contrib_f1_kd, contrib_C0_ku, contrib_f1_ku

!===================== Diagnostico em termos.dat =======================
      !Ponto do dominio escolhido para o diagnostico
      INTEGER :: i_dbg, j_dbg, k_dbg, l_dbg, p_dbg, q_dbg, r_dbg
      LOGICAL :: printou_termos
      !Diagnostico desligado em paralelo: o bloco de WRITE(18,...) dentro do
      !laco serializa as threads e escreve concorrentemente na mesma unidade.
      LOGICAL, PARAMETER :: DEBUG_TERMOS = .TRUE.

!===================== Cronometragem OpenMP ============================
      DOUBLE PRECISION :: TSTART, TEND, TSTART_EIG, TEND_EIG
      INTEGER :: NTHREADS, NTHREADS_MKL
      INTEGER, EXTERNAL :: MKL_GET_MAX_THREADS

!===================== Acompanhamento do progresso =====================
      INTEGER :: NFEITOS                            !Pares (i,j) ja concluidos (COMPARTILHADO)
      INTEGER :: NFEITOS_LOC                        !Copia local do contador (PRIVATE)
      INTEGER :: NTOTAL                             !Total de pares = Nmg*Nmz
      INTEGER :: NPASSO                             !Imprime a cada NPASSO pares concluidos
      DOUBLE PRECISION :: TAGORA, TDECOR, TRESTA    !Tempos para a estimativa
      DOUBLE PRECISION :: PCT                       !Percentual concluido

!===================== Alocacao dos arrays grandes no heap =============
        ALLOCATE(XMATRIX(NMA,NMA), ZMATRIX(NMA,NMA), c(NMG,NMZ))
        ALLOCATE(ALPHAR(NMA), ALPHAI(NMA), BETA(NMA))
        ALLOCATE(WR(NMA), WI(NMA))
        ALLOCATE(VR(NMA,NMA), VL(NMA,2*NMA), WORK(LWORK))

        open (unit = 10, file = "autovalores.dat",STATUS="UNKNOWN")
        open (unit = 12, file = "alfa.dat",STATUS="UNKNOWN")
        open (unit = 11, file = "autovetores.dat",STATUS="UNKNOWN")
        open (UNIT = 20, FILE = "inputs.dat", STATUS="UNKNOWN")
        open (unit = 14, file = 'erros.dat', status='unknown')
        open (unit = 16, file = "coeficientes.dat",STATUS="UNKNOWN")
        open (unit = 18, file = "termos.dat", STATUS="UNKNOWN")
    
        e = 0.0001d0
        PI = DACOS(-1.D0)       !3.14159265358979323846264338

      READ(20,*) NPARAM
          DO I = 1, NPARAM
           READ(20,*) nnz(I), nng(I), nnv(I)
          END DO
    CLOSE(20)
        
     !Parâmetros
          !Massas
          Mtot = 2.3d0
          m1 = 1.0d0
          m2 = 2.3d0
          m = (m1 + m2)/2
          mu = 1.8d0
          kappa = sqrt(m**2 - 0.25*Mtot**2)

          gam0 = 10.0d0
        
          WRITE(10, '(A, I0, A, I0, A, I0)') "NMA: ", nma, " NMG: ", nmg, " NMZ: ", nmz
          WRITE(10, '(7(A, F20.10))') &
         & "Mtot: ", Mtot, " mu: ", mu, " kappa: ", kappa, " e: ", e, " gam0: ", gam0, &
         & " m1: ", m1, " m2: ", m2
          WRITE(10, *) "-----Splines Vetorial Massas Diferentes------------------------"

          WRITE(12, *) "-----Splines Vetorial Massas Diferentes------------------------"
          WRITE(12, '(A, I0, A, I0, A, I0)') "NMA: ", nma, " NMG: ", nmg, " NMZ: ", nmz
          WRITE(12, '(7(A, F20.10))') &
         & "Mtot: ", Mtot, " mu: ", mu, " kappa: ", kappa, " e: ", e, " gam0: ", gam0, &
         & " m1: ", m1, " m2: ", m2
          
            iw = 14
            N_intervalZ = (NMZ-1)/2
            N_intervalG = (NMG-1)/2
            NCOL = 2

        !Ponto do dominio (z, z', gamma, gamma', v) onde os termos serao impressos
            i_dbg = 10       !indice de gamma   (gv)
            j_dbg = 10       !indice de z       (zv)
            k_dbg = 10       !indice da spline em gamma
            l_dbg = 10      !indice da spline em z
            p_dbg = 10       !indice de Gauss em gamma'
            q_dbg = 10     !indice de Gauss em z'
            r_dbg = 10    !indice de Gauss em v

        do ii = 1, NPARAM
            print*, ii
      
          !Número de pontos de Gauss para integração em cada variável
          Nz = Nnz(ii)
          Ng = Nng (ii)
          Nv = Nnv (ii)

          printou_termos = .FALSE.

      !------------------ Cabecalho de termos.dat: parametros ------------------
          WRITE(18,*) "========================================================"
          WRITE(18,'(A,I0)') " CONJUNTO DE PARAMETROS ii = ", ii
          WRITE(18,*) " Splines Vetorial - Massas Diferentes"
          WRITE(18,*) "========================================================"
          WRITE(18,*) "--- Massas e parametros fisicos ---"
          WRITE(18,'(A,F20.10)') " Mtot  = ", Mtot
          WRITE(18,'(A,F20.10)') " m1    = ", m1
          WRITE(18,'(A,F20.10)') " m2    = ", m2
          WRITE(18,'(A,F20.10)') " m     = ", m
          WRITE(18,'(A,F20.10)') " mu    = ", mu
          WRITE(18,'(A,F20.10)') " kappa = ", kappa
          WRITE(18,'(A,F20.10)') " gam0  = ", gam0
          WRITE(18,'(A,F20.10)') " e     = ", e
          WRITE(18,'(A,F20.10)') " PI    = ", PI
          WRITE(18,*) "--- Numero de splines (malhas) ---"
          WRITE(18,'(A,I0)') " NMZ (splines em z)     = ", NMZ
          WRITE(18,'(A,I0)') " NMG (splines em gamma) = ", NMG
          WRITE(18,'(A,I0)') " NMA = NMG*NMZ          = ", NMA
          WRITE(18,'(A,I0)') " N_intervalZ            = ", N_intervalZ
          WRITE(18,'(A,I0)') " N_intervalG            = ", N_intervalG
          WRITE(18,'(A,I0)') " NCOL (colocacao)       = ", NCOL
          WRITE(18,*) "--- Pontos de Gauss (integracao) ---"
          WRITE(18,'(A,I0)') " Nz (Gauss em z')       = ", Nz
          WRITE(18,'(A,I0)') " Ng (Gauss em gamma')   = ", Ng
          WRITE(18,'(A,I0)') " Nv (Gauss em v)        = ", Nv
          WRITE(18,*) "--------------------------------------------------------"

        !--------Construção das malhas-----------------------------
        call G1D(IW,-1.d0, N_intervalZ, 1.0d0, 1.d0, X)
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

        gv(1) = 0.d0
        gv(nmg) = 3.d0

    !Preparação das Splines
        call SPLGR1 (zv,Nmz)
        call SPLGR2 (gv,Nmg)
        
        
    !Montagem da Matriz

        !Pesos e absissas de Gauss-Legendre para cada variável
        CALL legauss(0.d0,1.d0,Nz,X,dX,1.d-15)
        CALL legauss(0.d0,3.d0,Ng,Y,dY,1.d-15)
        CALL legauss(0.d0,1.d0,Nv,W,dW,1.d-15)
           
        
        zmatrix  = 0.d0
        xmatrix  = 0.d0

        NTHREADS = 1
!$      NTHREADS = OMP_GET_MAX_THREADS()
        NTHREADS_MKL = MKL_GET_MAX_THREADS()
        WRITE(*,'(A,I0,A,I0,A)') " Montagem das matrizes com ", NTHREADS, &
                                  " thread(s) OpenMP (MKL_NUM_THREADS=", &
                                  NTHREADS_MKL, ")..."
        TSTART = 0.d0
!$      TSTART = OMP_GET_WTIME()

        NFEITOS = 0
        NTOTAL  = Nmg*Nmz
        !Cerca de 20 linhas de progresso ao longo de toda a montagem
        NPASSO  = MAX(1, NTOTAL/20)

!=======================================================================
! PARALELIZACAO OpenMP
!
!   COLLAPSE(2)      : funde os lacos i e j num espaco de Nmg*Nmz
!                      iteracoes, para saturar todas as threads.
!   SCHEDULE(DYNAMIC): a CPU e hibrida (P-cores + E-cores, velocidades
!                      diferentes); o balanceamento estatico deixaria os
!                      P-cores esperando pelos E-cores na barreira final.
!   Sem REDUCTION    : cada par (i,j) escreve na linha index1 = (j-1)*Nmg+i,
!                      que e unica -> nao ha duas threads escrevendo na
!                      mesma posicao de zmatrix/xmatrix. Deixa-las SHARED
!                      evita replicar as matrizes por thread.
!   PRIVATE          : todo escalar/vetor ESCRITO dentro do laco. Atencao a
!                      SPLz e SPLg, sobrescritos por SPLMD1/SPLMD2 a cada
!                      ponto de integracao.
!=======================================================================
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP   PRIVATE(i, j, k, l, p, q, r, index1, index2, &
!$OMP           g, z, gp, dgp, v, dv, zq, dzq, &
!$OMP           SPLz, SPLg, &
!$OMP           D0, Du, Dd, f1_ku, Co_ku, Co_kd, &
!$OMP           contrib_escu, contrib_escd, &
!$OMP           NFEITOS_LOC, TAGORA, TDECOR, TRESTA, PCT) &
!$OMP   COLLAPSE(2) SCHEDULE(DYNAMIC)
        !Loop para cada elemento da Matriz
        do i=1,Nmg
           do j=1,Nmz
              g=gv(i)
              z=zv(j)
              index1 = (j-1)*Nmg + i           !Juntei cada ponto (gi, zj) num vetor coluna de dimensão Nmg*Nmz
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
                      D0 = 0.25d0*(4.d0*g + Mtot**2*(z**2 - 1.d0) + 2.d0*m1**2*(z + 1.d0) - 2.d0*m2**2*(z - 1.d0))

                      Du = 0.25d0 * ( &
                            Mtot**2 * (-v) * (z + 1.0d0) * (zq + 1.0d0) * ((v - 1.0d0) * z - v * zq + 1.0d0) &
                            + v * ( m1**2 * (z + 1.0d0) + 2.0d0 * m1 * m2 * (z + 1.0d0) + &
                                    m2**2 * (4.0d0 * v * z - 4.0d0 * v * zq - 3.0d0 * z + 4.0d0 * zq + 1.0d0) - &
                                    4.0d0 * g * (v - 1.0d0) * (zq + 1.0d0) + 4.0d0 * gp * (z + 1.0d0) ) &
                            - 4.0d0 * mu**2 * (v - 1.0d0) * (z + 1.0d0) )

                       !f1_ku = 0.25d0*(z + 1.d0)*(Mtot**2*v*(z + 1.d0)*(zq + 1.d0)*((v - 1.d0)*z - v*zq + 1.d0) &
                                 !- v*(m1**2*(z + 1.d0) + 2.d0*m1*m2*(z + 1.d0) + m2**2*(4.d0*v*z - 4.d0*v*zq &
                                 !- 3.d0*z + 4.d0*zq + 1.d0) - 4.d0*g*(v - 1.d0)*(zq + 1.d0) + 4.d0*gp*(z + 1.d0)) &
                                 !+ 4.d0*mu**2*(v - 1.d0)*(z + 1.d0))

                        Co_ku = (-1.d0)*0.25d0*(z + 1.d0)*(Mtot**2*(z + 1.d0)*((v*zq + v - 2.d0)*((v - 2.d0)*z - v*zq) + 4.d0) &
                                 + 4.d0*(v - 2.d0)*(m2**2*(-v*z + v*zq + 2.d0*z) + g*(v*zq + v - 2.d0)))

                      ! Numerador = Co_ku - 2*f1_ku
                       !Numerador = (1.d0+z)**2

                        contrib_escu = 1.d0 / (32*PI**2*D0) * (v**2 / (Du**2)) * &
                        ((1.d0 + z)**2)*splg(k)*splz(l)*dzq*dgp*dv

                        contrib_C0_ku = 1.d0 /&
                         (32*PI**2*D0) * (v**2 / (Du**2)) * &
                         Co_ku * splg(k)*splz(l)*dzq*dgp*dv

                        contrib_f1_ku = 2.d0 /&
                        (32*PI**2*D0) * (v**2 / (Du)) * &
                        (z+1)* splg(k)*splz(l)*dzq*dgp*dv

                        zmatrix (index1, index2) = zmatrix (index1, index2)+ contrib_C0_ku + contrib_f1_ku
                        !zmatrix (index1, index2) = zmatrix (index1, index2)+ contrib_escu

      !------- Impressao dos termos para o ponto escolhido do dominio -------
                    IF (DEBUG_TERMOS .AND. .NOT. printou_termos &
                        .AND. i.EQ.i_dbg .AND. j.EQ.j_dbg &
                        .AND. k.EQ.k_dbg .AND. l.EQ.l_dbg .AND. p.EQ.p_dbg &
                        .AND. q.EQ.q_dbg .AND. r.EQ.r_dbg) THEN
!$OMP CRITICAL (TERMOS_DBG)
                      WRITE(18,*) ""
                      WRITE(18,*) "=== PONTO DO DOMINIO (indices) ==="
                      WRITE(18,'(A,7(I5))') " i, j, k, l, p, q, r = ", i, j, k, l, p, q, r
                      WRITE(18,'(A,I6,A,I6)') " index1 = ", index1, "   index2 = ", index2
                      WRITE(18,*) ""
                      WRITE(18,*) "=== VARIAVEIS DO PONTO ==="
                      WRITE(18,'(A,ES24.15)') " z       = ", z
                      WRITE(18,'(A,ES24.15)') " g       = ", g
                      WRITE(18,'(A,ES24.15)') " gp      = ", gp
                      WRITE(18,'(A,ES24.15)') " dgp     = ", dgp
                      WRITE(18,'(A,ES24.15)') " v       = ", v
                      WRITE(18,'(A,ES24.15)') " dv      = ", dv
                      WRITE(18,*) ""
                      WRITE(18,*) "=== RAMO SUPERIOR: z' de z ate 1 ==="
                      WRITE(18,'(A,ES24.15)') " zq      = ", zq
                      WRITE(18,'(A,ES24.15)') " dzq     = ", dzq
                      WRITE(18,'(A,ES24.15)') " D0      = ", D0
                      WRITE(18,'(A,ES24.15)') " Du      = ", Du
                      WRITE(18,'(A,ES24.15)') " f1_ku   = ", f1_ku
                      WRITE(18,'(A,ES24.15)') " Co_ku   = ", Co_ku
                      WRITE(18,'(A,ES24.15)') " Numerador (Co_ku + 2*f1_ku) = ", Co_ku + 2*f1_ku
                      WRITE(18,'(A,ES24.15)') " splg(k) = ", splg(k)
                      WRITE(18,'(A,ES24.15)') " splz(l) = ", splz(l)
                      WRITE(18,'(A,ES24.15)') " (1+z)**2                = ", (1+z)**2
                      WRITE(18,'(A,ES24.15)') " 32*PI**2*D0             = ", 32*PI**2*D0
                      WRITE(18,'(A,ES24.15)') " v**2/Du**2              = ", v**2/(Du**2)
                      WRITE(18,'(A,ES24.15)') " contribuicao a zmatrix  = ", contrib_escu
                      WRITE(18,*) ""
                      WRITE(18,*) "--- contrib_C0_ku e contrib_f1_ku (ramo superior) ---"
                      WRITE(18,'(A,ES24.15)') " contrib_C0_ku (com splines/jacob.)    = ", contrib_C0_ku
                      WRITE(18,'(A,ES24.15)') " contrib_C0_ku (sem splines/jacob.)    = ", &
                          1.d0 / (32*PI**2*D0) * (v**2 / (Du**2)) * Co_ku
                      WRITE(18,'(A,ES24.15)') " contrib_f1_ku (com splines/jacob.)    = ", contrib_f1_ku
                      WRITE(18,'(A,ES24.15)') " contrib_f1_ku (sem splines/jacob.)    = ", &
                          2.d0 / (32*PI**2*D0) * (v**2 / (Du)) * (z+1)
                      WRITE(18,'(A,ES24.15)') " soma contrib_C0_ku + contrib_f1_ku    = ", &
                          contrib_C0_ku + contrib_f1_ku
!$OMP END CRITICAL (TERMOS_DBG)
                    END IF

                         !TesteLog =  zmatrix (index1, index2)
                          !if (isnan(TesteLog)) then
                                !write(14,*) 'AVISO: NaN! i=', i, 'j=', j, 'k=', k, 'l=', l, 'p=', p, 'q=', q
                                !write(14,*) '  g =', g, '  gp =', gp, '  z =', z, '  zq =', zq
                                !print*,"a"
                            !end if

        ! theta z’ variando de -1 até z
                     dzq=(z+1.d0)*dx(q)
                     zq=(z+1.d0)*x(q)-1.d0
                     call SPLMD1 (zv,Nmz,zq,SPLz)

        !Termos do Kernel      
                        
                        D0 = 0.25d0*(4.d0*g + Mtot**2*(z**2 - 1.d0) + 2.d0*m1**2*(z + 1.d0) - 2.d0*m2**2*(z - 1.d0))
                        Dd = 0.25d0 * ( &
                            Mtot**2 * v * (z - 1.0d0) * (zq - 1.0d0) * ((v - 1.0d0) * z - v * zq - 1.0d0) &
                            + v * ( m1**2 * (-4.0d0 * v * z + 4.0d0 * (v - 1.0d0) * zq + 3.0d0 * z + 1.0d0) &
                                    - 2.0d0 * m1 * m2 * (z - 1.0d0) - m2**2 * z + m2**2 &
                                    + 4.0d0 * (gp + g * (v - 1.0d0) * (zq - 1.0d0)) - 4.0d0 * gp * z ) &
                            + 4.0d0 * mu**2 * (v - 1.0d0) * (z - 1.0d0) )

                        !f1_kd = 0.25d0*(z - 1.d0)*(Mtot**2*v*(z - 1.d0)*(zq - 1.d0)*((v - 1.d0)*z - v*zq - 1.d0) &
                                 !+ v*(m1**2*(-4.d0*v*z + 4.d0*(v - 1.d0)*zq + 3.d0*z + 1.d0) - 2.d0*m1*m2*(z - 1.d0) &
                                 !- m2**2*z + m2**2 + 4.d0*(gp + g*(v - 1.d0)*(zq - 1.d0)) - 4.d0*gp*z) &
                                 !+ 4.d0*mu**2*(v - 1.d0)*(z - 1.d0))

                        Co_kd =(-1.d0)*0.25d0*(z - 1.d0)*(Mtot**2*(z - 1.d0)*((v*(zq - 1.d0) + 2.d0)*((v - 2.d0)*z - v*zq) + 4.d0) &
                             + 4.d0*(v - 2.d0)*(m1**2*(-v*z + v*zq + 2.d0*z) + g*(v*(zq - 1.d0) + 2.d0)))

                         !Numerador = Co_kd - 2*f1_kd
                         !Numerador = (1.d0-z)**2
                        
                        contrib_escd = 1.d0 / (32*PI**2*D0) * (v**2 / (Dd**2)) * &
                        ((1.d0 - z)**2)*splg(k)*splz(l)*dzq*dgp*dv

                        contrib_C0_kd = 1.d0 / (32*PI**2*D0) * (v**2 / (Dd**2)) * &
                        Co_kd* splg(k)*splz(l)*dzq*dgp*dv

                        contrib_f1_kd = 2.d0 / (32*PI**2*D0) * (v**2 / (Dd)) * &
                        (1.d0-z) * splg(k)*splz(l)*dzq*dgp*dv

                        zmatrix (index1, index2) = zmatrix (index1, index2)+ contrib_C0_kd + contrib_f1_kd
                        !zmatrix (index1, index2) = zmatrix (index1, index2)+ contrib_escd

      !------- Impressao dos termos para o ponto escolhido do dominio -------
                    IF (DEBUG_TERMOS .AND. .NOT. printou_termos &
                        .AND. i.EQ.i_dbg .AND. j.EQ.j_dbg &
                        .AND. k.EQ.k_dbg .AND. l.EQ.l_dbg .AND. p.EQ.p_dbg &
                        .AND. q.EQ.q_dbg .AND. r.EQ.r_dbg) THEN
!$OMP CRITICAL (TERMOS_DBG)
                      WRITE(18,*) ""
                      WRITE(18,*) "=== RAMO INFERIOR: z' de -1 ate z ==="
                      WRITE(18,'(A,ES24.15)') " zq      = ", zq
                      WRITE(18,'(A,ES24.15)') " dzq     = ", dzq
                      WRITE(18,'(A,ES24.15)') " D0      = ", D0
                      WRITE(18,'(A,ES24.15)') " Dd      = ", Dd
                      WRITE(18,'(A,ES24.15)') " Co_kd   = ", Co_kd
                      WRITE(18,'(A,ES24.15)') " Numerador (Co_kd)           = ", Co_kd
                      WRITE(18,'(A,ES24.15)') " splg(k) = ", splg(k)
                      WRITE(18,'(A,ES24.15)') " splz(l) = ", splz(l)
                      WRITE(18,'(A,ES24.15)') " (1-z)**2                = ", (1-z)**2
                      WRITE(18,'(A,ES24.15)') " 32*PI**2*D0             = ", 32*PI**2*D0
                      WRITE(18,'(A,ES24.15)') " v**2/Dd**2              = ", v**2/(Dd**2)
                      WRITE(18,'(A,ES24.15)') " contribuicao a zmatrix  = ", contrib_escd
                      WRITE(18,*) ""
                      WRITE(18,*) "--- contrib_C0_kd e contrib_f1_kd (ramo inferior) ---"
                      WRITE(18,'(A,ES24.15)') " contrib_C0_kd (com splines/jacob.)    = ", contrib_C0_kd
                      WRITE(18,'(A,ES24.15)') " contrib_C0_kd (sem splines/jacob.)    = ", &
                          1.d0 / (32*PI**2*D0) * (v**2 / (Dd**2)) * Co_kd
                      WRITE(18,'(A,ES24.15)') " contrib_f1_kd (com splines/jacob.)    = ", contrib_f1_kd
                      WRITE(18,'(A,ES24.15)') " contrib_f1_kd (sem splines/jacob.)    = ", &
                          2.d0 / (32*PI**2*D0) * (v**2 / (Dd)) * (1.d0-z)
                      WRITE(18,'(A,ES24.15)') " soma contrib_C0_kd + contrib_f1_kd    = ", &
                          contrib_C0_kd + contrib_f1_kd
                      WRITE(18,*) ""
                      WRITE(18,'(A,ES24.15)') " zmatrix(index1,index2) acumulada = ", &
                                               zmatrix(index1,index2)

                      printou_termos = .TRUE.
!$OMP END CRITICAL (TERMOS_DBG)
                    END IF


                             end do
                        end do
        !Lado Esquerdo
                        call SPLMD1 (zv,Nmz,z,SPLz)
                        xmatrix (index1, index2) = xmatrix (index1, index2) + &
                        1.0d0 / ((g +gp + ((1-z**2)*kappa**2) + m**2*z**2)**2)*splg(k)*splz(l)*dgp

      !------- Lado esquerdo no ponto escolhido do dominio -------
                    IF (DEBUG_TERMOS .AND. i.EQ.i_dbg .AND. j.EQ.j_dbg &
                        .AND. k.EQ.k_dbg .AND. l.EQ.l_dbg) THEN
!$OMP CRITICAL (TERMOS_DBG)
                      WRITE(18,*) ""
                      WRITE(18,*) "=== LADO ESQUERDO (xmatrix) ==="
                      WRITE(18,'(A,ES24.15)') " gp (ultimo p = Ng)      = ", gp
                      WRITE(18,'(A,ES24.15)') " dgp                     = ", dgp
                      WRITE(18,'(A,ES24.15)') " denominador             = ", &
                          (g + gp + ((1-z**2)*kappa**2) + m**2*z**2)**2
                      WRITE(18,'(A,ES24.15)') " splg(k)                 = ", splg(k)
                      WRITE(18,'(A,ES24.15)') " splz(l) (em z)          = ", splz(l)
                      WRITE(18,'(A,ES24.15)') " xmatrix(index1,index2)  = ", xmatrix(index1,index2)
                      WRITE(18,*) "--------------------------------------------------------"
!$OMP END CRITICAL (TERMOS_DBG)
                    END IF
                        end do

               enddo
           enddo

      !------------------ Acompanhamento do progresso ------------------
      ! Um par (i,j) acabou de ser concluido. ATOMIC garante que o
      ! incremento nao se perca quando duas threads terminam ao mesmo
      ! tempo; e muito mais barato que CRITICAL e roda uma vez por par,
      ! nao por iteracao interna, entao nao pesa no desempenho.
      ! A copia para NFEITOS_LOC evita ler a variavel compartilhada
      ! (que outra thread pode estar alterando) na hora de imprimir.
!$OMP ATOMIC CAPTURE
              NFEITOS = NFEITOS + 1
              NFEITOS_LOC = NFEITOS
!$OMP END ATOMIC

              ! So a thread que fechou um multiplo de NPASSO imprime,
              ! para nao poluir a saida com Nmg*Nmz linhas.
              IF (MOD(NFEITOS_LOC,NPASSO).EQ.0 .OR. NFEITOS_LOC.EQ.NTOTAL) THEN
                 TAGORA = 0.d0
!$               TAGORA = OMP_GET_WTIME()
                 TDECOR = TAGORA - TSTART
                 PCT    = 100.d0*DBLE(NFEITOS_LOC)/DBLE(NTOTAL)
                 !Estimativa do tempo restante por extrapolacao linear
                 TRESTA = 0.d0
                 IF (NFEITOS_LOC.GT.0) &
                    TRESTA = TDECOR*DBLE(NTOTAL-NFEITOS_LOC)/DBLE(NFEITOS_LOC)
                 WRITE(*,'(A,I6,A,I6,A,F6.2,A,F9.1,A,F9.1,A)') &
                    "   progresso: ", NFEITOS_LOC, " / ", NTOTAL, &
                    "  (", PCT, "%)   decorrido ", TDECOR, &
                    " s   resta ~", TRESTA, " s"
                 FLUSH(6)
              END IF
           enddo
    end do
!$OMP END PARALLEL DO

!$      TEND = OMP_GET_WTIME()
!$      WRITE(*,'(A,F12.3,A)') " Tempo de montagem das matrizes = ", &
!$                             TEND-TSTART, " s"

          !Condicionar a matriz
            do i = 1, nma
                xmatrix (i,i) = xmatrix (i,i) + e
            end do

                    
        NTHREADS_MKL = MKL_GET_MAX_THREADS()
        WRITE(*,'(A,I0,A)') " DGGEV com ", NTHREADS_MKL, " thread(s) MKL..."

        TSTART_EIG = 0.d0
!$      TSTART_EIG = OMP_GET_WTIME()

        CALL DGGEV('N', 'V', NMA, ZMATRIX, NMA, XMATRIX, NMA, &
                   ALPHAR, ALPHAI, BETA, VL, NMA, VR, NMA, &
                   WORK, LWORK, INFO)

!$      TEND_EIG = OMP_GET_WTIME()
!$      WRITE(*,'(A,F12.3,A)') " Tempo do DGGEV                 = ", &
!$                             TEND_EIG-TSTART_EIG, " s"

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
        WRITE(12, '(A,I0,A,I0,A,I0,A)') "Numericamente_Nz",nz,"_Ng",ng,"_Nv",nv,".dat"
          
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
          c(i,j) = (vr(i + (j-1)*Nmg, 1))
        enddo
      enddo
      
      !Printar Matriz
      DO I = 1, NMG
         WRITE(16, '(9999ES16.8)') (c(I,J), J=1, NMZ)
      END DO
     
      end do

      DEALLOCATE(XMATRIX, ZMATRIX, c)
      DEALLOCATE(ALPHAR, ALPHAI, BETA, WR, WI)
      DEALLOCATE(VR, VL, WORK)

      CLOSE(16)
      CLOSE(18)
      CLOSE(10)
      CLOSE(12)
      CLOSE(14)
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
