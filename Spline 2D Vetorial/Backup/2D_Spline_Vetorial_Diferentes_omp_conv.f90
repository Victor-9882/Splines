!=======================================================================
!  SPLINE_VETORIAL_CONV
!
!  Versao voltada para a velocidade de convergencia, derivada de
!  2D_Spline_Vetorial_Diferentes_omp.f90.
!
!  DUAS mudancas estruturais ativas:
!
!  (1) QUADRATURA EM z' SUBDIVIDIDA POR INTERVALO DA MALHA
!      O integrando tem ramos distintos separados em z'=z e a base tem
!      derivada terceira descontinua nos nos. Gauss-Legendre so converge
!      rapido em integrando suave, entao aplicar um unico painel sobre
!      [z,1] (ou [-1,z]) destroi a taxa. Cada ramo e quebrado nos nos da
!      malha em z, com Nz pontos de Gauss por subintervalo, e o no que
!      contem z e cortado exatamente em z.
!      (Equivale ao teste SING / particao em elm_droite_mod do frances.)
!
!      EFEITO MEDIDO: com a malha fixa, variar Nz de 8 para 32 muda
!      lambda_1 apenas na 10a casa decimal. A direcao z -- a queixa
!      original -- ficou saturada: Nz=8 ja basta, onde a versao anterior
!      pedia dezenas de splines.
!
!  (2) BASE DE HERMITE CUBICA COM SUPORTE LOCAL
!      As splines cubicas naturais antigas (SPLGR1/SPLMD1) tinham
!      suporte global e impunham segunda derivada nula em z=+-1, o que
!      e fisicamente errado aqui (o kernel tem fator (1+-z)**2). A base
!      agora e a de Hermite de SJX: duas funcoes por no (valor e
!      derivada), cada uma nula fora dos dois intervalos vizinhos.
!      Nao impoe condicao de contorno espuria e refina localmente.
!      Com NIZ intervalos ha 2*(NIZ+1) funcoes de base em z.
!      Verificada: forma particao da unidade e interpola x**3 exatamente.
!
!  GAMMA: TRUNCAMENTO MANUAL EM 3, COMO NO CODIGO ORIGINAL
!      O mapeamento gamma = gam0*x/(1-x) que levava o dominio ao
!      infinito foi REVERTIDO a pedido. A malha em gamma volta a ser
!      [0, GMAX] com GMAX=3, e a quadratura em gamma' idem, exatamente
!      como em 2D_Spline_Vetorial_Diferentes_omp.f90.
!
!      ATENCAO -- custo do truncamento, MEDIDO neste arquivo:
!      gamma varia de fato em [0,infinito) e o peso de Nakanishi cai
!      como 1/(gamma+gamma'+...)**2, entao a cauda gamma' > GMAX nao e
!      desprezivel. Variando so GMAX (NIZ=8, NIG=6, m1=m2=1, Mtot=1.99,
!      Nz=Ng=Nv=12):
!         GMAX=  3 -> lambda_1 = 0.03985
!         GMAX=  6 -> lambda_1 = 0.04495
!         GMAX= 12 -> lambda_1 = 0.04992
!         GMAX= 25 -> lambda_1 = 0.05503
!         GMAX= 50 -> lambda_1 = 0.06002
!         GMAX=100 -> lambda_1 = 0.06616
!      Ou seja +66% de 3 para 100, ainda subindo sem saturar. Refinar a
!      malha converge para um valor que depende de GMAX, nao para o
!      valor fisico: com GMAX=3, lambda_1 nao e resultado final.
!      Para varrer isso basta mudar GMAX abaixo.
!      A variavel gam0 volta a nao ter uso (mantida so para registro
!      nos arquivos de saida, como no original).
!
!  Compilacao (igual a da versao anterior):
!    gfortran -O2 -fopenmp -o spline_conv \
!             2D_Spline_Vetorial_Diferentes_omp_conv.f90 -llapack -lblas
!=======================================================================
      PROGRAM SPLINE_VETORIAL_CONV
      USE OMP_LIB
      IMPLICIT NONE

!===================== Discretizacao ===================================
!     Numero de INTERVALOS das malhas. A base de Hermite tem 2 graus de
!     liberdade por no, logo 2*(N+1) funcoes para N intervalos.
      INTEGER, PARAMETER :: NIZ = 15                 !Intervalos em z
      INTEGER, PARAMETER :: NIG = 15                 !Intervalos em x (gamma mapeado)
      INTEGER, PARAMETER :: NMZ = 2*(NIZ+1)         !Funcoes de base em z
      INTEGER, PARAMETER :: NMG = 2*(NIG+1)         !Funcoes de base em gamma
      INTEGER, PARAMETER :: NMA = NMG*NMZ           !Dimensao do problema
      INTEGER, PARAMETER :: LWORK = 10*NMA

!===================== Malhas ==========================================
!     Nos das malhas, indexados de 0 a N como no codigo frances.
      DOUBLE PRECISION :: zv(0:NIZ)                 !Nos em z
      DOUBLE PRECISION :: gv(0:NIG)                 !Nos em gamma
      DOUBLE PRECISION :: zc(2*NIZ), gc(2*NIG)      !Pontos de colocacao de Gauss
      DOUBLE PRECISION :: zcol(NMZ), gcol(NMG)      !Colocacao completa (nos + Gauss)
      INTEGER, PARAMETER :: NCOL = 2                !Pontos de colocacao por intervalo
      INTEGER :: IW

!===================== Matrizes ========================================
      DOUBLE PRECISION, ALLOCATABLE :: XMATRIX(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: ZMATRIX(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: c(:,:)

!===================== Saida do DGGEV ==================================
      DOUBLE PRECISION, ALLOCATABLE :: ALPHAR(:), ALPHAI(:)
      DOUBLE PRECISION, ALLOCATABLE :: BETA(:)
      DOUBLE PRECISION, ALLOCATABLE :: WR(:), WI(:)
      DOUBLE PRECISION, ALLOCATABLE :: VR(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: VL(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
      INTEGER :: INFO

!===================== Quadratura ======================================
      DOUBLE PRECISION :: PI
      DOUBLE PRECISION :: XQ(1000), DXQ(1000)       !Gauss em [0,1] para z'
      DOUBLE PRECISION :: YQ(1000), DYQ(1000)       !Gauss em [-1,1] para x'
      DOUBLE PRECISION :: WQ(1000), DWQ(1000)       !Gauss em [0,1] para v
      INTEGER :: Nz, Ng, Nv

!===================== Parametros fisicos ==============================
      DOUBLE PRECISION :: Mtot, m1, m2, m, mu, kappa
      DOUBLE PRECISION :: gam0                      !Sem uso (mantida para registro na saida)
      DOUBLE PRECISION :: e                         !Regularizacao da diagonal
      DOUBLE PRECISION :: Alfa

!     Limite superior da malha e da quadratura em gamma.
!     Truncamento manual, como no codigo original. Ver a nota sobre
!     gamma no cabecalho: a cauda gamma > GMAX e descartada.
      DOUBLE PRECISION, PARAMETER :: GMAX = 3.d0

      CHARACTER :: kernel_type

      INTEGER :: NPARAM
      INTEGER :: Nnz(100), Nng(100), Nnv(100)
      INTEGER :: ii, i, j

      DOUBLE PRECISION :: TSTART_EIG, TEND_EIG
      INTEGER :: NTHREADS_MKL


        ALLOCATE(XMATRIX(NMA,NMA), ZMATRIX(NMA,NMA), c(NMG,NMZ))
        ALLOCATE(ALPHAR(NMA), ALPHAI(NMA), BETA(NMA))
        ALLOCATE(WR(NMA), WI(NMA))
        ALLOCATE(VR(NMA,NMA), VL(NMA,2*NMA), WORK(LWORK))

        OPEN (unit = 10, file = "autovalores_conv.dat", STATUS="UNKNOWN")
        OPEN (unit = 12, file = "alfa_conv.dat",        STATUS="UNKNOWN")
        OPEN (unit = 20, file = "inputs.dat",           STATUS="OLD")
        OPEN (unit = 14, file = "erros_conv.dat",       STATUS="UNKNOWN")
        OPEN (unit = 16, file = "coeficientes_conv.dat",STATUS="UNKNOWN")

        e  = 0.000001d0
        PI = DACOS(-1.D0)

        READ(20,*) NPARAM
        READ(20,*) kernel_type
        DO I = 1, NPARAM
          READ(20,*) nnz(I), nng(I), nnv(I)
        END DO
        CLOSE(20)

!===================== Parametros fisicos ==============================
        Mtot  = 2.3d0
        m1    = 1.0d0
        m2    = 2.3d0
        m     = (m1 + m2)/2
        mu    = 1.8d0
        kappa = sqrt(m**2 - 0.25*Mtot**2)

!       Sem efeito no calculo; so aparece nos arquivos de saida, como
!       acontecia no codigo original.
        gam0 = 10.0d0

        NTHREADS_MKL = 1
!$      NTHREADS_MKL = OMP_GET_MAX_THREADS()

        WRITE(10, '(A, I0, A, I0, A, I0)') "NMA: ", NMA, " NMG: ", NMG, " NMZ: ", NMZ
        WRITE(10, '(A, I0, A, I0)')        "Intervalos  NIZ: ", NIZ, "  NIG: ", NIG
        WRITE(10, '(7(A, F20.10))') &
       & "Mtot: ", Mtot, " mu: ", mu, " kappa: ", kappa, " e: ", e, " gam0: ", gam0, &
       & " m1: ", m1, " m2: ", m2
        WRITE(10, '(A,A)') " kernel_type: ", kernel_type
        WRITE(10, *) "-----Splines Vetorial - base Hermite, gamma truncado em GMAX-----"

        WRITE(12, *) "-----Splines Vetorial - base Hermite, gamma truncado em GMAX-----"
        WRITE(12, '(A, I0, A, I0, A, I0)') "NMA: ", NMA, " NMG: ", NMG, " NMZ: ", NMZ
        WRITE(12, '(7(A, F20.10))') &
       & "Mtot: ", Mtot, " mu: ", mu, " kappa: ", kappa, " e: ", e, " gam0: ", gam0, &
       & " m1: ", m1, " m2: ", m2
        WRITE(12, '(A,A)') " kernel_type: ", kernel_type

        iw = 14

!===================== Malha em z ======================================
!     Malha uniforme em [-1,1]. As bordas sao nos exatos; a base de
!     Hermite nao impoe nada sobre a segunda derivada ali, ao contrario
!     da spline natural, entao podemos usar -1 e +1 sem artificio.
        CALL G1D(IW, -1.d0, NIZ, 1.0d0, 1.d0, zv)
        CALL COLLOC(IW, NCOL, NIZ, zv, zc)

!===================== Malha em gamma ==================================
!     Malha uniforme em [0, GMAX], direta em gamma (sem mapeamento).
        CALL G1D(IW, 0.d0, NIG, 1.0d0, GMAX, gv)
        CALL COLLOC(IW, NCOL, NIG, gv, gc)

!===================== Pontos de colocacao =============================
!     Sao NMZ = 2*(NIZ+1) pontos: os dois nos extremos mais os 2*NIZ
!     pontos de Gauss internos, na ordem do codigo frances.
        CALL BUILD_COLLOC(zv, zc, NIZ, NCOL, zcol, NMZ)
        CALL BUILD_COLLOC(gv, gc, NIG, NCOL, gcol, NMG)

        WRITE(14,*) ""
        WRITE(14,*) "=== Pontos de colocacao em z ==="
        DO i = 1, NMZ
          WRITE(14,'(I5,2X,ES20.12)') i, zcol(i)
        END DO
        WRITE(14,*) ""
        WRITE(14,*) "=== Pontos de colocacao em gamma ==="
        DO i = 1, NMG
          WRITE(14,'(I5,2X,ES20.12)') i, gcol(i)
        END DO

      DO ii = 1, NPARAM
        PRINT *, "Conjunto de parametros ", ii

        Nz = Nnz(ii)
        Ng = Nng(ii)
        Nv = Nnv(ii)

!===================== Pesos de Gauss ==================================
!       XQ em [0,1]: usado para mapear cada SUBINTERVALO de z'.
!       YQ em [-1,1]: usado para mapear cada subintervalo de x'.
!       WQ em [0,1]: variavel v, sem mudanca.
        CALL legauss(0.d0,  1.d0, Nz, XQ, DXQ, 1.d-15)
        CALL legauss(-1.d0, 1.d0, Ng, YQ, DYQ, 1.d-15)
        CALL legauss(0.d0,  1.d0, Nv, WQ, DWQ, 1.d-15)

        CALL BUILD_MATRICES(zv, gv, zcol, gcol, NIZ, NIG, NMZ, NMG, NMA, &
                            m, m1, m2, mu, kappa, Mtot, PI, &
                            XQ, DXQ, YQ, DYQ, WQ, DWQ, Nz, Ng, Nv, &
                            ZMATRIX, XMATRIX, kernel_type)

!       Condicionamento da matriz da norma
        DO i = 1, NMA
          XMATRIX(i,i) = XMATRIX(i,i) + e
        END DO

        WRITE(*,'(A,I0,A)') " DGGEV com ", NTHREADS_MKL, " thread(s)..."

        TSTART_EIG = 0.d0
!$      TSTART_EIG = OMP_GET_WTIME()

        CALL DGGEV('N', 'V', NMA, ZMATRIX, NMA, XMATRIX, NMA, &
                   ALPHAR, ALPHAI, BETA, VL, NMA, VR, NMA, &
                   WORK, LWORK, INFO)

!$      TEND_EIG = OMP_GET_WTIME()
!$      WRITE(*,'(A,F12.3,A)') " Tempo do DGGEV                 = ", &
!$                             TEND_EIG-TSTART_EIG, " s"

        IF (INFO .NE. 0) THEN
          PRINT *, 'ERRO NO DGGEV: INFO = ', INFO
          STOP
        ENDIF

        DO I = 1, NMA
          IF (ABS(BETA(I)) .GT. 1.D-16) THEN
            WR(I) = ALPHAR(I) / BETA(I)
            WI(I) = ALPHAI(I) / BETA(I)
          ELSE
            WR(I) = 1.D+16
            WI(I) = 0.D0
          ENDIF
        END DO

        WRITE(10, '(A,I0,A,I0,A,I0,A)') "autovalores_Nz",Nz,"_Ng",Ng,"_Nv",Nv,".dat"
        WRITE(12, '(A,I0,A,I0,A,I0,A)') "Numericamente_Nz",Nz,"_Ng",Ng,"_Nv",Nv,".dat"

        DO I = 1, NMA
          WRITE(10,'(I4,2X,F20.12,2X,F20.12)') i, wr(i), wi(i)
        END DO

        WRITE(10, *) ""
        Alfa = 1.0d0 / (wr(1) * 16.0d0 * PI)
        WRITE(10, '(A, F20.10)') "Valor de Alfa: ", Alfa
        WRITE(12, '(A, F20.10)') "Valor de Alfa: ", Alfa

        WRITE(10, '(9999ES16.8)') (vr(J,1), J=1, NMA)

        DO j = 1, NMZ
          DO i = 1, NMG
            c(i,j) = vr(i + (j-1)*NMG, 1)
          END DO
        END DO

        DO I = 1, NMG
          WRITE(16, '(9999ES16.8)') (c(I,J), J=1, NMZ)
        END DO

      END DO

      DEALLOCATE(XMATRIX, ZMATRIX, c)
      DEALLOCATE(ALPHAR, ALPHAI, BETA, WR, WI)
      DEALLOCATE(VR, VL, WORK)

      CLOSE(16)
      CLOSE(10)
      CLOSE(12)
      CLOSE(14)

      END PROGRAM SPLINE_VETORIAL_CONV


!=======================================================================
!  BUILD_COLLOC
!  Monta a lista dos NM = 2*(N+1) pontos de colocacao de uma malha:
!  o no inicial, os 2*N pontos de Gauss internos e o no final. E a
!  mesma ordem usada nos lacos DO 1 / DO 2 de elm_droite_mod, e casa
!  com a contagem de funcoes da base de Hermite.
!=======================================================================
      SUBROUTINE BUILD_COLLOC(gv, gc, N, NCOL, col, NM)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, NCOL, NM
      DOUBLE PRECISION, INTENT(IN)  :: gv(0:N), gc(NCOL*N)
      DOUBLE PRECISION, INTENT(OUT) :: col(NM)
      INTEGER :: i

      col(1) = gv(0)
      DO i = 1, NCOL*N
        col(i+1) = gc(i)
      END DO
      col(NCOL*N+2) = gv(N)

      END SUBROUTINE BUILD_COLLOC


!=======================================================================
!  SJX_H — base de Hermite cubica com suporte local
!
!  Valor da funcao de base J (0 <= J <= 2*N+1) no ponto X, sobre a malha
!  G(0:N). Cada no i tem duas funcoes:
!     J = 2*i    -> vale 1 em x_i, derivada 0   (funcao de valor)
!     J = 2*i+1  -> vale 0 em x_i, derivada 1   (funcao de derivada)
!  Ambas sao nulas fora de [x_{i-1}, x_{i+1}], o que torna as matrizes
!  esparsas e desacopla o erro entre regioes distantes do dominio.
!
!  Transcrita de SJX do codigo frances (bsminf.f90), que a validou
!  contra a spline global com acordo melhor que 1.D-15.
!=======================================================================
      SUBROUTINE SJX_H(G, N, X, J, S)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, J
      DOUBLE PRECISION, INTENT(IN)  :: G(0:N), X
      DOUBLE PRECISION, INTENT(OUT) :: S
      INTEGER :: I
      DOUBLE PRECISION :: XI, XIM1, XIP1, HI, HIP1

      S = 0.D0
      I = J/2
      XI = G(I)

      IF (I .EQ. 0) THEN
        XIM1 = G(0)
        XIP1 = G(1)
        IF (X .LT. XIM1 .OR. X .GT. XIP1) RETURN
        HIP1 = XIP1 - XI
        IF (J .EQ. 2*I)   S = (2.D0/HIP1**3)*(X-XI+0.5D0*HIP1)*(X-XIP1)**2
        IF (J .EQ. 2*I+1) S = (X-XI)*(X-XIP1)**2/HIP1**2
        RETURN
      ENDIF

      IF (I .EQ. N) THEN
        XIM1 = G(N-1)
        XIP1 = G(N)
        IF (X .LT. XIM1 .OR. X .GT. XIP1) RETURN
        HI = XI - XIM1
        IF (J .EQ. 2*I)   S = -(2.D0/HI**3)*(X-XIM1)**2*(X-XI-0.5D0*HI)
        IF (J .EQ. 2*I+1) S = (X-XIM1)**2*(X-XI)/HI**2
        RETURN
      ENDIF

      XIM1 = G(I-1)
      XIP1 = G(I+1)
      IF (X .LT. XIM1 .OR. X .GE. XIP1) RETURN

      IF (X .LT. XI) THEN
        HI = XI - XIM1
        IF (J .EQ. 2*I)   S = -(2.D0/HI**3)*(X-XIM1)**2*(X-XI-0.5D0*HI)
        IF (J .EQ. 2*I+1) S = (X-XIM1)**2*(X-XI)/HI**2
      ELSE
        HIP1 = XIP1 - XI
        IF (J .EQ. 2*I)   S = (2.D0/HIP1**3)*(X-XI+0.5D0*HIP1)*(X-XIP1)**2
        IF (J .EQ. 2*I+1) S = (X-XI)*(X-XIP1)**2/HIP1**2
      ENDIF

      END SUBROUTINE SJX_H


!=======================================================================
!  SUPPORT_RANGE
!  Devolve os extremos [a,b] do suporte da funcao de base J.
!  Usado para restringir a quadratura ao trecho onde a base nao e nula,
!  que e o que torna a base local uma economia de verdade e nao so uma
!  troca de formula.
!=======================================================================
      SUBROUTINE SUPPORT_RANGE(G, N, J, IA, IB)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: N, J
      DOUBLE PRECISION, INTENT(IN) :: G(0:N)
      INTEGER, INTENT(OUT) :: IA, IB
      INTEGER :: I

      I  = J/2
      IA = I - 1
      IB = I + 1
      IF (IA .LT. 0) IA = 0
      IF (IB .GT. N) IB = N

      END SUBROUTINE SUPPORT_RANGE


!=======================================================================
!  BUILD_MATRICES
!
!  Monta ZMATRIX (kernel) e XMATRIX (norma) por colocacao.
!
!  Estrutura dos lacos, do lado de fora para dentro:
!     (i,j)   ponto de colocacao (x_i, z_j)  -> linha index1
!     (k,l)   funcao de base    (k em x, l em z) -> coluna index2
!     integracao em x' (gamma'), z' e v
!
!  As duas melhorias de quadratura aparecem aqui:
!   - a integral em gamma' cobre so o suporte da base k, subintervalo
!     a subintervalo;
!   - a integral em z' e quebrada nos nos da malha contidos no ramo, e
!     o no que contem z e cortado em z.
!=======================================================================
      SUBROUTINE BUILD_MATRICES(zv, gv, zcol, gcol, NIZ, NIG, NMZ, NMG, NMA, &
                                m, m1, m2, mu, kappa, Mtot, PI, &
                                XQ, DXQ, YQ, DYQ, WQ, DWQ, Nz, Ng, Nv, &
                                ZMATRIX, XMATRIX, kernel_type)
      USE OMP_LIB
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NIZ, NIG, NMZ, NMG, NMA, Nz, Ng, Nv
      DOUBLE PRECISION, INTENT(IN) :: zv(0:NIZ), gv(0:NIG)
      DOUBLE PRECISION, INTENT(IN) :: zcol(NMZ), gcol(NMG)
      DOUBLE PRECISION, INTENT(IN) :: m, m1, m2, mu, kappa, Mtot, PI
      DOUBLE PRECISION, INTENT(IN) :: XQ(Nz), DXQ(Nz)
      DOUBLE PRECISION, INTENT(IN) :: YQ(Ng), DYQ(Ng)
      DOUBLE PRECISION, INTENT(IN) :: WQ(Nv), DWQ(Nv)
      DOUBLE PRECISION, INTENT(OUT) :: ZMATRIX(NMA,NMA), XMATRIX(NMA,NMA)
      CHARACTER, INTENT(IN) :: kernel_type

      INTEGER :: i, j, k, l, p, q, r
      INTEGER :: index1, index2
      INTEGER :: iaG, ibG, iaZ, ibZ
      INTEGER :: isubG                              !subintervalo em gamma' (externo)
      INTEGER :: isubZ, nsub                        !subintervalo em z' (interno)
      DOUBLE PRECISION :: z, g
      DOUBLE PRECISION :: gp, dgp
      DOUBLE PRECISION :: zq, dzq, v, dv
      DOUBLE PRECISION :: aaG, bbG, dsx, ddx        !limites em gamma' (externo)
      DOUBLE PRECISION :: aa, bb                    !limites em z' (interno)
      DOUBLE PRECISION :: SPLz, SPLg, SPLz_col
      DOUBLE PRECISION :: peso, acc_z, acc_x
      DOUBLE PRECISION :: zbrk(0:NIZ+2)             !Nos do ramo + corte em z
      LOGICAL :: SING_U, SING_D
      DOUBLE PRECISION, PARAMETER :: TOL_SING = 1.D-10

      DOUBLE PRECISION, EXTERNAL :: KERNEL_UPPER_ESC, KERNEL_LOWER_ESC
      DOUBLE PRECISION, EXTERNAL :: KERNEL_UPPER_C0,  KERNEL_LOWER_C0
      DOUBLE PRECISION, EXTERNAL :: KERNEL_UPPER_F1,  KERNEL_LOWER_F1

      DOUBLE PRECISION :: TSTART, TEND
      INTEGER :: NTHREADS
      INTEGER :: NFEITOS, NFEITOS_LOC, NTOTAL, NPASSO
      DOUBLE PRECISION :: TAGORA, TDECOR, TRESTA, PCT

        ZMATRIX = 0.d0
        XMATRIX = 0.d0

        NTHREADS = 1
!$      NTHREADS = OMP_GET_MAX_THREADS()

        WRITE(*,'(A,I0,A)') " Montagem das matrizes com ", NTHREADS, " thread(s) OpenMP..."
        TSTART = 0.d0
!$      TSTART = OMP_GET_WTIME()

        NFEITOS = 0
        NTOTAL  = NMG*NMZ
        NPASSO  = MAX(1, NTOTAL/20)

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP   PRIVATE(i, j, k, l, p, q, r, index1, index2, &
!$OMP           z, g, gp, dgp, zq, dzq, v, dv, &
!$OMP           aaG, bbG, dsx, ddx, aa, bb, nsub, isubG, isubZ, zbrk, &
!$OMP           iaG, ibG, iaZ, ibZ, &
!$OMP           SPLz, SPLg, SPLz_col, peso, acc_z, acc_x, &
!$OMP           SING_U, SING_D, &
!$OMP           NFEITOS_LOC, TAGORA, TDECOR, TRESTA, PCT) &
!$OMP   COLLAPSE(2) SCHEDULE(DYNAMIC)
        DO i = 1, NMG
          DO j = 1, NMZ

            g = gcol(i)
            z = zcol(j)

!           Em z = +-1 o numerador e o denominador de um dos ramos se
!           anulam juntos (singularidade removivel). Como -1 e +1 sao
!           nos exatos da malha, a divisao literal ali e instavel;
!           pulamos o ramo, como faz Kernel_droite no codigo frances.
            SING_U = (DABS(z + 1.d0) .LT. TOL_SING)
            SING_D = (DABS(z - 1.d0) .LT. TOL_SING)

            index1 = (j-1)*NMG + i

            DO k = 1, NMG
              DO l = 1, NMZ

                index2 = (l-1)*NMG + k

!               Suporte das funcoes de base (indices 0-based em SJX_H)
                CALL SUPPORT_RANGE(gv, NIG, k-1, iaG, ibG)
                CALL SUPPORT_RANGE(zv, NIZ, l-1, iaZ, ibZ)

                acc_z = 0.d0
                acc_x = 0.d0

!===================================================================
!               LADO DIREITO: integral em gamma', z' e v
!===================================================================
!               Integral em gamma' restrita ao suporte da base k, um
!               subintervalo da malha por vez.
                DO isubG = iaG, ibG-1
                  aaG = gv(isubG)
                  bbG = gv(isubG+1)
                  dsx = 0.5d0*(bbG+aaG)
                  ddx = 0.5d0*(bbG-aaG)

                  DO p = 1, Ng
                    gp = ddx*YQ(p) + dsx
                    CALL SJX_H(gv, NIG, gp, k-1, SPLg)
                    IF (SPLg .EQ. 0.d0) CYCLE

                    dgp = ddx*DYQ(p)

!-------------------- Ramo superior: z' de z ate 1 --------------------
                    IF (.NOT. SING_U) THEN
!                     Quebra [z,1] nos nos da malha em z contidos no ramo.
                      CALL SPLIT_RANGE(zv, NIZ, z, 1.d0, zbrk, nsub)
                      DO isubZ = 1, nsub
                        aa = zbrk(isubZ-1)
                        bb = zbrk(isubZ)
                        IF (bb - aa .LT. 1.d-14) CYCLE
                        DO q = 1, Nz
                          zq  = (bb-aa)*XQ(q) + aa
                          dzq = (bb-aa)*DXQ(q)
                          CALL SJX_H(zv, NIZ, zq, l-1, SPLz)
                          IF (SPLz .EQ. 0.d0) CYCLE
                          DO r = 1, Nv
                            v  = WQ(r)
                            dv = DWQ(r)
                            peso = SPLg*SPLz*dzq*dgp*dv
                            IF (kernel_type .EQ. 'v') THEN
                              acc_z = acc_z &
                                - KERNEL_UPPER_C0(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)*peso &
                                + KERNEL_UPPER_F1(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)*peso
                            ELSE
                              acc_z = acc_z &
                                + KERNEL_UPPER_ESC(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)*peso
                            END IF
                          END DO
                        END DO
                      END DO
                    END IF

!-------------------- Ramo inferior: z' de -1 ate z --------------------
                    IF (.NOT. SING_D) THEN
                      CALL SPLIT_RANGE(zv, NIZ, -1.d0, z, zbrk, nsub)
                      DO isubZ = 1, nsub
                        aa = zbrk(isubZ-1)
                        bb = zbrk(isubZ)
                        IF (bb - aa .LT. 1.d-14) CYCLE
                        DO q = 1, Nz
                          zq  = (bb-aa)*XQ(q) + aa
                          dzq = (bb-aa)*DXQ(q)
                          CALL SJX_H(zv, NIZ, zq, l-1, SPLz)
                          IF (SPLz .EQ. 0.d0) CYCLE
                          DO r = 1, Nv
                            v  = WQ(r)
                            dv = DWQ(r)
                            peso = SPLg*SPLz*dzq*dgp*dv
                            IF (kernel_type .EQ. 'v') THEN
                              acc_z = acc_z &
                                - KERNEL_LOWER_C0(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)*peso &
                                + KERNEL_LOWER_F1(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)*peso
                            ELSE
                              acc_z = acc_z &
                                + KERNEL_LOWER_ESC(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)*peso
                            END IF
                          END DO
                        END DO
                      END DO
                    END IF

!===================================================================
!                   LADO ESQUERDO: norma, so integral em gamma'
!===================================================================
                    CALL SJX_H(zv, NIZ, z, l-1, SPLz_col)
                    acc_x = acc_x + SPLg*SPLz_col*dgp &
                          / ((g + gp + (1.d0-z**2)*kappa**2 + m**2*z**2)**2)

                  END DO
                END DO

                ZMATRIX(index1, index2) = acc_z
                XMATRIX(index1, index2) = acc_x

              END DO
            END DO

!$OMP ATOMIC CAPTURE
            NFEITOS = NFEITOS + 1
            NFEITOS_LOC = NFEITOS
!$OMP END ATOMIC

            IF (MOD(NFEITOS_LOC,NPASSO).EQ.0 .OR. NFEITOS_LOC.EQ.NTOTAL) THEN
              TAGORA = 0.d0
!$            TAGORA = OMP_GET_WTIME()
              TDECOR = TAGORA - TSTART
              PCT    = 100.d0*DBLE(NFEITOS_LOC)/DBLE(NTOTAL)
              TRESTA = 0.d0
              IF (NFEITOS_LOC.GT.0) &
                 TRESTA = TDECOR*DBLE(NTOTAL-NFEITOS_LOC)/DBLE(NFEITOS_LOC)
              WRITE(*,'(A,I6,A,I6,A,F6.2,A,F9.1,A,F9.1,A)') &
                 "   progresso: ", NFEITOS_LOC, " / ", NTOTAL, &
                 "  (", PCT, "%)   decorrido ", TDECOR, &
                 " s   resta ~", TRESTA, " s"
              FLUSH(6)
            END IF

          END DO
        END DO
!$OMP END PARALLEL DO

        TEND = 0.d0
!$      TEND = OMP_GET_WTIME()
!$      WRITE(*,'(A,F12.3,A)') " Tempo de montagem das matrizes = ", &
!$                             TEND-TSTART, " s"

      END SUBROUTINE BUILD_MATRICES


!=======================================================================
!  SPLIT_RANGE
!  Monta a lista de pontos de quebra para integrar de A ate B sobre a
!  malha G: os extremos A e B mais todos os nos internos estritamente
!  entre eles. Devolve BRK(0:NSUB) com NSUB subintervalos.
!
!  E isto que impede a quadratura de Gauss de atravessar um no (onde a
!  base de Hermite tem derivada terceira descontinua) ou o ponto z'=z
!  (onde o kernel troca de ramo). Sem isso a taxa de convergencia em z
!  cai de espectral para algebrica.
!=======================================================================
      SUBROUTINE SPLIT_RANGE(G, N, A, B, BRK, NSUB)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: N
      DOUBLE PRECISION, INTENT(IN)  :: G(0:N), A, B
      DOUBLE PRECISION, INTENT(OUT) :: BRK(0:N+2)
      INTEGER, INTENT(OUT) :: NSUB
      INTEGER :: i
      DOUBLE PRECISION, PARAMETER :: EPS = 1.D-14

      NSUB = 0
      BRK(0) = A
      DO i = 0, N
        IF (G(i) .GT. A + EPS .AND. G(i) .LT. B - EPS) THEN
          NSUB = NSUB + 1
          BRK(NSUB) = G(i)
        END IF
      END DO
      NSUB = NSUB + 1
      BRK(NSUB) = B

      END SUBROUTINE SPLIT_RANGE


! =======================================================================
!                    DENOMINADORES COMUNS DO KERNEL
!     (identicos aos da versao anterior -- a fisica nao mudou)
! =======================================================================

      DOUBLE PRECISION FUNCTION DENOM_D0(z, g, m1, m2, Mtot)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, g, m1, m2, Mtot

      DENOM_D0 = 0.25d0*(4.d0*g + Mtot**2*(z**2 - 1.d0) &
                 + 2.d0*m1**2*(z + 1.d0) - 2.d0*m2**2*(z - 1.d0))
      END FUNCTION DENOM_D0


      DOUBLE PRECISION FUNCTION DENOM_DU(z, zq, g, gp, v, m1, m2, mu, Mtot)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot

      DENOM_DU = 0.25d0 * ( &
            Mtot**2 * (-v) * (z + 1.0d0) * (zq + 1.0d0) * ((v - 1.0d0) * z - v * zq + 1.0d0) &
            + v * ( m1**2 * (z + 1.0d0) + 2.0d0 * m1 * m2 * (z + 1.0d0) + &
                    m2**2 * (4.0d0 * v * z - 4.0d0 * v * zq - 3.0d0 * z + 4.0d0 * zq + 1.0d0) - &
                    4.0d0 * g * (v - 1.0d0) * (zq + 1.0d0) + 4.0d0 * gp * (z + 1.0d0) ) &
            - 4.0d0 * mu**2 * (v - 1.0d0) * (z + 1.0d0) )
      END FUNCTION DENOM_DU


      DOUBLE PRECISION FUNCTION DENOM_DD(z, zq, g, gp, v, m1, m2, mu, Mtot)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot

      DENOM_DD = 0.25d0 * ( &
          Mtot**2 * v * (z - 1.0d0) * (zq - 1.0d0) * ((v - 1.0d0) * z - v * zq - 1.0d0) &
          + v * ( m1**2 * (-4.0d0 * v * z + 4.0d0 * (v - 1.0d0) * zq + 3.0d0 * z + 1.0d0) &
                  - 2.0d0 * m1 * m2 * (z - 1.0d0) - m2**2 * z + m2**2 &
                  + 4.0d0 * (gp + g * (v - 1.0d0) * (zq - 1.0d0)) - 4.0d0 * gp * z ) &
          + 4.0d0 * mu**2 * (v - 1.0d0) * (z - 1.0d0) )
      END FUNCTION DENOM_DD


! =======================================================================
!                    COEFICIENTES VETORIAIS C0
! =======================================================================

      DOUBLE PRECISION FUNCTION COEF_C0_U(z, zq, g, v, m2, Mtot)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, v, m2, Mtot

      COEF_C0_U = (-2 + v)*(1 + z)*(m2**2*((-2 + v)*z - v*zq) - g*(-2 + v + v*zq)) - &
                    (Mtot**2*(1 + z)**2*(4 + ((-2 + v)*z - v*zq)*(-2 + v + v*zq)))/4.
      END FUNCTION COEF_C0_U


      DOUBLE PRECISION FUNCTION COEF_C0_D(z, zq, g, v, m1, Mtot)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, v, m1, Mtot

      COEF_C0_D = (-2 + v)*(-1 + z)*(g*(-2 + v - v*zq) + m1**2*((-2 + v)*z - v*zq)) - &
            (Mtot**2*(-1 + z)**2*(4 + (2 + v*(-1 + zq))*((-2 + v)*z - v*zq)))/4.
      END FUNCTION COEF_C0_D


! =======================================================================
!                        KERNEL ESCALAR
! =======================================================================

      DOUBLE PRECISION FUNCTION KERNEL_UPPER_ESC(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Du
      DOUBLE PRECISION, EXTERNAL :: DENOM_D0, DENOM_DU

      D0 = DENOM_D0(z, g, m1, m2, Mtot)
      Du = DENOM_DU(z, zq, g, gp, v, m1, m2, mu, Mtot)

      KERNEL_UPPER_ESC = 1.d0 / (32*PI**2*D0) * (v**2 / (Du**2)) * (1.d0 + z)**2
      END FUNCTION KERNEL_UPPER_ESC


      DOUBLE PRECISION FUNCTION KERNEL_LOWER_ESC(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Dd
      DOUBLE PRECISION, EXTERNAL :: DENOM_D0, DENOM_DD

      D0 = DENOM_D0(z, g, m1, m2, Mtot)
      Dd = DENOM_DD(z, zq, g, gp, v, m1, m2, mu, Mtot)

      KERNEL_LOWER_ESC = 1.d0 / (32*PI**2*D0) * (v**2 / (Dd**2)) * (1.d0 - z)**2
      END FUNCTION KERNEL_LOWER_ESC


! =======================================================================
!                        KERNEL VETORIAL
! =======================================================================

      DOUBLE PRECISION FUNCTION KERNEL_UPPER_C0(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Du, Co_ku
      DOUBLE PRECISION, EXTERNAL :: DENOM_D0, DENOM_DU, COEF_C0_U

      D0    = DENOM_D0(z, g, m1, m2, Mtot)
      Du    = DENOM_DU(z, zq, g, gp, v, m1, m2, mu, Mtot)
      Co_ku = COEF_C0_U(z, zq, g, v, m2, Mtot)

      KERNEL_UPPER_C0 = 1.d0 / (32*PI**2*D0) * (v**2 / (Du**2)) * Co_ku
      END FUNCTION KERNEL_UPPER_C0


      DOUBLE PRECISION FUNCTION KERNEL_UPPER_F1(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Du
      DOUBLE PRECISION, EXTERNAL :: DENOM_D0, DENOM_DU

      D0 = DENOM_D0(z, g, m1, m2, Mtot)
      Du = DENOM_DU(z, zq, g, gp, v, m1, m2, mu, Mtot)

      KERNEL_UPPER_F1 = 2.d0 / (32*PI**2*D0) * (v**2 / Du) * (z + 1.d0)
      END FUNCTION KERNEL_UPPER_F1


      DOUBLE PRECISION FUNCTION KERNEL_LOWER_C0(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Dd, Co_kd
      DOUBLE PRECISION, EXTERNAL :: DENOM_D0, DENOM_DD, COEF_C0_D

      D0    = DENOM_D0(z, g, m1, m2, Mtot)
      Dd    = DENOM_DD(z, zq, g, gp, v, m1, m2, mu, Mtot)
      Co_kd = COEF_C0_D(z, zq, g, v, m1, Mtot)

      KERNEL_LOWER_C0 = 1.d0 / (32*PI**2*D0) * (v**2 / (Dd**2)) * Co_kd
      END FUNCTION KERNEL_LOWER_C0


      DOUBLE PRECISION FUNCTION KERNEL_LOWER_F1(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Dd
      DOUBLE PRECISION, EXTERNAL :: DENOM_D0, DENOM_DD

      D0 = DENOM_D0(z, g, m1, m2, Mtot)
      Dd = DENOM_DD(z, zq, g, gp, v, m1, m2, mu, Mtot)

      KERNEL_LOWER_F1 = 2.d0 / (32*PI**2*D0) * (v**2 / Dd) * (1.d0 - z)
      END FUNCTION KERNEL_LOWER_F1


! =======================================================================
!                   QUADRATURA E MALHAS (inalteradas)
! =======================================================================

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
