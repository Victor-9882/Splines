!=======================================================================
!  PLOT_FUNCOES_CONV
!
!  Versao de Plot_funcoes.f90 adaptada a
!  2D_Spline_Vetorial_Diferentes_omp_conv.f90.
!
!  O QUE MUDOU EM RELACAO A Plot_funcoes.f90
!
!  (1) BASE DE HERMITE NO LUGAR DA SPLINE NATURAL GLOBAL
!      SPLGR1/SPLMD1 e SPLGR2/SPLMD2 construiam tabelas em COMMON e
!      devolviam o vetor inteiro das bases de uma vez. A base agora e a
!      de Hermite local SJX_H, copiada verbatim do programa principal:
!      uma chamada por funcao de base, indice J de base 0.
!      Aqui isso vira EVAL_BASIS, que preenche o vetor spl(1:NM) com
!      spl(J+1) = SJX_H(...,J,...), mantendo o resto do codigo com a
!      mesma cara de antes.
!
!  (2) MALHAS COM INDICE DE 0 A N, E SEM OS ARTIFICIOS DAS BORDAS
!      A spline natural impunha segunda derivada nula nos extremos, o
!      que obrigava a afastar os nos das bordas:
!          zv(1)=-0.999999, zv(nmz)=0.999999, gv(1)=1.d-7
!      A base de Hermite nao impoe nada nos extremos, entao os nos sao
!      exatamente -1, +1 e 0, GMAX. Os limites tambem sumiram dos plots:
!      z=+-1 pode ser avaliado direto.
!
!  (3) DIMENSOES VEM DOS PARAMETERS, NAO DE inputs_plot.dat
!      No programa principal NIZ e NIG sao PARAMETER de tempo de
!      compilacao e fixam NMZ=2*(NIZ+1), NMG=2*(NIG+1). Nao ha como
!      escolher isso livremente aqui: os coeficientes so fazem sentido
!      sobre a MESMA malha que os gerou. Entao NIZ/NIG estao repetidos
!      abaixo e PRECISAM casar com os do programa principal.
!      O programa confere o tamanho real de coeficientes_conv.dat e
!      PARA com mensagem clara se nao bater, em vez de ler lixo.
!      (inputs_plot.dat deixou de ser lido; nao era mais compativel.)
!
!  (4) GAMMA TRUNCADO EM GMAX, COERENTE COM O PRINCIPAL
!      GMAX e um unico PARAMETER, usado na malha, nas quadraturas em
!      gamma' e no alcance dos plots. Antes o 3.d0 aparecia solto em
!      seis lugares diferentes.
!      Vale a mesma ressalva do cabecalho do principal: a cauda
!      gamma > GMAX e descartada, entao as distribuicoes normalizadas
!      dependem de GMAX.
!
!  (5) LE coeficientes_conv.dat / ESCREVE *_conv.dat
!      Nenhum arquivo da versao antiga e lido ou sobrescrito.
!
!  O QUE NAO MUDOU
!      A fisica: Psi, o propagador de Nakanishi, a troca xi=(1-z)/2 e
!      as eqs. (4), (5) e (6) estao como estavam.
!
!  Compilacao:
!    gfortran -O2 -o plot_conv Plot_funcoes_conv.f90
!=======================================================================
      PROGRAM PLOT_FUNCOES_CONV
      IMPLICIT NONE

!===================== Discretizacao ===================================
!     ATENCAO: estes dois PRECISAM ser iguais aos de
!     2D_Spline_Vetorial_Diferentes_omp_conv.f90. Se forem alterados la,
!     altere aqui tambem -- o programa avisa se nao casarem.
      INTEGER, PARAMETER :: NIZ = 15                !Intervalos em z
      INTEGER, PARAMETER :: NIG = 15                !Intervalos em gamma
      INTEGER, PARAMETER :: NMZ = 2*(NIZ+1)         !Funcoes de base em z
      INTEGER, PARAMETER :: NMG = 2*(NIG+1)         !Funcoes de base em gamma

!     Truncamento em gamma, igual ao do programa principal.
      DOUBLE PRECISION, PARAMETER :: GMAX = 3.d0

!     Menor gamma da malha logaritmica do plot transversal.
      DOUBLE PRECISION, PARAMETER :: GMIN_LOG = 1.d-4

!===================== Malhas ==========================================
!     Indexadas de 0 a N, como no programa principal.
      DOUBLE PRECISION :: zv(0:NIZ), gv(0:NIG)
      DOUBLE PRECISION :: splz(NMZ), splg(NMG)
      DOUBLE PRECISION :: c(NMG,NMZ)
      INTEGER :: IW

!===================== Quadratura ======================================
      DOUBLE PRECISION :: X(1000), DX(1000)         !em xi
      DOUBLE PRECISION :: Y(1000), DY(1000)         !em gamma' (interna)
      DOUBLE PRECISION :: W(1000), DW(1000)         !em gamma (externa)
      INTEGER          :: Nz, Ng, Nv

!===================== Parametros fisicos ==============================
      DOUBLE PRECISION :: m1, m2, m, Mtot, kappa
      DOUBLE PRECISION :: PI

!===================== Controle dos plots ==============================
      INTEGER          :: N_PLOT, N_PLOT_G, N_PLOT_Z
      DOUBLE PRECISION :: z_fixo, gamma_fixo
      DOUBLE PRECISION :: gamma_num, gamma_den, z_num, z_den

!===================== Contadores ======================================
      INTEGER :: i, j, p, q, r, k_plot, p_plot
      INTEGER :: NLIN, NCOLN, IOS

!===================== Acumuladores ====================================
      DOUBLE PRECISION :: g1_00, soma, g_val
      DOUBLE PRECISION :: D_num, D_den, psi_num, psi_den, psi_norm
      DOUBLE PRECISION :: gammap, dgp
      DOUBLE PRECISION :: gamma_plot, z_plot, xi_plot
      DOUBLE PRECISION :: gamma_ext, dgamma_ext, xi_ext, dxi_ext, z_ext
      DOUBLE PRECISION :: Norma, dist_long, dist_trans

      PI = DACOS(-1.D0)

!===================== Parametros fisicos ==============================
!     Mantidos iguais aos do programa principal.
      Mtot  = 2.3d0
      m1    = 1.0d0
      m2    = 2.3d0
      m     = (m1 + m2)/2
      kappa = sqrt(m**2 - 0.25*Mtot**2)

      Nz = 60
      Ng = 60
      Nv = 60

      OPEN (unit = 11, file = "autovetoresG_conv.dat", STATUS="UNKNOWN")
      OPEN (unit = 15, file = "autovetoresZ_conv.dat", STATUS="UNKNOWN")
      OPEN (unit = 14, file = "erros_plot_conv.dat",   STATUS="UNKNOWN")

!===================== Conferencia do arquivo de coeficientes ==========
!     Le a matriz so depois de confirmar que as dimensoes batem com
!     NIZ/NIG deste arquivo. Ler uma matriz de tamanho errado daria
!     numeros plausiveis e completamente falsos, entao vale parar.
      CALL CHECK_COEF("coeficientes_conv.dat", NLIN, NCOLN, IOS)

      IF (IOS .NE. 0) THEN
        WRITE(*,*) "ERRO: nao consegui abrir/ler coeficientes_conv.dat"
        WRITE(*,*) "Rode antes o 2D_Spline_Vetorial_Diferentes_omp_conv."
        STOP 1
      ENDIF

      IF (NLIN .NE. NMG .OR. NCOLN .NE. NMZ) THEN
        WRITE(*,*) "ERRO: coeficientes_conv.dat nao casa com NIZ/NIG."
        WRITE(*,'(A,I0,A,I0)') "  arquivo tem  linhas=", NLIN,   &
                               "  colunas=", NCOLN
        WRITE(*,'(A,I0,A,I0)') "  esperado     NMG   =", NMG,    &
                               "  NMZ    =", NMZ
        WRITE(*,'(A,I0,A,I0)') "  (deste arquivo: NIG=", NIG,    &
                               ", NIZ=", NIZ
        WRITE(*,*) " "
        WRITE(*,*) "Ajuste NIZ/NIG no topo deste arquivo para os mesmos"
        WRITE(*,*) "valores usados em 2D_Spline_Vetorial_..._conv.f90,"
        WRITE(*,*) "ou regere coeficientes_conv.dat."
        STOP 1
      ENDIF

      OPEN (unit = 16, file = "coeficientes_conv.dat", STATUS="OLD")
      DO i = 1, NMG
        READ(16, *) (c(i,j), j = 1, NMZ)
      END DO
      CLOSE(16)

      WRITE(*,'(A,I0,A,I0,A)') " Coeficientes lidos: ", NMG, " x ", NMZ, &
                               " (gamma x z)"

!===================== Malhas ==========================================
!     Identicas as do programa principal: uniformes, com os extremos
!     como nos exatos. Sem o recuo de 1.d-6 das bordas, que existia so
!     para fugir da condicao de contorno da spline natural.
      IW = 0
      CALL G1D(IW, -1.d0, NIZ, 1.0d0, 1.d0,  zv)
      CALL G1D(IW,  0.d0, NIG, 1.0d0, GMAX,  gv)

!===================== Normalizacao de referencia ======================
!     Mantido o comportamento do arquivo antigo: g1_00 fica em 1.
      CALL EVAL_BASIS(zv, NIZ, 0.d0, splz, NMZ)
      CALL EVAL_BASIS(gv, NIG, 0.d0, splg, NMG)

      g1_00 = 0.d0
      DO j = 1, NMZ
        DO i = 1, NMG
          g1_00 = g1_00 + c(i,j) * splg(i) * splz(j)
        END DO
      END DO
      WRITE(*,'(A,ES20.12)') " g(0,0) reconstruido = ", g1_00
      g1_00 = 1.d0

!=======================================================================
!  PLOT DE g x gamma, PARA z FIXO
!=======================================================================
      z_fixo = 0.8d0
      N_PLOT = 1000
      CALL EVAL_BASIS(zv, NIZ, z_fixo, splz, NMZ)

      DO p = 0, N_PLOT
        gamma_plot = (dble(p) / dble(N_PLOT)) * GMAX

        CALL EVAL_BASIS(gv, NIG, gamma_plot, splg, NMG)

        soma = 0.d0
        DO j = 1, NMZ
          DO i = 1, NMG
            soma = soma + c(i,j) * splg(i) * splz(j)
          END DO
        END DO

        WRITE(11, '(2(1X,ES25.17E3))') gamma_plot, soma/g1_00
      END DO

!=======================================================================
!  PLOT DE g x z, PARA gamma FIXO
!=======================================================================
      gamma_fixo = 2.d0
      N_PLOT = 1500
      CALL EVAL_BASIS(gv, NIG, gamma_fixo, splg, NMG)

      DO p = 0, N_PLOT
!       Distribuicao em cosseno: adensa pontos perto de z=+-1.
        z_plot = -cos( (dble(p) / dble(N_PLOT)) * PI )

        CALL EVAL_BASIS(zv, NIZ, z_plot, splz, NMZ)

        soma = 0.d0
        DO j = 1, NMZ
          DO i = 1, NMG
            soma = soma + c(i,j) * splg(i) * splz(j)
          END DO
        END DO

        WRITE(15, '(2(1X,ES25.17E3))') z_plot, soma/g1_00
      END DO

!=======================================================================
!  PSI NORMALIZADA: Psi(gamma, z_num) / Psi(0, z_den)
!
!  Psi(gamma,z) = Int_0^GMAX dgamma' g(gamma',z)
!                 / [gamma + gamma' + m^2 z^2 + (1-z^2) kappa^2]^2
!=======================================================================
      OPEN(unit = 17, file = "plot_psi_conv.dat", STATUS="UNKNOWN")

!     Quadratura em gamma', truncada em GMAX como no principal.
      CALL legauss(0.d0, GMAX, Ng, Y, DY, 1.d-15)

!---- denominador: Psi(0, z_den) ---------------------------------------
      gamma_den = 0.0d0
      z_den     = 0.5d0
      psi_den   = 0.d0

      z_den = 1.d0 - 2.d0*z_den            !xi -> z
      CALL EVAL_BASIS(zv, NIZ, z_den, splz, NMZ)

      DO p = 1, Ng
        gammap = Y(p)
        dgp    = DY(p)

        CALL EVAL_BASIS(gv, NIG, gammap, splg, NMG)

        g_val = 0.d0
        DO j = 1, NMZ
          DO i = 1, NMG
            g_val = g_val + c(i,j) * splg(i) * splz(j)
          END DO
        END DO

        D_den = gamma_den + gammap + (m**2)*(z_den**2) &
                + (1.d0 - z_den**2)*(kappa**2)

        psi_den = psi_den + (g_val * dgp) / (D_den**2)
      END DO

      IF (ABS(psi_den) .LT. 1.d-300) THEN
        WRITE(*,*) "ERRO: Psi de normalizacao deu zero."
        STOP 1
      ENDIF

!---- numerador: Psi(gamma, z_num) -------------------------------------
      z_num  = 0.2d0
      N_PLOT = 1000

      z_num = 1.d0 - 2.d0*z_num
      CALL EVAL_BASIS(zv, NIZ, z_num, splz, NMZ)

      DO k_plot = 0, N_PLOT
        gamma_plot = (dble(k_plot) / dble(N_PLOT)) * GMAX
        psi_num = 0.d0

        DO p = 1, Ng
          gammap = Y(p)
          dgp    = DY(p)

          CALL EVAL_BASIS(gv, NIG, gammap, splg, NMG)

          g_val = 0.d0
          DO j = 1, NMZ
            DO i = 1, NMG
              g_val = g_val + c(i,j) * splg(i) * splz(j)
            END DO
          END DO

          D_num = gamma_plot + gammap + (m**2)*(z_num**2) &
                  + (1.d0 - z_num**2)*(kappa**2)

          psi_num = psi_num + (g_val * dgp) / (D_num**2)
        END DO

        psi_norm = psi_num / psi_den

        WRITE(17, '(2(1X,ES25.17E3))') gamma_plot, psi_norm
      END DO

      CLOSE(17)

!=======================================================================
!  PSI NORMALIZADA x z, PARA gamma FIXO
!=======================================================================
      OPEN(unit = 18, file = "plot_psi_z_conv.dat", STATUS="UNKNOWN")

      gamma_num = 2.0d0
      N_PLOT    = 1500

      DO k_plot = 0, N_PLOT
        z_plot = -cos( (dble(k_plot) / dble(N_PLOT)) * PI )

        psi_num = 0.d0

        CALL EVAL_BASIS(zv, NIZ, z_plot, splz, NMZ)

        DO p = 1, Ng
          gammap = Y(p)
          dgp    = DY(p)

          CALL EVAL_BASIS(gv, NIG, gammap, splg, NMG)

          g_val = 0.d0
          DO j = 1, NMZ
            DO i = 1, NMG
              g_val = g_val + c(i,j) * splg(i) * splz(j)
            END DO
          END DO

          D_num = gamma_num + gammap + (m**2)*(z_plot**2) &
                  + (1.d0 - z_plot**2)*(kappa**2)

          psi_num = psi_num + (g_val * dgp) / (D_num**2)
        END DO

        psi_norm = psi_num / psi_den

        WRITE(18, '(2(1X,ES25.17E3))') z_plot, psi_norm
      END DO

      CLOSE(18)

!=======================================================================
!  DISTRIBUICOES DE MOMENTO  (eqs. 4, 5 e 6)
!
!  Com xi = (1-z)/2  <=>  z = 1 - 2*xi, de modo que xi in [0,1]
!  corresponde a z in [-1,1] e dz = -2 dxi.
!
!  (4) f1(gamma,xi)  = pi * |Psi(gamma,xi)|^2          -> plot 3D
!  (5) u(xi)         = pi * Int dgamma |Psi|^2         -> PDF (long.)
!  (6) D_perp(gamma) = pi * Int_0^1 dxi |Psi|^2        -> transversal
!
!  NORMALIZACAO:
!      pi * Int_0^1 dxi Int_0^GMAX dgamma |Psi|^2 = 1
!  imposta reescalando |Psi|^2 por 1/Norma. As integrais em gamma vao
!  ate GMAX, nao ate infinito -- ver a ressalva no cabecalho.
!=======================================================================
      N_PLOT_G = 200
      N_PLOT_Z = 200

!---- constante de normalizacao ----------------------------------------
      CALL legauss(0.d0, GMAX, Nv, W, DW, 1.d-15)   !em gamma
      CALL legauss(0.d0, 1.d0, Nz, X, DX, 1.d-15)   !em xi

      Norma = 0.d0
      DO r = 1, Nz
        xi_ext  = X(r)
        dxi_ext = DX(r)

!       z=+-1 e um no valido da base de Hermite; sem recuo artificial.
        z_ext = 1.d0 - 2.d0*xi_ext

        CALL EVAL_BASIS(zv, NIZ, z_ext, splz, NMZ)

        DO q = 1, Nv
          gamma_ext  = W(q)
          dgamma_ext = DW(q)

          psi_num = 0.d0
          DO p = 1, Ng
            gammap = Y(p)
            dgp    = DY(p)

            CALL EVAL_BASIS(gv, NIG, gammap, splg, NMG)

            g_val = 0.d0
            DO j = 1, NMZ
              DO i = 1, NMG
                g_val = g_val + c(i,j) * splg(i) * splz(j)
              END DO
            END DO

            D_num = gamma_ext + gammap + (m**2)*(z_ext**2) &
                    + (1.d0 - z_ext**2)*(kappa**2)

            psi_num = psi_num + (g_val * dgp) / (D_num**2)
          END DO

          psi_norm = psi_num / psi_den

          Norma = Norma + PI * (psi_norm**2) * dgamma_ext * dxi_ext
        END DO
      END DO

      WRITE(*,'(A,ES25.17E3)') &
        " Norma = pi*Int dxi Int dgamma |Psi|^2 (antes do reescalonamento) = ", Norma

      IF (Norma .LE. 0.d0) THEN
        WRITE(*,*) "ERRO: Norma nao positiva; nada a normalizar."
        STOP 1
      ENDIF

!---- (a) eq. 4 -> plot 3D:  f1(gamma,xi) = pi |Psi|^2 / Norma ---------
      OPEN(unit = 21, file = "plot_psi2_3d_conv.dat", STATUS="UNKNOWN")

      DO k_plot = 0, N_PLOT_Z
        xi_plot = dble(k_plot) / dble(N_PLOT_Z)
        z_plot  = 1.d0 - 2.d0*xi_plot

        CALL EVAL_BASIS(zv, NIZ, z_plot, splz, NMZ)

        DO p_plot = 0, N_PLOT_G
          gamma_plot = (dble(p_plot) / dble(N_PLOT_G)) * GMAX

          psi_num = 0.d0
          DO p = 1, Ng
            gammap = Y(p)
            dgp    = DY(p)

            CALL EVAL_BASIS(gv, NIG, gammap, splg, NMG)

            g_val = 0.d0
            DO j = 1, NMZ
              DO i = 1, NMG
                g_val = g_val + c(i,j) * splg(i) * splz(j)
              END DO
            END DO

            D_num = gamma_plot + gammap + (m**2)*(z_plot**2) &
                    + (1.d0 - z_plot**2)*(kappa**2)

            psi_num = psi_num + (g_val * dgp) / (D_num**2)
          END DO

          psi_norm = psi_num / psi_den

!         Linha em branco entre blocos, para o gnuplot fazer a malha.
          WRITE(21, '(3(1X,ES25.17E3))') gamma_plot, xi_plot, &
            PI * (psi_norm**2) / Norma
        END DO
        WRITE(21, *) ""
      END DO

      CLOSE(21)

!---- (b) eq. 5 -> PDF longitudinal u(xi) ------------------------------
      OPEN(unit = 22, file = "plot_dist_long_conv.dat", STATUS="UNKNOWN")

      CALL legauss(0.d0, GMAX, Nv, W, DW, 1.d-15)

      N_PLOT = 400
      DO k_plot = 0, N_PLOT
        xi_plot = dble(k_plot) / dble(N_PLOT)
        z_plot  = 1.d0 - 2.d0*xi_plot

        CALL EVAL_BASIS(zv, NIZ, z_plot, splz, NMZ)

        dist_long = 0.d0
        DO q = 1, Nv
          gamma_ext  = W(q)
          dgamma_ext = DW(q)

          psi_num = 0.d0
          DO p = 1, Ng
            gammap = Y(p)
            dgp    = DY(p)

            CALL EVAL_BASIS(gv, NIG, gammap, splg, NMG)

            g_val = 0.d0
            DO j = 1, NMZ
              DO i = 1, NMG
                g_val = g_val + c(i,j) * splg(i) * splz(j)
              END DO
            END DO

            D_num = gamma_ext + gammap + (m**2)*(z_plot**2) &
                    + (1.d0 - z_plot**2)*(kappa**2)

            psi_num = psi_num + (g_val * dgp) / (D_num**2)
          END DO

          psi_norm = psi_num / psi_den

          dist_long = dist_long + PI * (psi_norm**2) * dgamma_ext
        END DO

        WRITE(22, '(2(1X,ES25.17E3))') xi_plot, dist_long/Norma
      END DO

      CLOSE(22)

!---- (c) eq. 6 -> distribuicao transversal D_perp(gamma) --------------
      OPEN(unit = 23, file = "plot_dist_trans_conv.dat", STATUS="UNKNOWN")

      CALL legauss(0.d0, 1.d0, Nz, X, DX, 1.d-15)

!     Malha LOGARITMICA em gamma.
!     D_perp e fortemente picada na origem: cai ~7 ordens de grandeza
!     de gamma=0 ate GMAX, com ~99% da massa em gamma < 0.15. Numa
!     malha uniforme o pico fica com pouquissimos pontos, e quem
!     integrar o arquivo depois pelo trapezio erra dezenas de por
!     cento (medido: 1.02 e 1.30 em vez de 1). Amostrando em
!     log(gamma) o pico aparece resolvido.
!     A primeira linha e gamma=0; as demais vao de GMIN_LOG a GMAX.
      N_PLOT = 800
      DO k_plot = 0, N_PLOT
        IF (k_plot .EQ. 0) THEN
          gamma_plot = 0.d0
        ELSE
          gamma_plot = GMIN_LOG * (GMAX/GMIN_LOG) &
                       ** (dble(k_plot-1) / dble(N_PLOT-1))
        ENDIF

        dist_trans = 0.d0
        DO q = 1, Nz
          xi_ext  = X(q)
          dxi_ext = DX(q)

          z_ext = 1.d0 - 2.d0*xi_ext

          CALL EVAL_BASIS(zv, NIZ, z_ext, splz, NMZ)

          psi_num = 0.d0
          DO p = 1, Ng
            gammap = Y(p)
            dgp    = DY(p)

            CALL EVAL_BASIS(gv, NIG, gammap, splg, NMG)

            g_val = 0.d0
            DO j = 1, NMZ
              DO i = 1, NMG
                g_val = g_val + c(i,j) * splg(i) * splz(j)
              END DO
            END DO

            D_num = gamma_plot + gammap + (m**2)*(z_ext**2) &
                    + (1.d0 - z_ext**2)*(kappa**2)

            psi_num = psi_num + (g_val * dgp) / (D_num**2)
          END DO

          psi_norm = psi_num / psi_den

          dist_trans = dist_trans + PI * (psi_norm**2) * dxi_ext
        END DO

        WRITE(23, '(2(1X,ES25.17E3))') gamma_plot, dist_trans/Norma
      END DO

      CLOSE(23)

      CLOSE(11)
      CLOSE(15)
      CLOSE(14)

      WRITE(*,*) " "
      WRITE(*,*) "Arquivos gerados:"
      WRITE(*,*) "  autovetoresG_conv.dat      g x gamma  (z fixo)"
      WRITE(*,*) "  autovetoresZ_conv.dat      g x z      (gamma fixo)"
      WRITE(*,*) "  plot_psi_conv.dat          Psi x gamma"
      WRITE(*,*) "  plot_psi_z_conv.dat        Psi x z"
      WRITE(*,*) "  plot_psi2_3d_conv.dat      f1(gamma,xi)  (eq. 4)"
      WRITE(*,*) "  plot_dist_long_conv.dat    u(xi)         (eq. 5)"
      WRITE(*,*) "  plot_dist_trans_conv.dat   D_perp(gamma) (eq. 6)"

      END PROGRAM PLOT_FUNCOES_CONV


!=======================================================================
!  EVAL_BASIS
!
!  Preenche spl(1:NM) com o valor das NM funcoes da base de Hermite no
!  ponto X. Faz o papel que SPLMD1/SPLMD2 faziam na versao antiga, mas
!  sem tabela em COMMON: chama SJX_H uma vez por funcao.
!
!  Cuidado com o indice: SJX_H usa J de 0 a 2*N+1 (convencao francesa),
!  enquanto os coeficientes c(i,j) sao de base 1. A correspondencia e
!      spl(J+1) <-> funcao de base J
!  que e exatamente a usada em BUILD_MATRICES no programa principal
!  (chamadas com k-1 e l-1).
!
!  Como cada funcao e nula fora de dois intervalos, so ~4 entradas de
!  spl saem diferentes de zero. O laco completo custa pouco na escala
!  destes plots e mantem o codigo simples.
!=======================================================================
      SUBROUTINE EVAL_BASIS(G, N, X, SPL, NM)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, NM
      DOUBLE PRECISION, INTENT(IN)  :: G(0:N), X
      DOUBLE PRECISION, INTENT(OUT) :: SPL(NM)
      INTEGER :: J
      DOUBLE PRECISION :: S

      SPL = 0.d0
      DO J = 0, NM-1
        CALL SJX_H(G, N, X, J, S)
        SPL(J+1) = S
      END DO

      END SUBROUTINE EVAL_BASIS


!=======================================================================
!  CHECK_COEF
!
!  Descobre quantas linhas e quantas colunas o arquivo de coeficientes
!  realmente tem, antes de tentar le-lo na matriz.
!
!  Serve para nao ler silenciosamente uma matriz de tamanho errado: se
!  NIZ/NIG aqui nao baterem com os do programa principal, os numeros
!  lidos seriam plausiveis e completamente sem sentido. As colunas sao
!  contadas pela primeira linha, que e o formato fixo escrito la
!  ('9999ES16.8' numa linha por valor de gamma).
!=======================================================================
      SUBROUTINE CHECK_COEF(FNAME, NLIN, NCOLN, IOS)
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: FNAME
      INTEGER, INTENT(OUT) :: NLIN, NCOLN, IOS
      CHARACTER(LEN=200000) :: LINHA
      INTEGER :: IST, I, L
      LOGICAL :: DENTRO

      NLIN  = 0
      NCOLN = 0
      IOS   = 0

      OPEN(unit=77, file=FNAME, status='OLD', iostat=IST)
      IF (IST .NE. 0) THEN
        IOS = IST
        RETURN
      ENDIF

      DO
        READ(77, '(A)', iostat=IST) LINHA
        IF (IST .NE. 0) EXIT
        IF (LEN_TRIM(LINHA) .EQ. 0) CYCLE

        NLIN = NLIN + 1

!       Conta os campos da primeira linha nao vazia: cada transicao de
!       branco para nao-branco e um numero.
        IF (NLIN .EQ. 1) THEN
          DENTRO = .FALSE.
          L = LEN_TRIM(LINHA)
          DO I = 1, L
            IF (LINHA(I:I) .NE. ' ' .AND. LINHA(I:I) .NE. CHAR(9)) THEN
              IF (.NOT. DENTRO) THEN
                NCOLN  = NCOLN + 1
                DENTRO = .TRUE.
              ENDIF
            ELSE
              DENTRO = .FALSE.
            ENDIF
          END DO
        ENDIF
      END DO

      CLOSE(77)

      END SUBROUTINE CHECK_COEF


!=======================================================================
!  SJX_H — base de Hermite cubica com suporte local
!
!  COPIA VERBATIM de 2D_Spline_Vetorial_Diferentes_omp_conv.f90.
!  Precisa ser identica: e a mesma base que gerou os coeficientes.
!  Se for alterada la, copie de novo para ca.
!
!  Valor da funcao de base J (0 <= J <= 2*N+1) no ponto X, sobre a malha
!  G(0:N). Cada no i tem duas funcoes:
!     J = 2*i    -> vale 1 em x_i, derivada 0   (funcao de valor)
!     J = 2*i+1  -> vale 0 em x_i, derivada 1   (funcao de derivada)
!  Ambas sao nulas fora de [x_{i-1}, x_{i+1}].
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


!=====================================================================
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
