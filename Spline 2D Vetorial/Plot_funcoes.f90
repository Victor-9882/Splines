    IMPLICIT NONE

    ! -----------------------------------------------------------------
    ! BLOCOS COMMON
    !   PI       : constante pi, compartilhada com as rotinas
    !   ALPHAINT : malhas e pesos de quadratura (Gauss-Legendre) usados
    !              tanto aqui quanto pelas rotinas legauss/G1D/COLLOC
    ! -----------------------------------------------------------------
    DOUBLE PRECISION :: PI
    COMMON/PARAM/PI

    DOUBLE PRECISION :: X, DX, Y, DY, W, DW
    INTEGER          :: Nz, Ng, Nv
    COMMON/ALPHAINT/X(1000), DX(1000), Y(1000), DY(1000) &
        , W(1000), DW(1000), Nz, Ng, Nv

    ! -----------------------------------------------------------------
    ! DIMENSOES DO PROBLEMA
    !   NMZ/NMG : numero de splines em z e em gamma
    !   NMA     : dimensao total do produto tensorial (NMG*NMZ)
    !   NPARAM  : quantos pares (Nng,Nnz) serao lidos de inputs.dat
    ! -----------------------------------------------------------------
    INTEGER :: NMG, NMZ, NMA, NPARAM
    INTEGER :: N_intervalZ, N_intervalG
    INTEGER :: Nnz(100), Nng(100)

    ! -----------------------------------------------------------------
    ! MALHAS E BASES DE SPLINE
    !   zv/gv     : nos das malhas em z e em gamma
    !   splz/splg : valores das bases de spline no ponto avaliado
    !   XG/YG     : pontos de colocacao de Gauss em z e em gamma
    !   c         : coeficientes da expansao, lidos de coeficientes.dat
    ! -----------------------------------------------------------------
    DOUBLE PRECISION, ALLOCATABLE :: zv(:), gv(:)
    DOUBLE PRECISION, ALLOCATABLE :: splz(:), splg(:)
    DOUBLE PRECISION, ALLOCATABLE :: XG(:), YG(:)
    DOUBLE PRECISION, ALLOCATABLE :: c(:,:)

    ! -----------------------------------------------------------------
    ! PARAMETROS FISICOS
    !   m1, m2 : massas dos constituintes; m e a media
    !   Mtot   : massa total do sistema ligado
    !   kappa  : sqrt(m^2 - Mtot^2/4), momento de ligacao
    ! -----------------------------------------------------------------
    DOUBLE PRECISION :: m1, m2, m, Mtot, kappa

    ! -----------------------------------------------------------------
    ! CONTROLE DOS PLOTS
    !   iw                     : unidade de saida das rotinas de malha
    !   N_PLOT/N_PLOT_G/_Z     : numero de pontos de cada varredura
    !   max_gamma_visualizacao : truncamento superior em gamma
    !   z_fixo/gamma_fixo      : cortes usados nos plots de g
    !   gamma_num/gamma_den    : gamma do numerador/denominador de Psi
    !   z_num/z_den            : z do numerador/denominador de Psi
    ! -----------------------------------------------------------------
    INTEGER          :: iw, N_PLOT, N_PLOT_G, N_PLOT_Z
    DOUBLE PRECISION :: max_gamma_visualizacao
    DOUBLE PRECISION :: z_fixo, gamma_fixo
    DOUBLE PRECISION :: gamma_num, gamma_den, z_num, z_den

    ! -----------------------------------------------------------------
    ! CONTADORES DE LOOP
    !   i, j       : indices das bases splg(i)/splz(j)
    !   ii         : varre os pares (Nng,Nnz) de inputs.dat
    !   p          : quadratura interna em gamma'
    !   q, r       : quadraturas externas em gamma e em xi
    !   k_plot, p_plot : indices das varreduras dos plots
    ! -----------------------------------------------------------------
    INTEGER :: i, j, ii, p, q, r, k_plot, p_plot

    ! -----------------------------------------------------------------
    ! ACUMULADORES E VARIAVEIS DE TRABALHO
    !   g1_00      : valor de referencia de g, usado para normalizar
    !   soma       : reconstrucao de g(gamma,z) a partir das splines
    !   g_val      : g(gamma',z) dentro da integral do propagador
    !   D_num/D_den: denominador do propagador ao quadrado
    !   psi_num/den: Psi no ponto do plot e no ponto de normalizacao
    !   psi_norm   : razao Psi/Psi_den
    !   gammap/dgp : ponto e peso da quadratura interna em gamma'
    !   *_ext      : ponto e peso das quadraturas externas
    !   Norma      : pi*Int dxi Int dgamma |Psi|^2, para reescalar
    !   dist_long/trans : marginais u(xi) e D_perp(gamma)
    ! -----------------------------------------------------------------
    DOUBLE PRECISION :: g1_00, soma, g_val
    DOUBLE PRECISION :: D_num, D_den, psi_num, psi_den, psi_norm
    DOUBLE PRECISION :: gammap, dgp
    DOUBLE PRECISION :: gamma_plot, z_plot, xi_plot
    DOUBLE PRECISION :: gamma_ext, dgamma_ext, xi_ext, dxi_ext, z_ext
    DOUBLE PRECISION :: Norma, dist_long, dist_trans


    open (unit = 10, file = "autovalores.dat",STATUS="UNKNOWN")
    open (unit = 11, file = "autovetoresG.dat",STATUS="UNKNOWN")
    open (unit = 15, file = "autovetoresZ.dat",STATUS="UNKNOWN")
    open (unit = 16, file = "coeficientes.dat",STATUS="UNKNOWN")
    open (unit = 12, file = "alfa.dat",STATUS="UNKNOWN")
    open (unit = 13, file = "plotalfa.dat",STATUS="UNKNOWN")
    open (unit = 14, file = 'erros.dat', status='unknown')
    open (UNIT = 20, FILE = "inputs.dat", STATUS="UNKNOWN")

    PI = DACOS(-1.D0)       !3.14159265358979323846264338

        !Parâmetros
        !Massas
        Mtot = 2.3d0
        m1 = 1.0d0
        m2 = 2.3d0
        m = (m1+m2)/2
        kappa = sqrt(m**2 - 0.25*Mtot**2)

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

        !ALOCAR VARIÁVEIS
        ALLOCATE( zv(NMZ + 1), gv(NMG+1) )
        ALLOCATE( splz(NMZ), splg(NMG))
        ALLOCATE( c(NMA, NMA) )
        ALLOCATE( XG(NMA), YG(NMA) )

    ! Ler matriz de coeficientes
      DO I = 1, NMG
         READ(16, *) (c(I,J), J=1, NMZ)
      END DO

      ! 3. Feche o arquivo
      CLOSE(16)

        iw = 14
        N_intervalZ = (NMZ-1)/2
        N_intervalG = (NMG-1)/2

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


        g1_00 = 0.d0
        call SPLMD1(zv, Nmz, 0.d0, splz)
        call SPLMD2(gv, Nmg, 0.d0, splg)

        do j = 1, Nmz
            do i = 1, Nmg
                g1_00 = g1_00 + c(i,j) * splg(i) * splz(j)
            end do
        end do

       !g1_00 = abs(g1_00)
        g1_00 = 1.d0
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
                    write(11, '(ES25.17E3, 1X, ES25.17E3)') gamma_plot, soma/g1_00
                else
                    ! Imprime normalmente (o descritor ES já deixa um espaço natural para positivos)
                    write(11, '(2ES25.17E3)') gamma_plot, soma/g1_00
                end if
          end do

      ! Plot da função variando z para gamma fixo
      gamma_fixo = 2.d0
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
              write(15, '(ES25.17E3, 1X, ES25.17E3)') z_plot, soma/g1_00
          else
              ! Imprime normalmente
              write(15, '(2ES25.17E3)') z_plot, soma/g1_00
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
          gammap = Y(p)
          dgp    = dY(p)

          ! Avalia a spline em gama'
          call SPLMD2(gv, Nmg, gammap, splg)

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
              gammap = Y(p)
              dgp    = dY(p)

              call SPLMD2(gv, Nmg, gammap, splg)

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
      gamma_num = 2.0d0
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
              gammap = Y(p)
              dgp    = dY(p)

              ! Avalia a spline em gamma'
              call SPLMD2(gv, Nmg, gammap, splg)
              
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
          ! Espaços fixos (1X) evitam que um sinal negativo cole no número
          ! anterior, o que impede o xmgrace de ler a coluna corretamente
          write(17, '(1X, ES25.17E3, 1X, ES25.17E3)') z_plot, psi_norm
      end do
      
      close(17)

      ! =============================================================
      ! DISTRIBUICOES DE MOMENTO  (eqs. 4, 5 e 6)
      !
      ! Psi(gamma,xi) = Int_0^inf dgamma' g(gamma',z)
      !                 / [gamma + gamma' + m^2 z^2 + (1-z^2) kappa^2]^2
      ! com a mudanca de variavel  xi = (1-z)/2  <=>  z = 1 - 2*xi,
      ! de modo que xi in [0,1] corresponde a z in [-1,1] e dz = -2 dxi.
      !
      ! (4) f1(gamma,xi)  = pi * |Psi(gamma,xi)|^2          -> plot 3D
      ! (5) u(xi)         = pi * Int_0^inf dgamma |Psi|^2   -> PDF (long.)
      ! (6) D_perp(gamma) = pi * Int_0^1   dxi   |Psi|^2    -> transversal
      !
      ! NORMALIZACAO (eq. da mensagem):
      !     pi * Int_0^1 dxi Int_0^inf dgamma |Psi(xi,gamma)|^2 = 1
      ! Impomos isso reescalando Psi por 1/sqrt(Norma), onde
      !     Norma = pi * Int_0^1 dxi Int_0^{gamma_max} dgamma |Psi|^2
      ! calculada com a MESMA Psi usada nos tres plots. Assim f1 integra
      ! 1 no plano (gamma,xi) e u(xi), D_perp(gamma) sao marginais
      ! consistentes: Int_0^1 u dxi = Int_0^inf D_perp dgamma = 1.
      ! =============================================================

      ! Pontos/pesos de Gauss em gamma' para a integral interna do propagador
      CALL legauss(0.d0, 3.d0, Ng, Y, dY, 1.d-15)

      ! Malha dos plots
      N_PLOT_G = 200
      N_PLOT_Z = 200
      max_gamma_visualizacao = 3.d0

      ! ------------------------------------------------------------
      ! CONSTANTE DE NORMALIZACAO
      !   Norma = pi * Int_0^1 dxi Int_0^{gamma_max} dgamma |Psi|^2
      ! Quadratura de Gauss nas duas variaveis: xi direto em [0,1]
      ! (o jacobiano dz = -2 dxi ja esta embutido em integrar em xi).
      ! ------------------------------------------------------------
      CALL legauss(0.d0, 3.d0, Nv, W, dW, 1.d-15)                  ! em gamma
      CALL legauss(0.d0, 1.d0, Nz, X, dX, 1.d-15)                  ! em xi

      Norma = 0.d0
      do r = 1, Nz
          xi_ext  = X(r)
          dxi_ext = dX(r)

          z_ext = 1.d0 - 2.d0*xi_ext
          if (z_ext .le. -1.d0) z_ext = -0.999999d0
          if (z_ext .ge.  1.d0) z_ext =  0.999999d0

          call SPLMD1(zv, Nmz, z_ext, splz)

          do q = 1, Nv
              gamma_ext  = W(q)
              dgamma_ext = dW(q)

              psi_num = 0.d0
              do p = 1, Ng
                  gammap = Y(p)
                  dgp    = dY(p)

                  call SPLMD2(gv, Nmg, gammap, splg)

                  g_val = 0.d0
                  do j = 1, Nmz
                      do i = 1, Nmg
                          g_val = g_val + c(i,j) * splg(i) * splz(j)
                      end do
                  end do

                  D_num = gamma_ext + gammap + (m**2)*(z_ext**2) &
                          + (1.d0 - z_ext**2)*(kappa**2)

                  psi_num = psi_num + (g_val * dgp) / (D_num**2)
              end do

              psi_norm = psi_num / psi_den

              Norma = Norma + PI * (psi_norm**2) * dgamma_ext * dxi_ext
          end do
      end do

      WRITE(*,'(A,ES25.17E3)') &
          " Norma = pi*Int dxi Int dgamma |Psi|^2 (antes do reescalonamento) = ", Norma

      ! ------------------------------------------------------------
      ! (a) EQ. 4 -> PLOT 3D:  f1(gamma,xi) = pi |Psi(gamma,xi)|^2
      !     Colunas: gamma, xi, f1
      ! ------------------------------------------------------------
      open(unit = 21, file = "plot_psi2_3d.dat", STATUS="UNKNOWN")

      do k_plot = 0, N_PLOT_Z
          xi_plot = dble(k_plot) / dble(N_PLOT_Z)

          z_plot = 1.d0 - 2.d0*xi_plot
          if (z_plot .le. -1.d0) z_plot = -0.999999d0
          if (z_plot .ge.  1.d0) z_plot =  0.999999d0

          call SPLMD1(zv, Nmz, z_plot, splz)

          do p_plot = 0, N_PLOT_G
              gamma_plot = (dble(p_plot) / dble(N_PLOT_G)) * max_gamma_visualizacao

              psi_num = 0.d0
              do p = 1, Ng
                  gammap = Y(p)
                  dgp    = dY(p)

                  call SPLMD2(gv, Nmg, gammap, splg)

                  g_val = 0.d0
                  do j = 1, Nmz
                      do i = 1, Nmg
                          g_val = g_val + c(i,j) * splg(i) * splz(j)
                      end do
                  end do

                  D_num = gamma_plot + gammap + (m**2)*(z_plot**2) &
                          + (1.d0 - z_plot**2)*(kappa**2)

                  psi_num = psi_num + (g_val * dgp) / (D_num**2)
              end do

              psi_norm = psi_num / psi_den

              ! f1 = pi |Psi|^2, com |Psi|^2 ja reescalado por 1/Norma
              ! (linhas em branco separam blocos p/ gnuplot)
              write(21, '(3(1X,ES25.17E3))') gamma_plot, xi_plot, &
                  PI * (psi_norm**2) / Norma
          end do
          write(21, *) ""
      end do

      close(21)

      ! ------------------------------------------------------------
      ! (b) EQ. 5 -> PDF LONGITUDINAL:
      !     u(xi) = pi * Int_0^inf dgamma |Psi(gamma,xi)|^2
      !     Normalizada: Int_0^1 u(xi) dxi = 1
      !     Colunas: xi, u(xi)
      ! ------------------------------------------------------------
      open(unit = 22, file = "plot_dist_long.dat", STATUS="UNKNOWN")

      ! Pontos de Gauss para a integral EXTERNA em gamma
      CALL legauss(0.d0, 3.d0, Nv, W, dW, 1.d-15)

      N_PLOT = 400
      do k_plot = 0, N_PLOT
          xi_plot = dble(k_plot) / dble(N_PLOT)

          z_plot = 1.d0 - 2.d0*xi_plot
          if (z_plot .le. -1.d0) z_plot = -0.999999d0
          if (z_plot .ge.  1.d0) z_plot =  0.999999d0

          call SPLMD1(zv, Nmz, z_plot, splz)

          dist_long = 0.d0
          do q = 1, Nv
              gamma_ext = W(q)
              dgamma_ext = dW(q)

              psi_num = 0.d0
              do p = 1, Ng
                  gammap = Y(p)
                  dgp    = dY(p)

                  call SPLMD2(gv, Nmg, gammap, splg)

                  g_val = 0.d0
                  do j = 1, Nmz
                      do i = 1, Nmg
                          g_val = g_val + c(i,j) * splg(i) * splz(j)
                      end do
                  end do

                  D_num = gamma_ext + gammap + (m**2)*(z_plot**2) &
                          + (1.d0 - z_plot**2)*(kappa**2)

                  psi_num = psi_num + (g_val * dgp) / (D_num**2)
              end do

              psi_norm = psi_num / psi_den

              dist_long = dist_long + PI * (psi_norm**2) * dgamma_ext
          end do

          write(22, '(2(1X,ES25.17E3))') xi_plot, dist_long/Norma
      end do

      close(22)

      ! ------------------------------------------------------------
      ! (c) EQ. 6 -> DISTRIBUICAO TRANSVERSAL:
      !     D_perp(gamma) = pi * Int_0^1 dxi |Psi(gamma,xi)|^2
      !     Normalizada: Int_0^inf D_perp(gamma) dgamma = 1
      !     Colunas: gamma, D_perp(gamma)
      ! ------------------------------------------------------------
      open(unit = 23, file = "plot_dist_trans.dat", STATUS="UNKNOWN")

      ! Pontos de Gauss para a integral EXTERNA em xi, direto em [0,1]
      CALL legauss(0.d0, 1.d0, Nz, X, dX, 1.d-15)

      N_PLOT = 800
      do k_plot = 0, N_PLOT
          ! Varredura uniforme em gamma ate o truncamento gamma_max
          gamma_plot = (dble(k_plot) / dble(N_PLOT)) * max_gamma_visualizacao

          dist_trans = 0.d0
          do q = 1, Nz
              xi_ext  = X(q)
              dxi_ext = dX(q)

              z_ext = 1.d0 - 2.d0*xi_ext
              if (z_ext .le. -1.d0) z_ext = -0.999999d0
              if (z_ext .ge.  1.d0) z_ext =  0.999999d0

              call SPLMD1(zv, Nmz, z_ext, splz)

              psi_num = 0.d0
              do p = 1, Ng
                  gammap = Y(p)
                  dgp    = dY(p)

                  call SPLMD2(gv, Nmg, gammap, splg)

                  g_val = 0.d0
                  do j = 1, Nmz
                      do i = 1, Nmg
                          g_val = g_val + c(i,j) * splg(i) * splz(j)
                      end do
                  end do

                  D_num = gamma_plot + gammap + (m**2)*(z_ext**2) &
                          + (1.d0 - z_ext**2)*(kappa**2)

                  psi_num = psi_num + (g_val * dgp) / (D_num**2)
              end do

              psi_norm = psi_num / psi_den

              dist_trans = dist_trans + PI * (psi_norm**2) * dxi_ext
          end do

          write(23, '(2(1X,ES25.17E3))') gamma_plot, dist_trans/Norma
      end do

      close(23)


      DEALLOCATE(splz, splg, c, XG, YG, gv, zv)


      end do

      close (15)
      CLOSE(10)
      CLOSE(12)
      close(13)
      close (14)
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











        
