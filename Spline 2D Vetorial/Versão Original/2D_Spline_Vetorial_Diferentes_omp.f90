      PROGRAM SPLINE_VETORIAL
      USE OMP_LIB
      IMPLICIT NONE

!=======================================================================
! Todos os parametros vem do arquivo de entrada (padrao: inputs.dat,
! ou o primeiro argumento da linha de comando), em formato NAMELIST:
!
!   &saida  pasta = '...' /          -> pasta onde tudo sera salvo
!   &teste  ...  /                   -> um bloco por teste (N blocos)
!
! Cada bloco &teste herda os valores do bloco anterior; so o primeiro
! precisa ser completo. Para cada teste sao gravados, dentro da pasta,
! os arquivos <nome>_*.dat, alem de resumo.dat (uma linha por teste) e
! uma copia do arquivo de entrada.
!=======================================================================

!===================== Dimensoes do problema (lidas por teste) =========
      INTEGER :: NMG                                !Numero de splines em gamma
      INTEGER :: NMZ                                !Numero de splines em z
      INTEGER :: NMA                                !Dimensao do problema de autovalores (= NMG*NMZ)
      INTEGER :: LWORK                              !Tamanho do buffer de trabalho do LAPACK (= 10*NMA)
      INTEGER, PARAMETER :: NMAX_SPL   = 499        !Limite das rotinas SPLGR* (NP1=500)
      INTEGER, PARAMETER :: NMAX_GAUSS = 1000       !Tamanho dos vetores do COMMON /ALPHAINT/
      INTEGER, PARAMETER :: NMAX_TESTES = 500       !Numero maximo de blocos &teste

!===================== Malhas e bases de splines =======================
      DOUBLE PRECISION, ALLOCATABLE :: zv(:), gv(:) !Nos das malhas em z e em gamma (= pontos de colocacao)
      INTEGER :: IW                                 !Unidade de saida dos diagnosticos de malha
      INTEGER :: N_intervalZ, N_intervalG           !Numero de intervalos das malhas
      INTEGER :: NCOL                               !Pontos de colocacao por intervalo

!===================== Matrizes do problema generalizado ===============
      DOUBLE PRECISION, ALLOCATABLE :: XMATRIX(:,:)  !Lado esquerdo (norma)
      DOUBLE PRECISION, ALLOCATABLE :: ZMATRIX(:,:)  !Lado direito (kernel)
      DOUBLE PRECISION, ALLOCATABLE :: c(:,:)        !Coeficientes cij do autovetor fundamental
      DOUBLE PRECISION, ALLOCATABLE :: XG(:), YG(:)  !Pontos de colocacao gerados por COLLOC

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
      DOUBLE PRECISION :: e                         !Regularizacao da diagonal de XMATRIX
      DOUBLE PRECISION :: Alfa                      !Constante de acoplamento extraida de wr(1)

!===================== Parametros das malhas ===========================
      DOUBLE PRECISION :: A_z, A_g                  !Razao geometrica dos intervalos (G1D) em z e gamma
      DOUBLE PRECISION :: gamma_max                 !Limite superior em gamma (malha e quadratura)
      DOUBLE PRECISION :: z_borda                   !Nos extremos em z: zv(1)=-z_borda, zv(NMZ)=+z_borda

!===================== Selecao do kernel ===============================
      !'v' -> kernel vetorial (C0 + f1);  qualquer outro -> kernel escalar
      CHARACTER :: kernel_type

!===================== Entrada (NAMELIST) ==============================
      CHARACTER(LEN=40)  :: nome                    !Prefixo dos arquivos do teste (vazio = automatico)
      CHARACTER(LEN=200) :: pasta                   !Pasta de saida de toda a rodada
      LOGICAL :: DEBUG_TERMOS                       !Grava <nome>_termos.dat
      NAMELIST /saida/ pasta
      NAMELIST /teste/ nome, kernel_type, NMG, NMZ, Nz, Ng, Nv, &
                       Mtot, m1, m2, mu, e, A_z, A_g, gamma_max, z_borda, &
                       DEBUG_TERMOS
      !Valor que marca "nao definido" para checar se o primeiro bloco esta completo
      DOUBLE PRECISION, PARAMETER :: NAO_DEF = -1.D300

!===================== Controle da rodada ==============================
      CHARACTER(LEN=256) :: arq_entrada             !Arquivo de entrada
      CHARACTER(LEN=256) :: prefixo                 !<pasta>/<nome>_
      CHARACTER(LEN=40)  :: nomes(NMAX_TESTES)      !Nomes ja usados (detecta repeticao)
      CHARACTER(LEN=256) :: msg_io
      INTEGER :: NTESTES                            !Numero de blocos &teste no arquivo
      INTEGER :: ios
      LOGICAL :: existe
      INTEGER :: data_hora(8)
      INTEGER(KIND=8) :: T_INI, T_FIM, T_TAXA       !Cronometragem de cada teste
      DOUBLE PRECISION :: T_TESTE
      DOUBLE PRECISION :: Alfa_ant                  !Alfa do teste anterior (para alfas.dat)
      LOGICAL :: tem_ant                            !Ja existe um teste anterior bem-sucedido

!===================== Variaveis do laco externo =======================
      INTEGER :: ii                                 !Indice do teste
      INTEGER :: i, j                               !Indices genericos

!===================== Diagnostico em termos.dat =======================
      !Ponto do dominio escolhido para o diagnostico
      INTEGER :: i_dbg, j_dbg, k_dbg, l_dbg, p_dbg, q_dbg, r_dbg

!===================== Cronometragem OpenMP ============================
      DOUBLE PRECISION :: TSTART_EIG, TEND_EIG

        PI = DACOS(-1.D0)       !3.14159265358979323846264338

!===================== Leitura do arquivo de entrada ===================
        CALL GET_COMMAND_ARGUMENT(1, arq_entrada)
        IF (LEN_TRIM(arq_entrada) .EQ. 0) arq_entrada = 'inputs.dat'

        OPEN(UNIT=20, FILE=TRIM(arq_entrada), STATUS='OLD', ACTION='READ', IOSTAT=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'ERRO: nao foi possivel abrir ', TRIM(arq_entrada)
           STOP 1
        END IF

        pasta = ''
        READ(20, NML=saida, IOSTAT=ios, IOMSG=msg_io)
        IF (ios .NE. 0 .OR. LEN_TRIM(pasta) .EQ. 0) THEN
           PRINT *, 'ERRO: bloco &saida com "pasta" ausente ou invalido em ', TRIM(arq_entrada)
           IF (ios .GT. 0) PRINT *, TRIM(msg_io)
           STOP 1
        END IF

        !1a passada: valida TODOS os testes antes de gastar tempo de calculo
        REWIND(20)
        CALL ZERAR_PARAMETROS()
        NTESTES = 0
        DO
           nome = ''
           READ(20, NML=teste, IOSTAT=ios, IOMSG=msg_io)
           IF (ios .LT. 0) EXIT
           IF (ios .GT. 0) THEN
              PRINT '(A,I0,A)', ' ERRO ao ler o bloco &teste numero ', NTESTES+1, ':'
              PRINT *, TRIM(msg_io)
              STOP 1
           END IF
           NTESTES = NTESTES + 1
           IF (NTESTES .GT. NMAX_TESTES) THEN
              PRINT '(A,I0)', ' ERRO: numero de testes maior que NMAX_TESTES = ', NMAX_TESTES
              STOP 1
           END IF
           CALL VALIDAR_E_NOMEAR(NTESTES)
           DO i = 1, NTESTES-1
              IF (nomes(i) .EQ. nome) THEN
                 PRINT '(A,I0,A,I0,A,A)', ' ERRO: testes ', i, ' e ', NTESTES, &
                       ' tem o mesmo nome: ', TRIM(nome)
                 STOP 1
              END IF
           END DO
           nomes(NTESTES) = nome
        END DO

        IF (NTESTES .EQ. 0) THEN
           PRINT *, 'ERRO: nenhum bloco &teste encontrado em ', TRIM(arq_entrada)
           STOP 1
        END IF

!===================== Pasta de saida ==================================
        !Nunca sobrescreve uma rodada anterior
        INQUIRE(FILE=TRIM(pasta)//'/resumo.dat', EXIST=existe)
        IF (existe) THEN
           PRINT *, 'ERRO: ', TRIM(pasta), '/resumo.dat ja existe.'
           PRINT *, '      Escolha outra "pasta" no bloco &saida para nao sobrescrever resultados.'
           STOP 1
        END IF

        CALL EXECUTE_COMMAND_LINE('mkdir -p "'//TRIM(pasta)//'"', EXITSTAT=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'ERRO: nao foi possivel criar a pasta ', TRIM(pasta)
           STOP 1
        END IF
        CALL EXECUTE_COMMAND_LINE('cp "'//TRIM(arq_entrada)//'" "'//TRIM(pasta)//'/inputs.dat"')

        OPEN(UNIT=30, FILE=TRIM(pasta)//'/resumo.dat', STATUS='NEW')
        CALL DATE_AND_TIME(VALUES=data_hora)
        WRITE(30,'(A,I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
              '# Rodada iniciada em ', data_hora(1:3), data_hora(5:7)
        WRITE(30,'(A,A)') '# Entrada: ', TRIM(arq_entrada)
        WRITE(30,'(A)') '# Parametros completos de cada teste em <nome>_parametros.dat'
        WRITE(30,'(A)') '# Norma = pi*Int dxi Int dgamma |Psi|^2 (antes do reescalonamento):'
        WRITE(30,'(A)') '#   preenchida pelo programa de plot; NaN enquanto o plot nao rodar'
        WRITE(30,'(A1,A40,5A6,A3,6A12,2A24,A12,A24)') '#', 'nome', 'NMG', 'NMZ', 'Nz', 'Ng', 'Nv', 'k', &
              'Mtot', 'm1', 'm2', 'mu', 'A_z', 'A_g', 'Alfa', 'wr(1)', 'tempo(s)', 'Norma'
        FLUSH(30)

        !So os alfas, para os testes de convergencia (colunas numericas primeiro)
        OPEN(UNIT=31, FILE=TRIM(pasta)//'/alfas.dat', STATUS='NEW')
        WRITE(31,'(A,I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
              '# Rodada iniciada em ', data_hora(1:3), data_hora(5:7)
        WRITE(31,'(A)') '# dAlfa = Alfa - Alfa do teste anterior (0 no primeiro teste)'
        WRITE(31,'(A1,A3,5A6,2A24,2X,A)') '#', 'i', 'NMG', 'NMZ', 'Nz', 'Ng', 'Nv', &
              'Alfa', 'dAlfa', 'nome'
        FLUSH(31)
        tem_ant = .FALSE.
        Alfa_ant = 0.d0

        PRINT '(A,I0,A,A)', ' ', NTESTES, ' teste(s) validado(s). Saida em: ', TRIM(pasta)

            NCOL = 2

        !Ponto do dominio (z, z', gamma, gamma', v) onde os termos serao impressos
            i_dbg = 10       !indice de gamma   (gv)
            j_dbg = 10       !indice de z       (zv)
            k_dbg = 10       !indice da spline em gamma
            l_dbg = 10      !indice da spline em z
            p_dbg = 10       !indice de Gauss em gamma'
            q_dbg = 10     !indice de Gauss em z'
            r_dbg = 10    !indice de Gauss em v

!===================== 2a passada: executa os testes ===================
        REWIND(20)
        CALL ZERAR_PARAMETROS()

        do ii = 1, NTESTES
          nome = ''
          READ(20, NML=teste)
          CALL VALIDAR_E_NOMEAR(ii)

          WRITE(*,*)
          WRITE(*,'(A,I0,A,I0,A,A)') " ===== Teste ", ii, " de ", NTESTES, ": ", TRIM(nome)
          CALL SYSTEM_CLOCK(T_INI, T_TAXA)

     !Parâmetros derivados
          m = (m1 + m2)/2
          kappa = sqrt(m**2 - 0.25*Mtot**2)
          NMA   = NMG*NMZ
          LWORK = 10*NMA

          ALLOCATE(XMATRIX(NMA,NMA), ZMATRIX(NMA,NMA), c(NMG,NMZ))
          ALLOCATE(ALPHAR(NMA), ALPHAI(NMA), BETA(NMA))
          ALLOCATE(WR(NMA), WI(NMA))
          ALLOCATE(VR(NMA,NMA), VL(NMA,2*NMA), WORK(LWORK))
          ALLOCATE(zv(NMZ+1), gv(NMG+1), XG(NMZ), YG(NMG))

          prefixo = TRIM(pasta)//'/'//TRIM(nome)//'_'
          open (unit = 10, file = TRIM(prefixo)//"autovalores.dat",  STATUS="REPLACE")
          open (unit = 12, file = TRIM(prefixo)//"alfa.dat",         STATUS="REPLACE")
          open (unit = 14, file = TRIM(prefixo)//"erros.dat",        STATUS="REPLACE")
          open (unit = 16, file = TRIM(prefixo)//"coeficientes.dat", STATUS="REPLACE")
          IF (DEBUG_TERMOS) &
          open (unit = 18, file = TRIM(prefixo)//"termos.dat",       STATUS="REPLACE")

      !------------------ Parametros exatos do teste --------------------
      !Este arquivo pode ser usado de novo como bloco &teste de entrada
          open (unit = 19, file = TRIM(prefixo)//"parametros.dat", STATUS="REPLACE")
          CALL DATE_AND_TIME(VALUES=data_hora)
          WRITE(19,'(A,I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
                '! Executado em ', data_hora(1:3), data_hora(5:7)
          WRITE(19,'(A,A,A,I0)') '! Entrada: ', TRIM(arq_entrada), ', bloco &teste numero ', ii
          WRITE(19, NML=teste)
          WRITE(19,'(A)') '! Derivados:'
          WRITE(19,'(A,I0)')      '!   NMA   = ', NMA
          WRITE(19,'(A,ES24.15)') '!   m     = ', m
          WRITE(19,'(A,ES24.15)') '!   kappa = ', kappa
          CLOSE(19)

          WRITE(10, '(A, I0, A, I0, A, I0)') "NMA: ", nma, " NMG: ", nmg, " NMZ: ", nmz
          WRITE(10, '(6(A, F20.10))') &
         & "Mtot: ", Mtot, " mu: ", mu, " kappa: ", kappa, " e: ", e, &
         & " m1: ", m1, " m2: ", m2
          WRITE(10, '(A,A)') " kernel_type: ", kernel_type
          WRITE(10, *) "-----Splines Vetorial Massas Diferentes------------------------"

          WRITE(12, *) "-----Splines Vetorial Massas Diferentes------------------------"
          WRITE(12, '(A, I0, A, I0, A, I0)') "NMA: ", nma, " NMG: ", nmg, " NMZ: ", nmz
          WRITE(12, '(6(A, F20.10))') &
         & "Mtot: ", Mtot, " mu: ", mu, " kappa: ", kappa, " e: ", e, &
         & " m1: ", m1, " m2: ", m2
          WRITE(12, '(A,A)') " kernel_type: ", kernel_type

            iw = 14
            N_intervalZ = (NMZ-1)/2
            N_intervalG = (NMG-1)/2

      !------------------ Cabecalho de termos.dat: parametros ------------------
          IF (DEBUG_TERMOS) THEN
          WRITE(18,*) "========================================================"
          WRITE(18,'(A,I0,A,A)') " TESTE ii = ", ii, "  nome = ", TRIM(nome)
          WRITE(18,*) " Splines Vetorial - Massas Diferentes"
          WRITE(18,*) "========================================================"
          WRITE(18,*) "--- Massas e parametros fisicos ---"
          WRITE(18,'(A,F20.10)') " Mtot  = ", Mtot
          WRITE(18,'(A,F20.10)') " m1    = ", m1
          WRITE(18,'(A,F20.10)') " m2    = ", m2
          WRITE(18,'(A,F20.10)') " m     = ", m
          WRITE(18,'(A,F20.10)') " mu    = ", mu
          WRITE(18,'(A,F20.10)') " kappa = ", kappa
          WRITE(18,'(A,F20.10)') " e     = ", e
          WRITE(18,'(A,F20.10)') " PI    = ", PI
          WRITE(18,'(A,A)')      " kernel_type = ", kernel_type
          WRITE(18,*) "--- Numero de splines (malhas) ---"
          WRITE(18,'(A,I0)') " NMZ (splines em z)     = ", NMZ
          WRITE(18,'(A,I0)') " NMG (splines em gamma) = ", NMG
          WRITE(18,'(A,I0)') " NMA = NMG*NMZ          = ", NMA
          WRITE(18,*) "--- Pontos de Gauss (integracao) ---"
          WRITE(18,'(A,I0)') " Nz (Gauss em z')       = ", Nz
          WRITE(18,'(A,I0)') " Ng (Gauss em gamma')   = ", Ng
          WRITE(18,'(A,I0)') " Nv (Gauss em v)        = ", Nv
          WRITE(18,*) "--------------------------------------------------------"
          END IF

        call G1D(IW,-1.d0, N_intervalZ, A_z, 1.d0, X)
        call COLLOC(IW,2,N_intervalZ,X,XG)

        do i = 1, 2*N_intervalZ
          zv(i+1) = XG(i)
        end do

        zv(1)  = -z_borda
        zv(nmz)=  z_borda

        call G1D(IW,0.d0, N_intervalG, A_g, gamma_max, Y)
        call COLLOC(IW,2,N_intervalG,Y,YG)

        do i=1, 2*N_intervalG
           gv(i+1)=YG(i)
        enddo

        gv(1) = 0.d0
        gv(nmg) = gamma_max

    !Preparação das Splines
        call SPLGR1 (zv,Nmz)
        call SPLGR2 (gv,Nmg)

        !Pesos e absissas de Gauss-Legendre para cada variável
        !X(q) e a variavel NORMALIZADA em [0,1]; o z' fisico em [-1,1] e obtido
        !pelos mapeamentos afins de cada ramo (ver BUILD_ZMATRIX).
        CALL legauss(0.d0,1.d0,Nz,X,dX,1.d-15)


        CALL legauss(0.d0,gamma_max,Ng,Y,dY,1.d-15)
        CALL legauss(0.d0,1.d0,Nv,W,dW,1.d-15)

    !Montagem das matrizes
        CALL BUILD_ZMATRIX(zv, gv, NMZ, NMG, NMA, m, m1, m2, mu, kappa, Mtot, PI, &
                           X, dX, Y, dY, W, dW, Nz, Ng, Nv, &
                           ZMATRIX, XMATRIX, kernel_type, &
                           DEBUG_TERMOS, i_dbg, j_dbg, k_dbg, l_dbg, p_dbg, q_dbg, r_dbg)

          !Condicionar a matriz
            do i = 1, nma
                xmatrix (i,i) = xmatrix (i,i) + e
            end do

        WRITE(*,'(A)') " DGGEV..."

        TSTART_EIG = 0.d0
!$      TSTART_EIG = OMP_GET_WTIME()

        CALL DGGEV('N', 'V', NMA, ZMATRIX, NMA, XMATRIX, NMA, &
                   ALPHAR, ALPHAI, BETA, VL, NMA, VR, NMA, &
                   WORK, LWORK, INFO)

!$      TEND_EIG = OMP_GET_WTIME()
!$      WRITE(*,'(A,F12.3,A)') " Tempo do DGGEV                 = ", &
!$                             TEND_EIG-TSTART_EIG, " s"

        ! Verificação de erro: registra a falha e segue para o proximo teste
        IF (INFO .NE. 0) THEN
            PRINT *, 'ERRO NO DGGEV: INFO = ', INFO, ' (teste ', TRIM(nome), ')'
            WRITE(30,'(A2,A,A,I0)') '# ', TRIM(nome), '  FALHOU: DGGEV INFO = ', INFO
            FLUSH(30)
            WRITE(31,'(A2,A,A,I0)') '# ', TRIM(nome), '  FALHOU: DGGEV INFO = ', INFO
            FLUSH(31)
        ELSE

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

      CALL SYSTEM_CLOCK(T_FIM)
      T_TESTE = DBLE(T_FIM - T_INI)/DBLE(T_TAXA)

      !Uma linha por teste: o programa de plot le os nomes daqui
      !A ultima coluna (Norma, 24 caracteres) e sobrescrita pelo programa de plot
      WRITE(30,'(1X,A40,5I6,2X,A1,6F12.6,2ES24.15,F12.2,A24)') nome, NMG, NMZ, Nz, Ng, Nv, &
            kernel_type, Mtot, m1, m2, mu, A_z, A_g, Alfa, wr(1), T_TESTE, 'NaN'
      FLUSH(30)

      IF (.NOT. tem_ant) Alfa_ant = Alfa
      WRITE(31,'(1X,I3,5I6,2ES24.15,2X,A)') ii, NMG, NMZ, Nz, Ng, Nv, &
            Alfa, Alfa - Alfa_ant, TRIM(nome)
      FLUSH(31)
      Alfa_ant = Alfa
      tem_ant = .TRUE.

      WRITE(*,'(A,F20.10,A,F10.2,A)') " Alfa = ", Alfa, "   (", T_TESTE, " s)"

        END IF

      CLOSE(10)
      CLOSE(12)
      CLOSE(14)
      CLOSE(16)
      IF (DEBUG_TERMOS) CLOSE(18)

      DEALLOCATE(XMATRIX, ZMATRIX, c)
      DEALLOCATE(ALPHAR, ALPHAI, BETA, WR, WI)
      DEALLOCATE(VR, VL, WORK)
      DEALLOCATE(zv, gv, XG, YG)

      end do

      CLOSE(20)
      CLOSE(30)
      CLOSE(31)
      WRITE(*,*)
      WRITE(*,'(A,A)') " Fim. Resumo em ", TRIM(pasta)//'/resumo.dat'

      CONTAINS

! -----------------------------------------------------------------------
! Marca todos os parametros como "nao definidos" antes de ler o arquivo.
! Os blocos &teste seguintes herdam o que nao redefinirem.
! -----------------------------------------------------------------------
      SUBROUTINE ZERAR_PARAMETROS()
        kernel_type = ' '
        NMG = -1;  NMZ = -1
        Nz  = -1;  Ng  = -1;  Nv = -1
        Mtot = NAO_DEF;  m1 = NAO_DEF;  m2 = NAO_DEF;  mu = NAO_DEF
        e    = NAO_DEF
        A_z  = NAO_DEF;  A_g = NAO_DEF
        gamma_max = NAO_DEF;  z_borda = NAO_DEF
        DEBUG_TERMOS = .FALSE.
      END SUBROUTINE ZERAR_PARAMETROS

! -----------------------------------------------------------------------
! Confere o teste corrente e, se "nome" veio vazio, gera um automatico.
! Interrompe o programa com uma mensagem clara em caso de erro.
! -----------------------------------------------------------------------
      SUBROUTINE VALIDAR_E_NOMEAR(it)
        INTEGER, INTENT(IN) :: it
        LOGICAL :: ok

        ok = .TRUE.
        IF (kernel_type .EQ. ' ') CALL ERRO_TESTE(it, 'parametro nao definido: kernel_type', ok)
        IF (NMG .LT. 0)  CALL ERRO_TESTE(it, 'parametro nao definido: NMG', ok)
        IF (NMZ .LT. 0)  CALL ERRO_TESTE(it, 'parametro nao definido: NMZ', ok)
        IF (Nz .LT. 0)   CALL ERRO_TESTE(it, 'parametro nao definido: Nz', ok)
        IF (Ng .LT. 0)   CALL ERRO_TESTE(it, 'parametro nao definido: Ng', ok)
        IF (Nv .LT. 0)   CALL ERRO_TESTE(it, 'parametro nao definido: Nv', ok)
        IF (Mtot .EQ. NAO_DEF) CALL ERRO_TESTE(it, 'parametro nao definido: Mtot', ok)
        IF (m1 .EQ. NAO_DEF)   CALL ERRO_TESTE(it, 'parametro nao definido: m1', ok)
        IF (m2 .EQ. NAO_DEF)   CALL ERRO_TESTE(it, 'parametro nao definido: m2', ok)
        IF (mu .EQ. NAO_DEF)   CALL ERRO_TESTE(it, 'parametro nao definido: mu', ok)
        IF (e .EQ. NAO_DEF)    CALL ERRO_TESTE(it, 'parametro nao definido: e', ok)
        IF (A_z .EQ. NAO_DEF)  CALL ERRO_TESTE(it, 'parametro nao definido: A_z', ok)
        IF (A_g .EQ. NAO_DEF)  CALL ERRO_TESTE(it, 'parametro nao definido: A_g', ok)
        IF (gamma_max .EQ. NAO_DEF) CALL ERRO_TESTE(it, 'parametro nao definido: gamma_max', ok)
        IF (z_borda .EQ. NAO_DEF)   CALL ERRO_TESTE(it, 'parametro nao definido: z_borda', ok)

        IF (ok) THEN
          IF (NMG .LT. 3 .OR. NMG .GT. NMAX_SPL .OR. NMZ .LT. 3 .OR. NMZ .GT. NMAX_SPL) &
             CALL ERRO_TESTE(it, 'NMG e NMZ devem estar entre 3 e 499', ok)
          IF (Nz .LT. 1 .OR. Nz .GT. NMAX_GAUSS .OR. Ng .LT. 1 .OR. Ng .GT. NMAX_GAUSS &
              .OR. Nv .LT. 1 .OR. Nv .GT. NMAX_GAUSS) &
             CALL ERRO_TESTE(it, 'Nz, Ng e Nv devem estar entre 1 e 1000', ok)
          IF (((m1+m2)/2)**2 - 0.25d0*Mtot**2 .LE. 0.d0) &
             CALL ERRO_TESTE(it, 'kappa imaginario: precisa ((m1+m2)/2)**2 > Mtot**2/4', ok)
          IF (z_borda .LE. 0.d0 .OR. z_borda .GT. 1.d0) &
             CALL ERRO_TESTE(it, 'z_borda deve estar em (0, 1]', ok)
          IF (gamma_max .LE. 0.d0) CALL ERRO_TESTE(it, 'gamma_max deve ser positivo', ok)
          IF (A_z .LE. 0.d0 .OR. A_g .LE. 0.d0) CALL ERRO_TESTE(it, 'A_z e A_g devem ser positivos', ok)
        END IF

        IF (LEN_TRIM(nome) .EQ. 0) THEN
           WRITE(nome,'(A,I2.2,5(A,I0),A,A)') 'T', it, '_NMG', NMG, '_NMZ', NMZ, &
                 '_Nz', Nz, '_Ng', Ng, '_Nv', Nv, '_', kernel_type
        ELSE IF (SCAN(TRIM(nome), ' /') .GT. 0) THEN
           CALL ERRO_TESTE(it, 'nome nao pode ter espacos nem "/"', ok)
        END IF

        IF (.NOT. ok) THEN
           PRINT '(A,I0,A)', ' Corrija o bloco &teste numero ', it, ' (ou o primeiro bloco).'
           STOP 1
        END IF

      END SUBROUTINE VALIDAR_E_NOMEAR

      SUBROUTINE ERRO_TESTE(it, texto, ok)
        INTEGER, INTENT(IN) :: it
        CHARACTER(LEN=*), INTENT(IN) :: texto
        LOGICAL, INTENT(INOUT) :: ok
        PRINT '(A,I0,A,A)', ' ERRO no teste ', it, ': ', texto
        ok = .FALSE.
      END SUBROUTINE ERRO_TESTE

    END PROGRAM SPLINE_VETORIAL


! -----------------------------------------------------------------------
! BUILD_ZMATRIX
! Monta ZMATRIX (lado direito, kernel) e XMATRIX (lado esquerdo, norma).
! A avaliacao do kernel fica nas funcoes KERNEL_*; para trocar a fisica
! basta mexer nelas, sem tocar na estrutura dos lacos de quadratura.
! -----------------------------------------------------------------------
      SUBROUTINE BUILD_ZMATRIX(zv, gv, Nmz, Nmg, NMA, m, m1, m2, mu, kappa, Mtot, PI, &
                               X, dX, Y, dY, W, dW, Nz, Ng, Nv, &
                               ZMATRIX, XMATRIX, kernel_type, &
                               DEBUG_TERMOS, i_dbg, j_dbg, k_dbg, l_dbg, p_dbg, q_dbg, r_dbg)
      USE OMP_LIB
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: Nmz, Nmg, NMA, Nz, Ng, Nv
      DOUBLE PRECISION, INTENT(IN) :: zv(*), gv(*)
      DOUBLE PRECISION, INTENT(IN) :: m, m1, m2, mu, kappa, Mtot, PI
      DOUBLE PRECISION, INTENT(IN) :: X(Nz), dX(Nz)
      DOUBLE PRECISION, INTENT(IN) :: Y(Ng), dY(Ng)
      DOUBLE PRECISION, INTENT(IN) :: W(Nv), dW(Nv)
      DOUBLE PRECISION, INTENT(OUT) :: ZMATRIX(NMA,NMA), XMATRIX(NMA,NMA)
      CHARACTER, INTENT(IN) :: kernel_type
      LOGICAL, INTENT(IN) :: DEBUG_TERMOS
      INTEGER, INTENT(IN) :: i_dbg, j_dbg, k_dbg, l_dbg, p_dbg, q_dbg, r_dbg

!===================== Bases de splines no ponto corrente ==============
      DOUBLE PRECISION :: splz(Nmz), splg(Nmg)

!===================== Variaveis do loop de montagem ===================
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
      DOUBLE PRECISION :: Co_ku                     !Termo C0 do numerador no ramo superior
      DOUBLE PRECISION :: Co_kd                     !Termo C0 do numerador no ramo inferior
      DOUBLE PRECISION :: peso                      !Produto splines * jacobianos * pesos
      LOGICAL :: SING_U, SING_D                     !z=-1 (ramo Du) / z=+1 (ramo Dd): 0/0 removivel
      DOUBLE PRECISION, PARAMETER :: TOL_SING = 1.D-10 !Tolerancia em torno de z=+-1
      DOUBLE PRECISION :: contrib_escu, contrib_escd  !Contribuicoes escalares de cada ramo
      DOUBLE PRECISION :: contrib_C0_kd, contrib_f1_kd, contrib_C0_ku, contrib_f1_ku
      LOGICAL :: printou_termos

!===================== Funcoes de kernel ===============================
      DOUBLE PRECISION, EXTERNAL :: KERNEL_UPPER_ESC, KERNEL_LOWER_ESC
      DOUBLE PRECISION, EXTERNAL :: KERNEL_UPPER_C0,  KERNEL_LOWER_C0
      DOUBLE PRECISION, EXTERNAL :: KERNEL_UPPER_F1,  KERNEL_LOWER_F1
      DOUBLE PRECISION, EXTERNAL :: DENOM_D0, DENOM_DU, DENOM_DD
      DOUBLE PRECISION, EXTERNAL :: COEF_C0_U, COEF_C0_D

!===================== Cronometragem OpenMP ============================
      DOUBLE PRECISION :: TSTART, TEND
      INTEGER :: NTHREADS

!===================== Acompanhamento do progresso =====================
      INTEGER :: NFEITOS                            !Pares (i,j) ja concluidos (COMPARTILHADO)
      INTEGER :: NFEITOS_LOC                        !Copia local do contador (PRIVATE)
      INTEGER :: NTOTAL                             !Total de pares = Nmg*Nmz
      INTEGER :: NPASSO                             !Imprime a cada NPASSO pares concluidos
      DOUBLE PRECISION :: TAGORA, TDECOR, TRESTA    !Tempos para a estimativa
      DOUBLE PRECISION :: PCT                       !Percentual concluido

        zmatrix  = 0.d0
        xmatrix  = 0.d0
        printou_termos = .FALSE.

        NTHREADS = 1
!$      NTHREADS = OMP_GET_MAX_THREADS()

        WRITE(*,'(A,I0,A)') " Montagem das matrizes com ", NTHREADS, &
                             " thread(s) OpenMP..."
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
!$OMP           g, z, gp, dgp, v, dv, zq, dzq, peso, &
!$OMP           SPLz, SPLg, &
!$OMP           D0, Du, Dd, Co_ku, Co_kd, SING_U, SING_D, &
!$OMP           contrib_escu, contrib_escd, &
!$OMP           contrib_C0_kd, contrib_f1_kd, contrib_C0_ku, contrib_f1_ku, &
!$OMP           NFEITOS_LOC, TAGORA, TDECOR, TRESTA, PCT) &
!$OMP   COLLAPSE(2) SCHEDULE(DYNAMIC)
        !Loop para cada elemento da Matriz
        do i=1,Nmg
           do j=1,Nmz
              g=gv(i)
              z=zv(j)
              !Du ~ (z+1) e Dd ~ (z-1): o numerador tem o mesmo fator, entao
              !a razao e uma singularidade removivel (0/0) em z=-1 e z=+1.
              !Como agora zv(1)=-1 e zv(nmz)=+1 sao nos exatos da malha, a
              !divisao literal por Du/Dd nesses pontos e numericamente instavel
              !(cancelamento catastrofico); a contribuicao correta e finita
              !mas pulamos o ramo aqui, como faz o codigo frances em
              !Kernel_droite (if(dabs(u∓1)<1.D-10) return).
              SING_U = (DABS(z + 1.d0) .LT. TOL_SING)
              SING_D = (DABS(z - 1.d0) .LT. TOL_SING)
              index1 = (j-1)*Nmg + i           !Juntei cada ponto (gi, zj) num vetor coluna de dimensão Nmg*Nmz
              do k=1,Nmg
                do l=1, Nmz
                    index2 = (l-1)*Nmg + k    !Juntei cada Iteração das Splines Sg e Sz da integração em um vetor coluna de dimenção Nmg*Nmz

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

        ! theta z’ variando de z até 1
                     dzq=(1.d0-z)*dX(q)
                     zq=(1.d0-z)*X(q)+z

                     call SPLMD1 (zv,Nmz,zq,SPLz)
                     call SPLMD2 (gv,Nmg,gp,SPLg)

        !Peso comum a todas as contribuicoes (splines + jacobianos + pesos de Gauss)
                     peso = splg(k)*splz(l)*dzq*dgp*dv

        !----------------- Ramo superior: avaliacao do kernel -----------------
                     IF (kernel_type .EQ. 'v') THEN
                        contrib_C0_ku = KERNEL_UPPER_C0(z, zq, g, gp, v, m1, m2, mu, Mtot, PI) * peso
                        contrib_f1_ku = KERNEL_UPPER_F1(z, zq, g, gp, v, m1, m2, mu, Mtot, PI) * peso
                        contrib_escu  = 0.d0
                        zmatrix (index1, index2) = zmatrix (index1, index2) &
                                                   - contrib_C0_ku + contrib_f1_ku
                     ELSE
                        contrib_escu  = KERNEL_UPPER_ESC(z, zq, g, gp, v, m1, m2, mu, Mtot, PI) * peso
                        contrib_C0_ku = 0.d0
                        contrib_f1_ku = 0.d0
                        zmatrix (index1, index2) = zmatrix (index1, index2) + contrib_escu
                     END IF

      !------- Impressao dos termos para o ponto escolhido do dominio -------
                    IF (DEBUG_TERMOS .AND. .NOT. SING_U .AND. .NOT. printou_termos &
                        .AND. i.EQ.i_dbg .AND. j.EQ.j_dbg &
                        .AND. k.EQ.k_dbg .AND. l.EQ.l_dbg .AND. p.EQ.p_dbg &
                        .AND. q.EQ.q_dbg .AND. r.EQ.r_dbg) THEN
                      !Denominadores e coeficientes so sao recalculados aqui,
                      !fora do caminho quente, para o relatorio de diagnostico.
                      D0    = DENOM_D0(z, g, m1, m2, Mtot)
                      Du    = DENOM_DU(z, zq, g, gp, v, m1, m2, mu, Mtot)
                      Co_ku = COEF_C0_U(z, zq, g, v, m2, Mtot)
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
                      WRITE(18,'(A,ES24.15)') " Co_ku   = ", Co_ku
                      WRITE(18,'(A,ES24.15)') " splg(k) = ", splg(k)
                      WRITE(18,'(A,ES24.15)') " splz(l) = ", splz(l)
                      WRITE(18,'(A,ES24.15)') " (1+z)**2                = ", (1+z)**2
                      WRITE(18,'(A,ES24.15)') " 32*PI**2*D0             = ", 32*PI**2*D0
                      WRITE(18,'(A,ES24.15)') " v**2/Du**2              = ", v**2/(Du**2)
                      WRITE(18,'(A,ES24.15)') " contribuicao escalar    = ", contrib_escu
                      WRITE(18,*) ""
                      WRITE(18,*) "--- contrib_C0_ku e contrib_f1_ku (ramo superior) ---"
                      WRITE(18,'(A,ES24.15)') " contrib_C0_ku (com splines/jacob.)    = ", contrib_C0_ku
                      WRITE(18,'(A,ES24.15)') " contrib_C0_ku (sem splines/jacob.)    = ", &
                          KERNEL_UPPER_C0(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
                      WRITE(18,'(A,ES24.15)') " contrib_f1_ku (com splines/jacob.)    = ", contrib_f1_ku
                      WRITE(18,'(A,ES24.15)') " contrib_f1_ku (sem splines/jacob.)    = ", &
                          KERNEL_UPPER_F1(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
                      WRITE(18,'(A,ES24.15)') " soma contrib_C0_ku + contrib_f1_ku    = ", &
                          contrib_C0_ku + contrib_f1_ku
!$OMP END CRITICAL (TERMOS_DBG)
                    END IF

        ! theta z’ variando de -1 até z
                     dzq=(z+1.d0)*dx(q)
                     zq=(z+1.d0)*x(q)-1.d0
                     call SPLMD1 (zv,Nmz,zq,SPLz)

                     peso = splg(k)*splz(l)*dzq*dgp*dv

        !----------------- Ramo inferior: avaliacao do kernel -----------------
                     IF (kernel_type .EQ. 'v') THEN
                        contrib_C0_kd = KERNEL_LOWER_C0(z, zq, g, gp, v, m1, m2, mu, Mtot, PI) * peso
                        contrib_f1_kd = KERNEL_LOWER_F1(z, zq, g, gp, v, m1, m2, mu, Mtot, PI) * peso
                        contrib_escd  = 0.d0
                        zmatrix (index1, index2) = zmatrix (index1, index2) &
                                                   - contrib_C0_kd + contrib_f1_kd
                     ELSE
                        contrib_escd  = KERNEL_LOWER_ESC(z, zq, g, gp, v, m1, m2, mu, Mtot, PI) * peso
                        contrib_C0_kd = 0.d0
                        contrib_f1_kd = 0.d0
                        zmatrix (index1, index2) = zmatrix (index1, index2) + contrib_escd
                     END IF

      !------- Impressao dos termos para o ponto escolhido do dominio -------
                    IF (DEBUG_TERMOS .AND. .NOT. SING_D .AND. .NOT. printou_termos &
                        .AND. i.EQ.i_dbg .AND. j.EQ.j_dbg &
                        .AND. k.EQ.k_dbg .AND. l.EQ.l_dbg .AND. p.EQ.p_dbg &
                        .AND. q.EQ.q_dbg .AND. r.EQ.r_dbg) THEN
                      D0    = DENOM_D0(z, g, m1, m2, Mtot)
                      Dd    = DENOM_DD(z, zq, g, gp, v, m1, m2, mu, Mtot)
                      Co_kd = COEF_C0_D(z, zq, g, v, m1, Mtot)
!$OMP CRITICAL (TERMOS_DBG)
                      WRITE(18,*) ""
                      WRITE(18,*) "=== RAMO INFERIOR: z' de -1 ate z ==="
                      WRITE(18,'(A,ES24.15)') " zq      = ", zq
                      WRITE(18,'(A,ES24.15)') " dzq     = ", dzq
                      WRITE(18,'(A,ES24.15)') " D0      = ", D0
                      WRITE(18,'(A,ES24.15)') " Dd      = ", Dd
                      WRITE(18,'(A,ES24.15)') " Co_kd   = ", Co_kd
                      WRITE(18,'(A,ES24.15)') " splg(k) = ", splg(k)
                      WRITE(18,'(A,ES24.15)') " splz(l) = ", splz(l)
                      WRITE(18,'(A,ES24.15)') " (1-z)**2                = ", (1-z)**2
                      WRITE(18,'(A,ES24.15)') " 32*PI**2*D0             = ", 32*PI**2*D0
                      WRITE(18,'(A,ES24.15)') " v**2/Dd**2              = ", v**2/(Dd**2)
                      WRITE(18,'(A,ES24.15)') " contribuicao escalar    = ", contrib_escd
                      WRITE(18,*) ""
                      WRITE(18,*) "--- contrib_C0_kd e contrib_f1_kd (ramo inferior) ---"
                      WRITE(18,'(A,ES24.15)') " contrib_C0_kd (com splines/jacob.)    = ", contrib_C0_kd
                      WRITE(18,'(A,ES24.15)') " contrib_C0_kd (sem splines/jacob.)    = ", &
                          KERNEL_LOWER_C0(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
                      WRITE(18,'(A,ES24.15)') " contrib_f1_kd (com splines/jacob.)    = ", contrib_f1_kd
                      WRITE(18,'(A,ES24.15)') " contrib_f1_kd (sem splines/jacob.)    = ", &
                          KERNEL_LOWER_F1(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
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

        TEND = 0.d0
!$      TEND = OMP_GET_WTIME()
!$      WRITE(*,'(A,F12.3,A)') " Tempo de montagem das matrizes = ", &
!$                             TEND-TSTART, " s"

      END SUBROUTINE BUILD_ZMATRIX


! =======================================================================
!                    DENOMINADORES COMUNS DO KERNEL
! =======================================================================

! -----------------------------------------------------------------------
! D0 — denominador comum aos dois ramos (funcao apenas de z e gamma)
! -----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DENOM_D0(z, g, m1, m2, Mtot)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, g, m1, m2, Mtot

      DENOM_D0 = 0.25d0*(4.d0*g + Mtot**2*(z**2 - 1.d0) &
                 + 2.d0*m1**2*(z + 1.d0) - 2.d0*m2**2*(z - 1.d0))
      END FUNCTION DENOM_D0


! -----------------------------------------------------------------------
! Du — denominador do ramo superior (z' em [z, 1])
! -----------------------------------------------------------------------
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


! -----------------------------------------------------------------------
! Dd — denominador do ramo inferior (z' em [-1, z])
! -----------------------------------------------------------------------
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

! -----------------------------------------------------------------------
! C0 do ramo superior
! -----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION COEF_C0_U(z, zq, g, v, m2, Mtot)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, v, m2, Mtot

      !COEF_C0_U = (-1.d0)*0.25d0*(z + 1.d0) &
          !*(Mtot**2*(z + 1.d0)*((v*zq + v - 2.d0)*((v - 2.d0)*z - v*zq) + 4.d0) &
            !+ 4.d0*(v - 2.d0)*(m2**2*(-v*z + v*zq + 2.d0*z) + g*(v*zq + v - 2.d0)))

      COEF_C0_U = (-2 + v)*(1 + z)*(m2**2*((-2 + v)*z - v*zq) - g*(-2 + v + v*zq)) - & 
                    (Mtot**2*(1 + z)**2*(4 + ((-2 + v)*z - v*zq)*(-2 + v + v*zq)))/4.
      END FUNCTION COEF_C0_U


! -----------------------------------------------------------------------
! C0 do ramo inferior
! -----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION COEF_C0_D(z, zq, g, v, m1, Mtot)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, v, m1, Mtot

     ! COEF_C0_D = (-1.d0)*0.25d0*(z - 1.d0) &
         ! *(Mtot**2*(z - 1.d0)*((v*(zq - 1.d0) + 2.d0)*((v - 2.d0)*z - v*zq) + 4.d0) &
           ! + 4.d0*(v - 2.d0)*(m1**2*(-v*z + v*zq + 2.d0*z) + g*(v*(zq - 1.d0) + 2.d0)))

      COEF_C0_D = (-2 + v)*(-1 + z)*(g*(-2 + v - v*zq) + m1**2*((-2 + v)*z - v*zq)) - &
            (Mtot**2*(-1 + z)**2*(4 + (2 + v*(-1 + zq))*((-2 + v)*z - v*zq)))/4.
      END FUNCTION COEF_C0_D


! =======================================================================
!                        KERNEL ESCALAR
! Numerador (1 +- z)**2; nao ha termo f1. Mesma estrutura dos kernels
! vetoriais, de modo que o builder trate os dois casos igualmente.
! =======================================================================

! -----------------------------------------------------------------------
! KERNEL_UPPER_ESC — integrando escalar para z' em [z, 1]
! -----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION KERNEL_UPPER_ESC(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Du
      DOUBLE PRECISION, EXTERNAL :: DENOM_D0, DENOM_DU

      D0 = DENOM_D0(z, g, m1, m2, Mtot)
      Du = DENOM_DU(z, zq, g, gp, v, m1, m2, mu, Mtot)

      KERNEL_UPPER_ESC = 1.d0 / (32*PI**2*D0) * (v**2 / (Du**2)) * (1.d0 + z)**2
      END FUNCTION KERNEL_UPPER_ESC


! -----------------------------------------------------------------------
! KERNEL_LOWER_ESC — integrando escalar para z' em [-1, z]
! -----------------------------------------------------------------------
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
! O numerador se separa em duas partes, tratadas por funcoes distintas:
!   C0 -> pesa 1/Du**2 (ou 1/Dd**2)
!   f1 -> pesa 1/Du    (ou 1/Dd),  com fator 2
! O builder soma -(C0 + f1) a zmatrix.
! =======================================================================

! -----------------------------------------------------------------------
! KERNEL_UPPER_C0 — parte C0 do kernel vetorial, z' em [z, 1]
! -----------------------------------------------------------------------
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


! -----------------------------------------------------------------------
! KERNEL_UPPER_F1 — parte f1 do kernel vetorial, z' em [z, 1]
! -----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION KERNEL_UPPER_F1(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Du
      DOUBLE PRECISION, EXTERNAL :: DENOM_D0, DENOM_DU

      D0 = DENOM_D0(z, g, m1, m2, Mtot)
      Du = DENOM_DU(z, zq, g, gp, v, m1, m2, mu, Mtot)

      KERNEL_UPPER_F1 = 2.d0 / (32*PI**2*D0) * (v**2 / Du) * (z + 1.d0)
      END FUNCTION KERNEL_UPPER_F1


! -----------------------------------------------------------------------
! KERNEL_LOWER_C0 — parte C0 do kernel vetorial, z' em [-1, z]
! -----------------------------------------------------------------------
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


! -----------------------------------------------------------------------
! KERNEL_LOWER_F1 — parte f1 do kernel vetorial, z' em [-1, z]
! -----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION KERNEL_LOWER_F1(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Dd
      DOUBLE PRECISION, EXTERNAL :: DENOM_D0, DENOM_DD

      D0 = DENOM_D0(z, g, m1, m2, Mtot)
      Dd = DENOM_DD(z, zq, g, gp, v, m1, m2, mu, Mtot)

      KERNEL_LOWER_F1 = 2.d0 / (32*PI**2*D0) * (v**2 / Dd) * (1.d0 - z)
      END FUNCTION KERNEL_LOWER_F1

    
    
    
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
