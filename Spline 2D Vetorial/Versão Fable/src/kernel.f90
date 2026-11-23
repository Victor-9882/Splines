! =======================================================================
!              KERNEL DO PROPAGADOR VETORIAL: P^(1) E P^(2)
!
! Incluido (INCLUDE) no inicio de 2D_Spline_Vetorial_Diferentes_zeta_omp.f90
! e no driver de testes (testes_zeta/). Todas as rotinas sao PURE: sem
! SAVE e sem estado global, podem ser chamadas de dentro do laco OpenMP.
! Ficam em MODULEs para que as interfaces PURE sejam explicitas.
!
! Lado direito da EBS:  -lambda^2/(32 pi^2) Int dz' Int dgamma' (P1 - (1-zeta) P2) g
!
!   MODULE KERNEL_P1      -> termo g_munu (P^(1)) e kernel escalar,
!                            usado por BUILD_ZMATRIX e TERMOS_P1
!   MODULE KERNEL_P2_ZETA -> termo (1-zeta) (P^(2)),
!                            usado por BUILD_ZMATRIX_P2 e TERMOS_P2
! =======================================================================


! =======================================================================
!                     P^(1): TERMO g_munu (zeta = 1)
!
! Ramo superior (z' em [z,1]) e inferior (z' em [-1,z]):
!   vetorial: -C0 * v^2/(32 pi^2 D0 Du^2) + 2 v^2 (1+z)/(32 pi^2 D0 Du)
!             (idem com Dd, Co_kd e (1-z) no ramo inferior)
!   escalar :  v^2 (1+-z)^2 / (32 pi^2 D0 Du^2)   (Dd no ramo inferior)
!
! Convencoes: Mtot -> M, g -> gamma, gp -> gamma', zq -> z'.
! =======================================================================
      MODULE KERNEL_P1
      IMPLICIT NONE
      CONTAINS

! =======================================================================
!                    DENOMINADORES COMUNS DO KERNEL
! =======================================================================

! -----------------------------------------------------------------------
! D0 — denominador comum aos dois ramos (funcao apenas de z e gamma)
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION DENOM_D0(z, g, m1, m2, Mtot)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, g, m1, m2, Mtot

      DENOM_D0 = 0.25d0*(4.d0*g + Mtot**2*(z**2 - 1.d0) &
                 + 2.d0*m1**2*(z + 1.d0) - 2.d0*m2**2*(z - 1.d0))
      END FUNCTION DENOM_D0


! -----------------------------------------------------------------------
! Du — denominador do ramo superior (z' em [z, 1])
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION DENOM_DU(z, zq, g, gp, v, m1, m2, mu, Mtot)
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
      PURE DOUBLE PRECISION FUNCTION DENOM_DD(z, zq, g, gp, v, m1, m2, mu, Mtot)
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
      PURE DOUBLE PRECISION FUNCTION COEF_C0_U(z, zq, g, v, m2, Mtot)
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
      PURE DOUBLE PRECISION FUNCTION COEF_C0_D(z, zq, g, v, m1, Mtot)
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
      PURE DOUBLE PRECISION FUNCTION KERNEL_UPPER_ESC(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Du

      D0 = DENOM_D0(z, g, m1, m2, Mtot)
      Du = DENOM_DU(z, zq, g, gp, v, m1, m2, mu, Mtot)

      KERNEL_UPPER_ESC = 1.d0 / (32*PI**2*D0) * (v**2 / (Du**2)) * (1.d0 + z)**2
      END FUNCTION KERNEL_UPPER_ESC


! -----------------------------------------------------------------------
! KERNEL_LOWER_ESC — integrando escalar para z' em [-1, z]
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION KERNEL_LOWER_ESC(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Dd

      D0 = DENOM_D0(z, g, m1, m2, Mtot)
      Dd = DENOM_DD(z, zq, g, gp, v, m1, m2, mu, Mtot)

      KERNEL_LOWER_ESC = 1.d0 / (32*PI**2*D0) * (v**2 / (Dd**2)) * (1.d0 - z)**2
      END FUNCTION KERNEL_LOWER_ESC


! =======================================================================
!                        KERNEL VETORIAL
! O numerador se separa em duas partes, tratadas por funcoes distintas:
!   C0 -> pesa 1/Du**2 (ou 1/Dd**2)
!   f1 -> pesa 1/Du    (ou 1/Dd),  com fator 2
! O builder soma -C0 + f1 a zmatrix.
! =======================================================================

! -----------------------------------------------------------------------
! KERNEL_UPPER_C0 — parte C0 do kernel vetorial, z' em [z, 1]
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION KERNEL_UPPER_C0(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Du, Co_ku

      D0    = DENOM_D0(z, g, m1, m2, Mtot)
      Du    = DENOM_DU(z, zq, g, gp, v, m1, m2, mu, Mtot)
      Co_ku = COEF_C0_U(z, zq, g, v, m2, Mtot)

      KERNEL_UPPER_C0 = 1.d0 / (32*PI**2*D0) * (v**2 / (Du**2)) * Co_ku
      END FUNCTION KERNEL_UPPER_C0


! -----------------------------------------------------------------------
! KERNEL_UPPER_F1 — parte f1 do kernel vetorial, z' em [z, 1]
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION KERNEL_UPPER_F1(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Du

      D0 = DENOM_D0(z, g, m1, m2, Mtot)
      Du = DENOM_DU(z, zq, g, gp, v, m1, m2, mu, Mtot)

      KERNEL_UPPER_F1 = 2.d0 / (32*PI**2*D0) * (v**2 / Du) * (z + 1.d0)
      END FUNCTION KERNEL_UPPER_F1


! -----------------------------------------------------------------------
! KERNEL_LOWER_C0 — parte C0 do kernel vetorial, z' em [-1, z]
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION KERNEL_LOWER_C0(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Dd, Co_kd

      D0    = DENOM_D0(z, g, m1, m2, Mtot)
      Dd    = DENOM_DD(z, zq, g, gp, v, m1, m2, mu, Mtot)
      Co_kd = COEF_C0_D(z, zq, g, v, m1, Mtot)

      KERNEL_LOWER_C0 = 1.d0 / (32*PI**2*D0) * (v**2 / (Dd**2)) * Co_kd
      END FUNCTION KERNEL_LOWER_C0


! -----------------------------------------------------------------------
! KERNEL_LOWER_F1 — parte f1 do kernel vetorial, z' em [-1, z]
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION KERNEL_LOWER_F1(z, zq, g, gp, v, m1, m2, mu, Mtot, PI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zq, g, gp, v, m1, m2, mu, Mtot, PI
      DOUBLE PRECISION :: D0, Dd

      D0 = DENOM_D0(z, g, m1, m2, Mtot)
      Dd = DENOM_DD(z, zq, g, gp, v, m1, m2, mu, Mtot)

      KERNEL_LOWER_F1 = 2.d0 / (32*PI**2*D0) * (v**2 / Dd) * (1.d0 - z)
      END FUNCTION KERNEL_LOWER_F1

      END MODULE KERNEL_P1


! =======================================================================
!                  P^(2): TERMO (1-zeta) DO PROPAGADOR
!
!   P2 = P2_reg + P2_delta
!
!   P2_reg   = -1/2 Int_0^1 dv v^2 Int_0^{1-v} du
!                [ N_ku (1+z)^3 th(z'-z) / (D0 Du2^3) + N_kd (1-z)^3 th(z-z') / (D0 Dd2^3) ]
!
!   Para nao haver 1/(1+z) nem 1/(1-z), o numerador e escrito com os
!   polos reescalados k~_u = (1+z) k_u, k~_d = (1-z) k_d. Como N e
!   quadratico em k^-, N~ = c^2 N (c = 1+z ou 1-z) e um polinomio em
!   (k~, c); dos tres fatores (1+-z)^3 sobra um explicito:
!
!   P2_reg   = -1/2 Int dv v^2 Int du [ (1+z) N~_u th(z'-z) / (D0 Du2^3)
!                                     + (1-z) N~_d th(z-z') / (D0 Dd2^3) ]
!
!   Du2, Dd2: DENOM_DU2, DENOM_DD2 (finitos e > 0; Du2 = -(1+z) f_u).
!
!   P2_delta = 2 delta(z'-z) Int_0^1 dv v^3 / ( a(v) [a(v) + (1-zeta)(1-v) mu^2] )
!
! Convencoes: Mtot -> M, g -> gamma, gp -> gamma', zp -> z'.
! =======================================================================
      MODULE KERNEL_P2_ZETA
      IMPLICIT NONE
      CONTAINS

! -----------------------------------------------------------------------
! Polos dos propagadores escalares em k^-, ja multiplicados por (1+-z):
!   k_u = M/2 - 2(gamma+m2^2)/(M(1+z)),   k_d = -M/2 + 2(gamma+m1^2)/(M(1-z))
!   kut = (1+z) k_u = (M^2(1+z) - 4(gamma+m2^2))/(2M)
!   kdt = (1-z) k_d = (M^2(z-1) + 4(gamma+m1^2))/(2M)
! Sem singularidade em z = -+1.
! -----------------------------------------------------------------------
      PURE SUBROUTINE POLOS_KMENOS(z, g, m1, m2, Mtot, kut, kdt)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: z, g, m1, m2, Mtot
      DOUBLE PRECISION, INTENT(OUT) :: kut, kdt

      kut =  (-4*g - 4*m2**2 + Mtot**2*(1 + z))/(2.*Mtot)
      kdt =  (4*g + 4*m1**2 + Mtot**2*(-1 + z))/(2.*Mtot)
      END SUBROUTINE POLOS_KMENOS


! -----------------------------------------------------------------------
! D0 = gamma + (1-z^2) kappa^2 + (Delta - z mbar)^2
! (mesmo valor de DENOM_D0, escrito na forma da derivacao)
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION D0_P2(z, g, m1, m2, Mtot)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, g, m1, m2, Mtot
      DOUBLE PRECISION :: mb, Dl, k2

      mb = 0.5d0*(m1 + m2)
      Dl = 0.5d0*(m2 - m1)
      k2 = mb**2 - 0.25d0*Mtot**2
      D0_P2 = g + (1.d0 - z**2)*k2 + (Dl - z*mb)**2
      END FUNCTION D0_P2


! -----------------------------------------------------------------------
! K_D = (M/2) v (1-v) (z'-z)   (coeficiente de k^- no denominador)
! KD_P2 e LD_P2 nao entram em P2_REG (que usa DENOM_DU2/DD2); ficam para
! diagnostico (TERMOS_P2) e testes: -(K_D k~ + c l_D) = D2.
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION KD_P2(z, zp, v, Mtot)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zp, v, Mtot

      KD_P2 = 0.5d0*Mtot*v*(1.d0 - v)*(zp - z)
      END FUNCTION KD_P2


! -----------------------------------------------------------------------
! l_D^(2) = -v(1-v)(gamma + M^2 z z'/4) - v(kappa^2 + M^2 v z'^2/4 + gamma')
!           - (1-v) mu^2 + mu^2 (1-zeta) u
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION LD_P2(z, zp, g, gp, v, u, m1, m2, mu, Mtot, zeta)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zp, g, gp, v, u, m1, m2, mu, Mtot, zeta
      DOUBLE PRECISION :: k2

      k2 = (0.5d0*(m1 + m2))**2 - 0.25d0*Mtot**2
      LD_P2 = -v*(1.d0 - v)*(g + 0.25d0*Mtot**2*z*zp) &
              - v*(k2 + 0.25d0*Mtot**2*v*zp**2 + gp) &
              - (1.d0 - v)*mu**2 + mu**2*(1.d0 - zeta)*u
      END FUNCTION LD_P2


! -----------------------------------------------------------------------
! Du2 — denominador do ramo superior (z' em [z, 1]) do termo (1-zeta)
!   Du2 = DENOM_DU - mu^2 (1-zeta) (1+z) u = -(K_D k~_u + (1+z) l_D)
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION DENOM_DU2(z, zl, g, gp, v, u, m1, m2, mu, Mtot, zeta)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zl, g, gp, v, u, m1, m2, mu, Mtot, zeta

      DENOM_DU2 = (-4*mu**2*(1 + z)*(-1 + u + v - u*zeta) &
                   - Mtot**2*v*(1 + z)*(1 + zl)*(1 + (-1 + v)*z - v*zl) &
                   + v*(4*gp*(1 + z) + m1**2*(1 + z) + 2*m1*m2*(1 + z) - 4*g*(-1 + v)*(1 + zl) &
                        + m2**2*(1 - 3*z + 4*v*z + 4*zl - 4*v*zl)))/4.d0
      END FUNCTION DENOM_DU2


! -----------------------------------------------------------------------
! Dd2 — denominador do ramo inferior (z' em [-1, z]) do termo (1-zeta)
!   Dd2 = DENOM_DD - mu^2 (1-zeta) (1-z) u = -(K_D k~_d + (1-z) l_D)
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION DENOM_DD2(z, zl, g, gp, v, u, m1, m2, mu, Mtot, zeta)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: z, zl, g, gp, v, u, m1, m2, mu, Mtot, zeta

      DENOM_DD2 = (4*mu**2*(-1 + z)*(-1 + u + v - u*zeta) &
                   + Mtot**2*v*(-1 + z)*(-1 + zl)*(-1 + (-1 + v)*z - v*zl) &
                   + v*(m2**2 - 2*m1*m2*(-1 + z) - 4*gp*z - m2**2*z + 4*(gp + g*(-1 + v)*(-1 + zl)) &
                        + m1**2*(1 + 3*z - 4*v*z + 4*(-1 + v)*zl)))/4.d0
      END FUNCTION DENOM_DD2


! -----------------------------------------------------------------------
! Numerador reescalado N~ = c^2 N(k^-) do termo (1-zeta), escrito com
! kmt = c k^- e o denominador D2 = Du2 ou Dd2 (c f = -D2), sendo
! c = 1+z no ramo superior e 1-z no inferior:
!   c k^2 = -z (M/2) kmt - c gamma,   c k.p = (M/2) kmt - c z M^2/4
!   N     = 4 d0 + f ( kp d1 + 4 d2 + M^2 d3 + k2 d4 + 12 f d5 )
!   N~    = 4 c^2 d0 - D2 ( c kp d1 + 4 c d2 + c M^2 d3 + c k2 d4 - 12 D2 d5 )
! d0 e quadratico e d2 linear em k^-, de modo que cada termo de N~ e um
! polinomio em (kmt, c): nenhuma divisao por c. Com c = 1 reduz a N(k^-).
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION NUM_P2(kmt, c, z, zp, g, v, D2, Mtot)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: kmt, c, z, zp, g, v, D2, Mtot
      DOUBLE PRECISION :: k2t, kpt, d0t, d1, d2t, d3, d4, d5, M

      M   = Mtot
      k2t = -z*0.5d0*M*kmt - c*g                 ! c k^2
      kpt = 0.5d0*M*kmt - c*0.25d0*z*M**2        ! c k.p

      ! c^2 d0, c d2
      d0t = -(v**2/16.d0) * ( 4.d0*((v-1.d0)*zp - 1.d0)*kpt + 4.d0*(v-2.d0)*k2t + c*M**2*zp*(v*zp - 2.d0) ) &
                          * ( 4.d0*((v-1.d0)*zp + 1.d0)*kpt + 4.d0*(v-2.d0)*k2t + c*M**2*zp*(v*zp + 2.d0) )
      d1  = -4.d0*(v-1.d0)*v*zp
      d2t = -(v/2.d0) * ( zp*(4.d0*(v-1.d0)*kpt + c*M**2*v*zp) + 4.d0*(v-2.d0)*k2t )
      d3  = 1.d0 - v**2*zp**2
      d4  = -4.d0*(v-1.d0)**2
      d5  = -1.d0

      NUM_P2 = 4.d0*d0t - D2*( kpt*d1 + 4.d0*d2t + c*M**2*d3 + k2t*d4 - 12.d0*D2*d5 )
      END FUNCTION NUM_P2


! -----------------------------------------------------------------------
! P2_REG — parte regular de P^(2) num ponto (gamma, z, gamma', z')
!   ramo = +1 : z' em [z, 1]  (theta(z'-z)), polo k_u
!   ramo = -1 : z' em [-1, z] (theta(z-z')), polo k_d
!   valor = -1/(2 D0) Int dv v^2 Int du c N~/D2^3
!   c = 1+z, D2 = Du2 (ramo +1)  ou  c = 1-z, D2 = Dd2 (ramo -1); N~ = c^2 N
! v: Gauss (W, dW) em [0,1] (o mesmo de P^(1)); u = (1-v) t, t: Gauss (T, dT).
! ok = .FALSE. se D2 <= 0 ou algum valor nao finito aparecer (nao mascara).
! Se c = 0 o ramo vale 0 (fator c explicito; N~ e D2 sao finitos), que e
! o limite exato; so acontece com z_borda = 1. O RETURN evita 0/0 quando
! tambem D2 = 0 (z' = z = -+1).
! -----------------------------------------------------------------------
      PURE SUBROUTINE P2_REG(ramo, g, z, gp, zp, m1, m2, mu, Mtot, zeta, &
                             W, dW, Nv, T, dT, Nu, valor, ok)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ramo, Nv, Nu
      DOUBLE PRECISION, INTENT(IN) :: g, z, gp, zp, m1, m2, mu, Mtot, zeta
      DOUBLE PRECISION, INTENT(IN) :: W(Nv), dW(Nv), T(Nu), dT(Nu)
      DOUBLE PRECISION, INTENT(OUT) :: valor
      LOGICAL, INTENT(OUT) :: ok
      INTEGER :: r, s
      DOUBLE PRECISION :: kut, kdt, kmt, c, D0, v, u, D2, N, soma, soma_u

      valor = 0.d0
      ok = .TRUE.
      CALL POLOS_KMENOS(z, g, m1, m2, Mtot, kut, kdt)
      IF (ramo .EQ. 1) THEN
         c   = 1.d0 + z
         kmt = kut
      ELSE
         c   = 1.d0 - z
         kmt = kdt
      END IF
      IF (c .EQ. 0.d0) RETURN
      D0 = D0_P2(z, g, m1, m2, Mtot)

      soma = 0.d0
      do r = 1, Nv !integral em v
         v  = W(r)
         soma_u = 0.d0
         do s = 1, Nu !integral em u
            u  = (1.d0 - v)*T(s)
            IF (ramo .EQ. 1) THEN
               D2 = DENOM_DU2(z, zp, g, gp, v, u, m1, m2, mu, Mtot, zeta)
            ELSE
               D2 = DENOM_DD2(z, zp, g, gp, v, u, m1, m2, mu, Mtot, zeta)
            END IF
            IF (.NOT. (D2 .GT. 0.d0)) ok = .FALSE.
            N  = NUM_P2(kmt, c, z, zp, g, v, D2, Mtot)
            soma_u = soma_u + dT(s)*c*N/D2**3
         end do
         soma = soma + dW(r)*v**2*(1.d0 - v)*soma_u
      end do

      valor = -0.5d0*soma/D0
      IF (.NOT. (ABS(valor) .LE. HUGE(valor))) ok = .FALSE.
      END SUBROUTINE P2_REG


! -----------------------------------------------------------------------
! P2_DELTA — coeficiente do termo local: P2_delta = delta(z'-z) * valor
!   valor = 2 Int_0^1 dv v^3 / ( a (a + (1-zeta)(1-v) mu^2) )
!   a(v)  = -v(1-v) gamma - v (kappa^2 + gamma' + z^2 M^2/4) - (1-v) mu^2
! (a integral em u ja foi feita analiticamente). ok = .FALSE. se a >= 0,
! a + (1-zeta)(1-v)mu^2 >= 0 ou valor nao finito.
! -----------------------------------------------------------------------
      PURE SUBROUTINE P2_DELTA(g, z, gp, m1, m2, mu, Mtot, zeta, W, dW, Nv, valor, ok)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Nv
      DOUBLE PRECISION, INTENT(IN) :: g, z, gp, m1, m2, mu, Mtot, zeta
      DOUBLE PRECISION, INTENT(IN) :: W(Nv), dW(Nv)
      DOUBLE PRECISION, INTENT(OUT) :: valor
      LOGICAL, INTENT(OUT) :: ok
      INTEGER :: r
      DOUBLE PRECISION :: k2, v, a, b, soma

      ok = .TRUE.
      k2 = (0.5d0*(m1 + m2))**2 - 0.25d0*Mtot**2
      soma = 0.d0
      do r = 1, Nv
         v = W(r)
         a = -v*(1.d0 - v)*g - v*(k2 + gp + 0.25d0*z**2*Mtot**2) - (1.d0 - v)*mu**2
         b = a + (1.d0 - zeta)*(1.d0 - v)*mu**2
         IF (.NOT. (a .LT. 0.d0 .AND. b .LT. 0.d0)) ok = .FALSE.
         soma = soma + dW(r)*v**3/(a*b)
      end do

      valor = 2.d0*soma
      IF (.NOT. (ABS(valor) .LE. HUGE(valor))) ok = .FALSE.
      END SUBROUTINE P2_DELTA

      END MODULE KERNEL_P2_ZETA
