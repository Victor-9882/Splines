! =======================================================================
!            TERMO (1-zeta) DO PROPAGADOR VETORIAL: P^(2)
!
! Incluido (INCLUDE) no inicio de 2D_Spline_Vetorial_Diferentes_zeta_omp.f90
! e no driver de testes (testes_zeta/). Todas as rotinas sao PURE: sem
! SAVE e sem estado global, podem ser chamadas de dentro do laco OpenMP.
! Ficam num MODULE para que as interfaces PURE sejam explicitas.
!
! Lado direito da EBS:  -lambda^2/(32 pi^2) Int dz' Int dgamma' (P1 - (1-zeta) P2) g
!
!   P2 = P2_reg + P2_delta
!
!   P2_reg   = -1/2 Int_0^1 dv v^2 Int_0^{1-v} du
!                [ N_ku (1+z)^3 th(z'-z) / (D0 Du2^3) + N_kd (1-z)^3 th(z-z') / (D0 Dd2^3) ]
!            = +1/2 Int dv v^2 Int du [ N_ku/(D0 f_u^3) th(z'-z) + N_kd/(D0 f_d^3) th(z-z') ]
!              (Du2 = -(1+z) f_u,  Dd2 = -(1-z) f_d)
!
!   P2_delta = 2 delta(z'-z) Int_0^1 dv v^3 / ( a(v) [a(v) + (1-zeta)(1-v) mu^2] )
!
! Convencoes: Mtot -> M, g -> gamma, gp -> gamma', zp -> z'.
! =======================================================================
      MODULE KERNEL_P2_ZETA
      IMPLICIT NONE
      CONTAINS

! -----------------------------------------------------------------------
! Polos dos propagadores escalares em k^-
!   k_u = M/2 - 2(gamma+m2^2)/(M(1+z)),   k_d = -M/2 + 2(gamma+m1^2)/(M(1-z))
! -----------------------------------------------------------------------
      PURE SUBROUTINE POLOS_KMENOS(z, g, m1, m2, Mtot, ku, kd)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: z, g, m1, m2, Mtot
      DOUBLE PRECISION, INTENT(OUT) :: ku, kd

      ku =  0.5d0*Mtot - 2.d0*(g + m2**2)/(Mtot*(1.d0 + z))
      kd = -0.5d0*Mtot + 2.d0*(g + m1**2)/(Mtot*(1.d0 - z))
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
! Numerador N(k^-) do termo (1-zeta), com k2, kp e f no mesmo k^-:
!   k^2 = -z (M/2) k^- - gamma,   k.p = (M/2) k^- - z M^2/4
!   N   = 4 d0 + f ( kp d1 + 4 d2 + M^2 d3 + k2 d4 + 12 f d5 )
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION NUM_P2(km, z, zp, g, v, f, Mtot)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: km, z, zp, g, v, f, Mtot
      DOUBLE PRECISION :: k2, kp, d0, d1, d2, d3, d4, d5, M

      M  = Mtot
      k2 = -z*0.5d0*M*km - g
      kp = 0.5d0*M*km - 0.25d0*z*M**2

      d0 = -(v**2/16.d0) * ( 4.d0*((v-1.d0)*zp - 1.d0)*kp + 4.d0*(v-2.d0)*k2 + M**2*zp*(v*zp - 2.d0) ) &
                         * ( 4.d0*((v-1.d0)*zp + 1.d0)*kp + 4.d0*(v-2.d0)*k2 + M**2*zp*(v*zp + 2.d0) )
      d1 = -4.d0*(v-1.d0)*v*zp
      d2 = -(v/2.d0) * ( zp*(4.d0*(v-1.d0)*kp + M**2*v*zp) + 4.d0*(v-2.d0)*k2 )
      d3 = 1.d0 - v**2*zp**2
      d4 = -4.d0*(v-1.d0)**2
      d5 = -1.d0

      NUM_P2 = 4.d0*d0 + f*( kp*d1 + 4.d0*d2 + M**2*d3 + k2*d4 + 12.d0*f*d5 )
      END FUNCTION NUM_P2


! -----------------------------------------------------------------------
! P2_REG — parte regular de P^(2) num ponto (gamma, z, gamma', z')
!   ramo = +1 : z' em [z, 1]  (theta(z'-z)), polo k_u
!   ramo = -1 : z' em [-1, z] (theta(z-z')), polo k_d
!   valor = 1/(2 D0) Int dv v^2 Int du N(k)/f(k)^3
! v: Gauss (W, dW) em [0,1] (o mesmo de P^(1)); u = (1-v) t, t: Gauss (T, dT).
! ok = .FALSE. se f >= 0 ou algum valor nao finito aparecer (nao mascara).
! Se 1+z = 0 (ramo +1) ou 1-z = 0 (ramo -1) o ramo vale 0, que e o limite
! exato (N/f^3 ~ 1/k^- -> 0); so acontece com z_borda = 1.
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
      DOUBLE PRECISION :: ku, kd, km, D0, v, u, KDc, lD, f, N, soma, soma_u

      valor = 0.d0
      ok = .TRUE.
      IF (ramo .EQ. 1) THEN
         IF (1.d0 + z .EQ. 0.d0) RETURN
      ELSE
         IF (1.d0 - z .EQ. 0.d0) RETURN
      END IF

      CALL POLOS_KMENOS(z, g, m1, m2, Mtot, ku, kd)
      IF (ramo .EQ. 1) THEN
         km = ku
      ELSE
         km = kd
      END IF
      D0 = D0_P2(z, g, m1, m2, Mtot)

      soma = 0.d0
      do r = 1, Nv
         v  = W(r)
         KDc = KD_P2(z, zp, v, Mtot)
         soma_u = 0.d0
         do s = 1, Nu
            u  = (1.d0 - v)*T(s)
            lD = LD_P2(z, zp, g, gp, v, u, m1, m2, mu, Mtot, zeta)
            f  = KDc*km + lD
            IF (.NOT. (f .LT. 0.d0)) ok = .FALSE.
            N  = NUM_P2(km, z, zp, g, v, f, Mtot)
            soma_u = soma_u + dT(s)*N/f**3
         end do
         soma = soma + dW(r)*v**2*(1.d0 - v)*soma_u
      end do

      valor = 0.5d0*soma/D0
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
