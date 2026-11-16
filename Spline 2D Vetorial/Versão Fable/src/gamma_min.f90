!=======================================================================
! gamma_min.f90
!
! Limitante inferior do suporte de g(gamma', z') na representacao de
! Nakanishi com o denominador escrito em termos de kappa^2:
!
!     gamma' >= gamma_min(z') = Delta**2 - 2 z' mbar Delta
!
! com mbar = (m1 + m2)/2 e Delta = (m2 - m1)/2 COM SINAL (convencao do
! codigo: kernel.f90, D0 = gamma + (1-z**2) kappa**2 + (Delta - z mbar)**2,
! ou seja, z = +1 <-> m1 e z = -1 <-> m2).
!
! Origem: o denominador "natural" de Nakanishi para massas distintas e
!     gamma'_nat + m1**2 (1+z')/2 + m2**2 (1-z')/2 - M**2/4 - k**2 - p.k z'
! com gamma'_nat >= 0 (ausencia de cortes). Escrevendo-o com kappa**2 =
! mbar**2 - M**2/4, gamma' = gamma'_nat + gamma_min(z'). Identidades:
!     gamma_min(z) + z**2 mbar**2 = (z mbar - Delta)**2 >= 0
!     gamma_min(+1) = m1**2 - mbar**2,  gamma_min(-1) = m2**2 - mbar**2
!     D0(gamma, z) = gamma + gamma_min(z) + (1-z**2) kappa**2 + z**2 mbar**2
! Para m2 > m1 coincide com a forma do artigo boson-fermion (eqs. 15-16),
! -2 z mbar |Delta| + Delta**2. Delta = 0 => gamma_min = 0 exatamente.
!
! Os programas trabalham com a variavel deslocada
!     gamma~ = gamma' - gamma_min(z')  em [0, gamma_max],
! na qual g e expandida em splines; gamma' = gamma~ + gamma_min(z') entra
! onde gamma' aparece explicitamente (kernel e denominador do lado esquerdo).
!
! Este arquivo e incluido (INCLUDE) pelo programa principal e pelo
! programa de plot: e a UNICA definicao de gamma_min.
!=======================================================================
      MODULE GAMMA_MIN_MOD
      IMPLICIT NONE
      CONTAINS

! -----------------------------------------------------------------------
! gamma_min(z) = Delta**2 - 2 z mbar Delta
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION GAMMA_MIN(z, mbar, Delta)
      DOUBLE PRECISION, INTENT(IN) :: z, mbar, Delta

      GAMMA_MIN = Delta**2 - 2.d0*z*mbar*Delta
      END FUNCTION GAMMA_MIN

! -----------------------------------------------------------------------
! Mesma funcao a partir de m1, m2, com chave liga/desliga:
! ativo = .FALSE. devolve 0 (codigo anterior, gamma' a partir de 0).
! -----------------------------------------------------------------------
      PURE DOUBLE PRECISION FUNCTION GAMMA_MIN_EF(z, m1, m2, ativo)
      DOUBLE PRECISION, INTENT(IN) :: z, m1, m2
      LOGICAL, INTENT(IN) :: ativo

      IF (ativo) THEN
         GAMMA_MIN_EF = GAMMA_MIN(z, 0.5d0*(m1 + m2), 0.5d0*(m2 - m1))
      ELSE
         GAMMA_MIN_EF = 0.d0
      END IF
      END FUNCTION GAMMA_MIN_EF

      END MODULE GAMMA_MIN_MOD
