# T1 — Limitante inferior γ_min(z)

## Decisão de convenção (documentada)
γ_min(z) = Δ² − 2 z m̄ Δ, com **Δ = (m2 − m1)/2 com sinal** e m̄ = (m1+m2)/2 (`src/gamma_min.f90`, única definição, incluída pelos dois programas).
- Motivo: no código D0 = γ − (1−z²)M²/4 + m1²(1+z)/2 + m2²(1−z)/2 = γ + (1−z²)κ² + (Δ − z m̄)² (kernel.f90:43, 279), ou seja z=+1 ↔ m1 e z=−1 ↔ m2. O denominador natural de Nakanishi é γ'_nat + m1²(1+z')/2 + m2²(1−z')/2 − M²/4 − k² − p·k z' com γ'_nat ≥ 0; escrito com κ² dá γ' = γ'_nat + γ_min(z').
- γ_min(+1) = m1² − m̄², γ_min(−1) = m2² − m̄²; γ_min(z) + z²m̄² = (z m̄ − Δ)² ≥ 0; D0 = γ + γ_min(z) + (1−z²)κ² + z²m̄².
- Para m2 > m1 (todos os casos usados até aqui) coincide com a forma do artigo, −2 z m̄ |Δ| + Δ². Para m1 > m2 a forma com |Δ| quebraria a simetria m1↔m2; a com sinal a preserva (teste abaixo).
- **Limite do teste de troca:** ele distingue "com sinal" de "|Δ|", mas não o sinal global do termo em z (Δ² + 2 z m̄ Δ também seria simétrico e também daria denominadores positivos). Esse sinal foi fixado pela estrutura de D0 acima. **Ponto para confirmar com o orientador.**

## Implementação (variável deslocada γ̃ = γ' − γ_min(z') ∈ [0, γ_max])
Base de splines, nós, colocação e quadratura em γ ficam iguais; `kernel.f90` não mudou (recebe γ' físico pelo argumento `gp`).
Nova chave no namelist `&teste`: `USAR_GMIN` (padrão `.TRUE.`; `.FALSE.` = código anterior) e `CHECAR_DENOM` (grava `denominadores.dat`).
Substituições em `2D_Spline_Vetorial_Diferentes_zeta_omp.f90` (M) e `Plot_funcoes_zeta.f90` (P):
- M `BUILD_ZMATRIX`: `gmz` (M:874); `gpf = gp + GAMMA_MIN_EF(zq,...)` nos dois ramos (M:906, 931); 6 chamadas do kernel com `gpf` (M:916-917, 922, 938-939, 944); lado esquerdo `g + gp + gmz + (1-z²)κ² + m̄²z²` (M:956); `gpf, gmz` no PRIVATE do OpenMP (M:885).
- M `BUILD_ZMATRIX_P2` (termo 1−ζ): `gmz` (M:1228), `gmq` nos dois ramos (M:1237, 1241), `P2_REG(..., Y(p) + gmq, ...)` (M:1247), `P2_DELTA(..., Y(p) + gmz, ...)` (M:1271), PRIVATE (M:1221).
- M diagnósticos `TERMOS_P1` (M:1040-1144) e `TERMOS_P2` (M:1345-1382): mesma substituição.
- M nova rotina `CHECA_DENOMINADORES` (M:~705-775); eco de Δ e γ_min(±1) em `parametros.dat` (M:300-301).
- P: `+ GAMMA_MIN_EF(z, m1, m2, usar_gmin)` nos 7 denominadores do propagador (P:362, 397, 449, 534, 581, 639, 696): ψ, norma, u(ξ), D⊥, f1. `usar_gmin = .FALSE.` antes de ler `parametros.dat` (P:198), para que resultados antigos sejam plotados como foram calculados.
- Nos arquivos `autovetoresG/Z` a abscissa γ é **γ̃** (γ' = γ̃ + γ_min(z)); nos de ψ, u, D⊥ o γ externo é o físico.

## Testes (malha B0: 16×16, Ng=Nz=32, Nv=Nu=24; dados em `Testes/T1_gamma_min/`)
Káon: m1 = 0.6044, m2 = 1.3956, Mtot = 1.3671, μ = 1.8 (γ_min ∈ [−0.635, +0.948]). Leve: m1 = 1.0, m2 = 1.1, Mtot = 2.04, μ = 0.5.

| teste | critério | resultado | |
|---|---|---|---|
| Δ=0, ζ=1: novo × baseline | \|Δλ/λ\| < 1e-12 | λ = 0.07716502377625574 nos dois; `coeficientes.dat` e `dist_long.dat` idênticos byte a byte | passou |
| Δ=0, ζ=0.5: novo × binário pré-T1 | idem | λ = 0.06629201497442426 nos dois | passou |
| `USAR_GMIN=F` × binário pré-T1 (káon, ζ=1 e ζ=0.5) | idêntico | λ = 0.03027549824636168 e 0.02445733876054305 nos dois | passou |
| denominadores > 0 (káon, leve, ζ=1 e 0.5) | 0 não positivos | káon: mín D0 = 0.357, LHS = 0.361, Du/(1+z) = 0.365, Dd/(1−z) = 0.367; leve: 0.063, 0.067, 0.067, 0.067; 6.3e6 pontos por ramo; P² (ζ=0.5) sem aborto | passou |
| troca m1↔m2 (káon) | λ invariante, g(γ,z)↔g(γ,−z) | \|Δλ/λ\| = 4.9e-13; max\|Δg\|/max\|g\| = 8.3e-9; ⟨ξ⟩+⟨ξ⟩' − 1 = −1.6e-8 | passou |
| troca m1↔m2 (leve) | idem | 4.9e-14; 2.6e-10; 7.0e-12 | passou |

Efeito de γ_min (ζ=1, mesma malha):

| caso | λ sem | λ com | Δλ/λ | ⟨ξ⟩ sem → com |
|---|---|---|---|---|
| káon ζ=1 | 0.0302754982 | 0.0303938301 | +3.9e-3 | 0.27443 → 0.28391 |
| káon ζ=0.5 | 0.0244573388 | 0.0245504882 | +3.8e-3 | 0.28132 → 0.28521 |
| leve ζ=1 | 0.0637989862 | 0.0637943387 | −7.3e-5 | 0.47438 → 0.47466 |

A diferença é esperada: sem γ_min, g é forçada a ter suporte em γ' ≥ 0 para todo z'. No lado leve (z'→+1, γ_min < 0) isso corta a parte γ' ∈ [γ_min, 0) do suporte verdadeiro; no lado pesado (γ_min > 0) a spline tenta representar um trecho onde g deveria ser nula.

## Pico espúrio
- **Na borda do lado leve o defeito some:** sem γ_min, g(γ̃=0, z) cresce até o seu máximo em z=+1 (g(0,+1)/max = 1.00) e ψ(γ=2, z=+1)/max = 0.114; com γ_min, 0.004 e 0.025. Com as massas trocadas o mesmo ocorre em z=−1.
- **Os picos de u(ξ) perto das bordas permanecem:** káon, máximos locais em ξ ≈ 0.05, 0.19, 0.79, 0.95 com ou sem γ_min (o de ξ≈0.95, z≈−0.9, vai de 1.07 para 1.19). A rodada antiga `2026-11-10_zeta_varredura_kaon` (20×20) tem os mesmos picos em ξ = 0.040, 0.150, 0.845, 0.960: a posição muda com a malha, o que indica artefato de discretização/corte e não de γ_min. Ver T2 (γ_max e NMZ).
- Figuras: `Testes/T1_gamma_min/comparacoes/{psi_z,g_z,u}_kaon_gmin.{agr,png}`.
