# Baseline (Fase 0) — referência de regressão

Código: cópia de `Versão Original` (árvore de trabalho de 16/11/2026) com o limite inferior de γ' de volta a `0.d0`
(o `0.22d0` não commitado deixava os nós fora de ordem: gv(2)=0.0704 < gv(1)=0.22). Commit `b95e405`, gfortran 13.3 `-O3 -fopenmp -m64` + MKL 2026.1, 8 threads.
Física: m1 = m2 = 1.0, Mtot = 1.99, μ = 0.5, ζ = 1 (Feynman), kernel vetorial; γ_max = 3, A_z = A_g = 1, z_borda = 0.9999999, e = 1e-6.

| | B0 (referência) | B1 (réplica de 10/nov) |
|---|---|---|
| malha | NMG=NMZ=16, Ng=Nz=32, Nv=24 | NMG=NMZ=20, Ng=Nz=40, Nv=30 |
| λ = wr_max | 0.07716502377625574 | 0.07723524144479198 |
| Alfa = 1/(16π λ) | 0.2578158719184969 | 0.2575814811261706 |
| ⟨ξ⟩ | 0.5000000000268 | 0.4999999999867 |
| ∫u dξ | 1.0000219 | — |
| max\|u(ξ)−u(1−ξ)\|/max u | 7.2e-11 | — |
| autovalores complexos | 104 de 256 | — |
| tempo | 16.6 s | 122 s (máquina ocupada; 85 s em 10/nov) |

- B1 reproduz a rodada `Resultados/2026-11-10_zeta_Massas_Iguais` (λ = 0.07723524144479216, Alfa = 0.2575814811261701): diferença relativa 2e-15 (MKL × LAPACK do sistema). Logo aquela rodada usou limite inferior 0.
- ⟨ξ⟩ não é calculado pelo Fortran: é ∫ξ u dξ / ∫u dξ por trapézio sobre `derived/dist_long.dat` (401 pontos), feito em `tools/run.py`.
- Dados: `Testes/baseline/<run_id>/` (B0: `m1-1_m2-1_M-1.99_mu-0.5_zeta-1_G16x16_q32-32-24-24_gmax-3`).
