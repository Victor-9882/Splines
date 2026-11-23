# Relatório final — orquestração IC kernel vetorial (16/11/2026)

## O que mudou
- `src/gamma_min.f90`: γ_min(z) = Δ² − 2 z m̄ Δ, Δ = (m2−m1)/2 com sinal; usado pelo programa principal e pelo de plot via γ̃ = γ' − γ_min(z'). Chaves novas em `&teste`: `USAR_GMIN` (padrão `.TRUE.`) e `CHECAR_DENOM`.
- Limite inferior de γ' de volta a `0.d0` na cópia (o `0.22d0` não commitado deixava os nós fora de ordem).
- `tools/run.py` (árvore `Testes/<teste>/<run_id>/{input,raw,derived,plots,logs}` + `resumo.csv`), `tools/agr.py`, `plot_run.py`, `plot_compare.py`, `conv_local.py`.
- `cluster/`: `plano.md`, `manifest.csv` (110 jobs, ~178 CPU·h), `job_array.slurm`, `collect.py`, `make_manifest.py`. Dry-run local ok. **Nada submetido.**
- Git: branches `fable/fase0`, `fable/t1-gamma-min`, `fable/t2-t3` (4 commits); sem merge em `main`, sem push. `Versão Original/` intacta.

## Resultados-chave
- Baseline (m1=m2=1, Mtot=1.99, μ=0.5, ζ=1, 16×16): λ = 0.0771650238, ⟨ξ⟩ = 0.5000. Réplica 20×20 reproduz a rodada de 10/nov (2e-15).
- T1: Δ=0 idêntico byte a byte; troca m1↔m2 preserva λ (4.9e-13) e leva g(γ,z) em g(γ,−z); nenhum denominador não positivo.
- Káon, ζ=1: λ 0.030275 → 0.030394 (+0.39%) e ⟨ξ⟩ 0.2744 → 0.2839 com γ_min.
- Pico espúrio: γ_min corrige a borda do lado leve (g(0,z=+1)/max de 1.00 para 0.004); os picos de u(ξ) perto das bordas continuam e mudam de posição com a malha.
- T2 (91 rodadas): λ converge em NMG, Ng, Nu; não em NMZ (4e-4 por passo) nem em Nv (1e-3 a 2e-3).
- **λ não converge em três parâmetros:** γ_max (káon: 7–11% por passo, logarítmico; iguais: 1–2%), regulador `e` (káon ~8% por duas décadas; iguais ~1%; sem patamar), `z_borda` (káon entre 0.028 e 0.083; iguais entre 0.0708 e 0.0772).
- ζ de 0 a 3: λ contínuo e monótono; ζ=1 reproduz o código anterior.
- Tempo: t ≈ 1.19e-8 s · (NMG·NMZ)²·Ng·Nz·Nv com 8 threads.

## Dúvidas para o orientador
1. Sinal de γ_min: fixado pela estrutura de D0 e por teste numérico (com o sinal oposto g não se anula em z=+1). O teste de troca de massas não distingue os dois sinais.
2. A dependência em γ_max é um corte ultravioleta genuíno do kernel vetorial na escada?
3. Como tratar a falta de patamar em `e` e a sensibilidade a `z_borda` (hipótese: X mal condicionada; linhas de colocação dos nós extremos não resolvidas pela quadratura).

## Pendências do Victor
- Informar o cluster (partição, tempo, módulos) e ajustar os campos `AJUSTAR` de `cluster/job_array.slurm`.
- Decidir o merge das branches `fable/*`.
- Sugestão não implementada: tirar o kernel do laço (k,l) em `BUILD_ZMATRIX` (ganho de ~100× no tempo).

Detalhes: `relatorios/T1_gamma_min.md`, `relatorios/T2_convergencia.md`, `cluster/plano.md`, `baseline.md`.
