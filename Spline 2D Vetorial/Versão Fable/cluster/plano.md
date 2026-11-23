# T3 — Plano de convergência em larga escala (cluster) — NADA FOI SUBMETIDO

## Pendências antes de submeter (decisão do Victor)
- Qual cluster: partição/fila, limite de tempo, núcleos por nó, módulos (gcc/gfortran, MKL, python3). Em `job_array.slurm` estão marcados com `AJUSTAR`.
- Recompilar no cluster: `make` (MKL) ou `make LAPACK=sistema`.

## O que rodar e por quê
A T2 mostrou que o erro dominante não é de malha: λ depende de γ_max, do regulador `e` e de `z_borda` (ver `orquestracao/relatorios/T2_convergencia.md`). O plano prioriza caracterizar esses três em malha melhor, e só então a varredura física. Malha de referência: NMG=24, NMZ=24, Ng=48, Nz=48, Nv=48, Nu=12, γ_max=3.

| grupo (pasta em `Testes/`) | conteúdo | jobs |
|---|---|---|
| `cluster_convergencia` | casos E e K, ζ=1; um de cada vez: NMG {20,24,32,40}, NMZ {20,24,32,40}, Ng {32,48,64,96}, Nz {32,48,64,96}, Nv {32,48,64,96,128} | 34 |
| `cluster_estabilidade` | E e K: γ_max {3,4.5,6,9,12} com NMG {24,36,48,72,96} (densidade constante); e {1e-14…1e-2}; z_borda {0.9…0.9999999} | 34 |
| `cluster_fisica` | ζ {0,0.25,0.5,0.75,1,1.5,2,3} × m2/m1 {1,1.5,2.309} (ligação do káon); Mtot {1,1.367,1.6,1.8,1.9,1.99} e μ {0.15,0.5,1,1.8} com ζ=1 e 0 | 42 |
| **total** | `cluster/manifest.csv` (idx → teste, run_id, parâmetros, estimativas) | **110** |

## Custo (16 núcleos por job; `python3 cluster/make_manifest.py --cpus N` recalcula)
- Tempo por job: t ≈ a·(NMG·NMZ)²·Ng·Nz·Nv·(8/threads), a = 1.19e-8 s ajustado a 79 rodadas locais da T2. Supõe escala linear em threads na montagem (laço OpenMP sobre os pontos de colocação) e mesma velocidade por núcleo da máquina local: tratar como ordem de grandeza (a T2 teve dispersão de 0.8 a 2.0).
- DGGEV: O(n³) com n = NMG·NMZ; desprezível até n ≈ 2300 (estimado em ~25 s).
- Mediana 0.06 h (4 min) por job; máximo 0.97 h (γ_max=12, NMG=96, n=2304). **Total ≈ 178 CPU·h** (convergência 51, estabilidade 86, física 41).
- Memória: ~6 matrizes n×n em dupla precisão: 0.2 GB típico, 0.44 GB no maior job. `--mem=4G` sobra.
- Sugestão de diretivas: `--cpus-per-task=16`, `--time=04:00:00`, `--array=1-110%20`.

## Arquivos
- `cluster/make_manifest.py` — gera `manifest.csv` (e `manifest_dryrun.csv` com `--dryrun`).
- `cluster/job_array.slurm` — uma linha do manifest por tarefa; `OMP_NUM_THREADS = MKL_NUM_THREADS = $SLURM_CPUS_PER_TASK`; confere que o run_id do manifest é o gerado por `tools/run.py`; grava só em `Testes/<teste>/<run_id>/` com `--no-csv --skip-existing` (sem arquivo compartilhado entre jobs; relançar o array só refaz o que falta).
- `cluster/collect.py Testes/<teste> --plots --png` — reconstrói `resumo.csv` a partir dos run_id e gera λ×N por parâmetro com `tools/plot_compare.py`.

## Dry-run local (feito em 16/11/2026)
`python3 cluster/make_manifest.py --dryrun` e, para i = 1, 2: `SLURM_ARRAY_TASK_ID=i SLURM_CPUS_PER_TASK=2 MANIFEST=cluster/manifest_dryrun.csv bash cluster/job_array.slurm`.
As duas entradas (malha 8×8, E com ζ=1 e K com ζ=0.5) rodaram em 0.6 s cada, geraram a árvore completa em `Testes/cluster_dryrun/` e foram agregadas por `collect.py` (2 run_id). Os 110 run_id do manifest definitivo foram conferidos contra `run.py --print-run-id` (0 divergências); índice inexistente é recusado com erro.

## Ordem sugerida
1. `cluster_estabilidade` primeiro (86 CPU·h): decide como tratar γ_max, e, z_borda.
2. Com isso resolvido, `cluster_convergencia` e `cluster_fisica`.
3. Se o kernel for tirado do laço (k,l) (ver T2, seção 4), o custo total cai para poucas CPU·h e as malhas podem ser bem mais finas.
