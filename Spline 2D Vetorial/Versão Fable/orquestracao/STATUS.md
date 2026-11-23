# STATUS — orquestração IC kernel vetorial

Pasta de trabalho: `Versão Fable/` · Referência somente leitura: `../Versão Original/` (intacta)
Regra alterada pelo Victor em 16/11/2026: executar o plano inteiro sem pedir aprovação por fase. Continuam valendo: nada apagado/sobrescrito/renomeado fora do que foi gerado aqui, nada fora da pasta de trabalho, nenhum job submetido.

| tarefa | estado | branch/commit | próximo passo |
|---|---|---|---|
| Fase 0 — mapa | concluída | `fable/fase0` 7cb9351 | — |
| Fase 0 — baseline | concluída (`baseline.md`) | `fable/fase0` b95e405 | — |
| T1 — γ_min(z) | concluída, testes ok (`relatorios/T1_gamma_min.md`) | `fable/t1-gamma-min` 4fcb8f0 | confirmar convenção de sinal com o orientador; merge em `main` fica com o Victor |
| T5 — saídas | concluída (`relatorios/T5_saidas.md`) | `tools/run.py` em `fable/fase0` | — |
| T4 — .agr | concluída (`relatorios/T4_agr.md`) | `tools/agr.py`, `plot_run.py`, `plot_compare.py` | — |
| T2 — convergência local | concluída (`relatorios/T2_convergencia.md`): 91 rodadas; λ **não** converge em γ_max, `e` e `z_borda` | `fable/t2-t3` | levar ao orientador os três pontos de não convergência |
| T3 — plano do cluster | concluída (`cluster/plano.md`, `manifest.csv` com 110 jobs, ~178 CPU·h, dry-run ok); **nada submetido** | `fable/t2-t3` | Victor informar cluster (partição, tempo, módulos) e ajustar `job_array.slurm` |
