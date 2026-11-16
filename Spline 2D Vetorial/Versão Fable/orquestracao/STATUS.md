# STATUS — orquestração IC kernel vetorial

Pasta de trabalho: `Versão Fable/` · Fontes de referência (somente leitura): `../Versão Original/`
Fontes atuais: `2D_Spline_Vetorial_Diferentes_zeta_omp.f90`, `kernel.f90`, `Plot_funcoes_zeta.f90`, `inputs.dat`

| tarefa | estado | branch/commit | próximo passo |
|---|---|---|---|
| Fase 0 — mapa do código | concluída (16/11) | — | `MAPA_CODIGO.md` + `relatorios/T0_mapa.md`; revisado pelo orquestrador (seção 8) |
| Fase 0 — baseline | proposta enviada; aguardando aprovação | — | decidir limite inferior de γ' (0.d0 do HEAD vs 0.22 não commitado) e malha |
| T1 — γ_min(z) | não iniciada | — | plano no checkpoint da Fase 0 |
| T5 — organização de saídas | não iniciada | — | plano no checkpoint da Fase 0 |
| T4 — .agr automático | não iniciada | — | depende de T5 |
| T2 — convergência local | não iniciada | — | depende de T1, T4, T5 |
| T3 — plano do cluster | não iniciada | — | depende de T2 |
