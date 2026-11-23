# T5 — Organização das saídas

O Fortran já prefixava tudo com `&saida pasta`; não foi preciso mexer nos `open`. A organização é feita por `tools/run.py`:
roda os dois programas dentro de `<run_dir>/work`, depois move os arquivos recém-gerados para `raw/` e `derived/` com nomes canônicos (sem o prefixo `<nome>_zeta<v>_`).

```
Testes/<teste>/<run_id>/
  input/params.nml  input/run_info.txt
  raw/      coeficientes.dat autovalores.dat alfa.dat alfas.dat parametros.dat erros.dat resumo.dat inputs.dat plot_parametros.dat [denominadores.dat termos.dat matrizes.bin]
  derived/  dist_long.dat dist_trans.dat psi.dat psi_z.dat psi2_3d.dat autovetoresG.dat autovetoresZ.dat observaveis.txt
  plots/agr/ plots/png/ logs/spline.log logs/plot.log
Testes/<teste>/resumo.csv
```
- run_id: `m1-<>_m2-<>_M-<>_mu-<>_zeta-<>_G<NMG>x<NMZ>_q<Ng>-<Nz>-<Nv>-<Nu>_gmax-<>` + sufixos só quando diferem do padrão (`Ag-`, `Az-`, `zb-`, `e-`, `gmin0`, `--tag`).
- `resumo.csv`: parâmetros + lambda, alfa, xi_medio, norma_u, assim_u, psi_z_borda_menos/mais, i_max, n_complexos, lambda_2, norma_psi, tempo_s, tempo_parede_s, threads, data. Escrita com trava (`flock`).
- `run_info.txt`: data, hostname, commit/branch (marca `src modificado`), flags de compilação (`bin/flags.txt`, gravado pelo Makefile), threads, tempo de parede, códigos de retorno, job SLURM.
- Nunca sobrescreve: run_id existente é recusado (`--tag` cria outro; `--skip-existing` pula). `--dry-run` não cria nada.
- Uso: `make` (com MKL carregado) e `python3 tools/run.py --teste <nome> --set NMG=20 zeta=0.5 [--plots]`.
- Nada de `Versão Original/` (nem `Resultados/`, `Testes/` antigos) foi movido, renomeado ou apagado.
