# Versão Fable: kernel vetorial, spline 2D

Resolve o problema de autovalores generalizado `Z v = λ X v` para a função de Nakanishi `g(γ, z)` do estado ligado de dois bósons de massas `m1`, `m2` (kernel vetorial `P^(1) − (1−ζ) P^(2)`, ζ = parâmetro de calibre). A constante de acoplamento sai de `Alfa = 1/(16 π λ_max)`.

Esta pasta é a cópia de trabalho. `../Versão Original/` é a referência e não é alterada. O histórico do que mudou e por quê está em [orquestracao/](orquestracao/) (`STATUS.md`, `RELATORIO_FINAL.md`, `MAPA_CODIGO.md`, `relatorios/`).

Em volta dos dois programas Fortran há um conjunto de scripts Python que automatiza o ciclo inteiro: **rodar → organizar saídas → extrair observáveis → tabelar (CSV) → gerar `.agr` → converter para PNG**, de um run isolado até varreduras de centenas de jobs no cluster.

---

## Sumário

1. [Mapa das ferramentas](#1-mapa-das-ferramentas)
2. [Requisitos](#2-requisitos)
3. [Compilar (`Makefile`)](#3-compilar-makefile)
4. [Receitas rápidas](#4-receitas-rápidas)
5. [Rodar um teste (`tools/run.py`)](#5-rodar-um-teste-toolsrunpy)
6. [Saídas, `observaveis.txt` e `resumo.csv`](#6-saídas-observaveistxt-e-resumocsv)
7. [Levas e varreduras (`tools/leva.py`, laço, `tools/conv_local.py`)](#7-levas-e-varreduras)
8. [Gráficos Xmgrace (`agr.py`, `plot_run.py`, `plot_compare.py`)](#8-gráficos-xmgrace)
9. [Conversão `.agr` → PNG / EPS / SVG](#9-conversão-agr--png--eps--svg)
10. [Cluster (`cluster/`)](#10-cluster)
11. [Modo direto: vários testes num só `inputs.dat`](#11-modo-direto-vários-testes-num-só-inputsdat)
12. [Parâmetros de `&teste` e `&plot`](#12-parâmetros-de-teste-e-plot)
13. [Resultados conhecidos e pontos em aberto](#13-resultados-conhecidos-e-pontos-em-aberto)

---

## 1. Mapa das ferramentas

| arquivo | o que faz | entrada | saída |
|---|---|---|---|
| `Makefile` | compila os dois programas e registra as flags | `src/*.f90` | `bin/spline_omp`, `bin/plot`, `bin/flags.txt` |
| `src/2D_Spline_Vetorial_Diferentes_zeta_omp.f90` | programa principal: monta `Z` e `X`, roda DGGEV, grava coeficientes e α | namelist (`&saida`, `&plot`, `&teste`…) | `*_coeficientes.dat`, `*_autovalores.dat`, `alfas.dat`, `resumo.dat`… |
| `src/Plot_funcoes_zeta.f90` | pós-processamento: u(ξ), ψ, D⊥, autovetores em γ e z | pasta com `resumo.dat` | `*_plot_*.dat`, `autovetoresG/Z.dat`, norma em `resumo.dat` |
| `src/kernel.f90`, `src/gamma_min.f90` | kernels P^(1), P^(2) e limitante γ_min(z) (incluídos pelos dois programas) | — | — |
| **`tools/leva.py`** | roda uma **leva inteira descrita num arquivo `.toml`** (grade de parâmetros, destino, gráficos, comparações, tabela) | `levas/*.toml` | runs + `<pasta>/leva_<nome>/` |
| `levas/exemplo.toml` | modelo comentado de arquivo de leva | — | — |
| **`tools/run.py`** | roda **um** teste, organiza a árvore, calcula observáveis, alimenta o CSV | `--teste`, `--set CHAVE=VALOR` | `Testes/<teste>/<run_id>/`, `resumo.csv` |
| **`tools/conv_local.py`** | varreduras de convergência prontas (um parâmetro por vez), tabelas e gráficos | — | `Testes/convergencia_inicial/` |
| **`tools/agr.py`** | biblioteca + CLI que escreve `.agr` no estilo de `template_plot_long.agr` | arquivos de 2 colunas | `.agr` |
| **`tools/plot_run.py`** | todos os gráficos de um run (u, g, ψ, D⊥) | `<run_dir>` | `plots/agr/*.agr`, `plots/png/*.png` |
| **`tools/plot_compare.py`** | `overlay` (curvas de vários runs) e `conv` (observável × parâmetro, erro em log) | runs ou `resumo.csv` | `.agr` (+ `.png`) |
| `tools/template_plot_long.agr` | modelo de estilo (fontes, cores, eixos) que o `agr.py` reproduz | — | — |
| **`cluster/make_manifest.py`** | gera a lista de jobs com estimativa de tempo/memória | — | `cluster/manifest.csv` |
| **`cluster/job_array.slurm`** | job array SLURM: uma linha do manifest por tarefa | `manifest.csv` | `Testes/cluster_*/<run_id>/` |
| **`cluster/collect.py`** | reconstrói o `resumo.csv` a partir dos runs e plota λ × N | `Testes/<teste>` | `resumo.csv`, `comparacoes/` |

Todos os scripts usam **só a biblioteca padrão do Python 3** (não precisa de numpy nem matplotlib). Os gráficos são sempre arquivos Xmgrace; o PNG é feito pelo `gracebat`.

## 2. Requisitos

| item | para quê | observação |
|---|---|---|
| `gfortran` | compilar | flags `-O3 -fopenmp -m64` |
| Intel oneAPI MKL (`/opt/intel/oneapi/setvars.sh`) | LAPACK/BLAS (padrão) | alternativa: `make LAPACK=sistema` |
| `python3` | todas as ferramentas | só biblioteca padrão |
| Xmgrace (`xmgrace`, `gracebat`) | abrir `.agr` / converter para PNG, EPS, SVG | só é preciso para ver ou converter os gráficos, não para gerar os `.agr` |

## 3. Compilar (`Makefile`)

```bash
cd "Spline 2D Vetorial/Versão Fable"
source /opt/intel/oneapi/setvars.sh     # MKL (padrão)
make                                    # gera bin/spline_omp e bin/plot
```

| comando | efeito |
|---|---|
| `make` | compila com MKL (precisa de `$MKLROOT`) |
| `make LAPACK=sistema` | usa `-llapack -lblas` do sistema |
| `make clean` | apaga binários, `.mod` e `bin/flags.txt` |

O `Makefile` grava a linha de compilação em `bin/flags.txt`; o `run.py` copia essa linha para o `run_info.txt` de cada rodada.

## 4. Receitas rápidas

```bash
# leva completa a partir de um arquivo (rodadas + gráficos + comparações + tabela)
cp levas/exemplo.toml levas/minha_leva.toml        # editar
python3 tools/leva.py levas/minha_leva.toml --listar
python3 tools/leva.py levas/minha_leva.toml

# um run com os parâmetros padrão (B0), já com gráficos .agr + .png
python3 tools/run.py --teste baseline --plots

# mudar parâmetros
python3 tools/run.py --teste conv --set NMG=20 NMZ=20 zeta=0.5 --plots

# conferir o que seria feito, sem criar nada
python3 tools/run.py --teste conv --set m2=1.5 --dry-run

# varredura em laço (rápida: sem pós-processamento)
for n in 16 20 24; do python3 tools/run.py --teste varre --set NMZ=$n --no-plot --skip-existing; done

# λ × NMZ e |Δλ| em log, direto do CSV, com PNG
python3 tools/plot_compare.py conv --csv Testes/varre/resumo.csv --param NMZ --y lambda --out Testes/varre/lambda.agr --png
python3 tools/plot_compare.py conv --csv Testes/varre/resumo.csv --param NMZ --y lambda --logerr --out Testes/varre/lambda_err.agr --png

# sobrepor u(ξ) de dois runs, legenda pelo ζ, com ⟨ξ⟩ marcado
python3 tools/plot_compare.py overlay --runs Testes/conv/<run1> Testes/conv/<run2> \
        --quantity u --legend-key zeta --vlines --out Testes/conv/u_zeta.agr --png

# converter um .agr qualquer para PNG
gracebat -hardcopy -hdevice PNG -printfile figura.png figura.agr
```

## 5. Rodar um teste (`tools/run.py`)

Executa `bin/spline_omp` e depois `bin/plot` para **uma** combinação de parâmetros (um `run_id`), move as saídas para uma árvore padronizada, calcula os observáveis e acrescenta uma linha ao `resumo.csv` do teste.

```bash
python3 tools/run.py --teste <nome> [--set CHAVE=VALOR ...] [opções]
```

**Parâmetros padrão (B0):** massas iguais (`m1=m2=1`, `Mtot=1.99`, `mu=0.5`), `zeta=1`, `NMG=NMZ=16`, `Ng=Nz=32`, `Nv=Nu=24`, `gamma_max=3`, `A_g=A_z=1`, `z_borda=0.9999999`, `e=1e-6`. Qualquer chave de `&teste` ou `&plot` pode ser trocada com `--set` (a lista completa está na §12). Valores lógicos aceitam `T/F`, `true/false`, `1/0`, `sim/não`; números aceitam notação `1d-6`.

| opção | efeito |
|---|---|
| `--teste NOME` | **obrigatória**: subpasta de `Testes/` onde o run é guardado |
| `--set K=V ...` | altera parâmetros em relação ao padrão (chave desconhecida → erro listando as válidas) |
| `--threads N` | threads OpenMP e MKL (padrão: `$OMP_NUM_THREADS`, senão 8) |
| `--no-plot` | não roda `bin/plot` (sem u, ψ, ⟨ξ⟩; bem mais rápido em varreduras de λ) |
| `--plots` | ao final, roda `tools/plot_run.py <run_dir> --png` (log em `logs/plot_run.log`) |
| `--tag X` | sufixo extra no `run_id` (para repetir uma combinação já existente) |
| `--skip-existing` | se o `run_id` já existe, pula sem erro (permite relançar varreduras interrompidas) |
| `--no-csv` | não escreve no `resumo.csv` (jobs do cluster; ver §10) |
| `--dry-run` | mostra `run_id`, pasta e parâmetros, e sai |
| `--print-run-id` | imprime só o `run_id` e sai (usado pelo `.slurm` para conferir o manifest) |
| `--bin PASTA` | pasta dos executáveis (padrão `bin/`) |
| `--raiz-testes PASTA` | raiz da árvore de saídas (padrão `Testes/`) |

O que o script faz por você:

- **Carrega o MKL sozinho**: se `$MKLROOT` não está definido, faz `source /opt/intel/oneapi/setvars.sh` só para o processo filho.
- **Nunca sobrescreve**: um `run_id` existente é recusado (use `--tag` ou `--skip-existing`).
- **Registra a procedência** em `input/run_info.txt`: data, host, commit e ramo do git (marca `src modificado` se `src/` tiver alterações não commitadas), flags de compilação, threads, tempo de parede, códigos de retorno, comando exato e `SLURM_JOB_ID`.
- **Imprime uma linha de resumo**: `ok | <run_id> | lambda=… alfa=… <xi>=… | 16.6s | N arquivos`. O código de saída é 0 se tudo deu certo e 2 se algum programa falhou.

**O `run_id`** codifica os parâmetros, por exemplo `m1-1_m2-1_M-1.99_mu-0.5_zeta-1_G16x16_q32-32-24-24_gmax-3` (`G` = NMG×NMZ, `q` = Ng-Nz-Nv-Nu). Sufixos `Ag-`, `Az-`, `zb-`, `e-`, `kern-` só aparecem quando esses parâmetros diferem do padrão, `gmin0` aparece se `usar_gmin=F`, e no fim vem a `--tag`, se houver.

## 6. Saídas, `observaveis.txt` e `resumo.csv`

```
Testes/<teste>/
    resumo.csv                       uma linha por run_id (escrita com trava de arquivo)
    <run_id>/
        input/    params.nml         namelist exata usada
                  run_info.txt       procedência (data, host, commit, flags, threads, tempo, rc)
        raw/      coeficientes.dat, autovalores.dat, alfa.dat, parametros.dat, erros.dat,
                  resumo.dat, alfas.dat, plot_parametros.dat
                  (+ termos.dat, matrizes.bin, denominadores.dat quando ligados)
        derived/  dist_long.dat (u(ξ)), dist_trans.dat (D⊥(γ)), psi.dat (ψ(γ)), psi_z.dat (ψ(z)),
                  psi2_3d.dat, autovetoresG.dat, autovetoresZ.dat, observaveis.txt
        plots/    agr/ e png/        (preenchidos por --plots ou tools/plot_run.py)
        logs/     spline.log, plot.log, plot_run.log
```

O `run.py` renomeia os arquivos com prefixo do Fortran (`run_zeta1_plot_dist_long.dat` → `derived/dist_long.dat`, etc.), para que todos os runs tenham os mesmos nomes.

### Observáveis calculados automaticamente

Gravados em `derived/observaveis.txt` (formato `chave = valor`, junto com os parâmetros) e copiados para o `resumo.csv`:

| chave | origem | significado |
|---|---|---|
| `lambda`, `alfa` | `resumo.dat` | maior autovalor real positivo e α = 1/(16πλ) |
| `lambda_2` | `autovalores.dat` | 2º maior autovalor real positivo (separação espectral) |
| `i_max` | `autovalores.dat` | índice do autovalor escolhido |
| `n_complexos` | `autovalores.dat` | nº de autovalores com parte imaginária ≠ 0 |
| `norma_psi` | `resumo.dat` (coluna gravada pelo `bin/plot`) | norma de ψ |
| `norma_u` | `dist_long.dat` (trapézio) | ∫u(ξ)dξ |
| `xi_medio` | `dist_long.dat` | ⟨ξ⟩ = ∫ξu/∫u |
| `assim_u` | `dist_long.dat` | max\|u(ξ) − u(1−ξ)\| / max\|u\| (simetria ξ ↔ 1−ξ) |
| `psi_z_borda_menos/mais` | `psi_z.dat` | ψ(z=∓1)/max\|ψ\| (detecta pico espúrio na borda) |
| `usar_gmin` | `parametros.dat` | `T`/`F` efetivo (`n/a` em binários antigos) |
| `tempo_s` | `resumo.dat` | tempo do cálculo medido pelo Fortran |

Colunas do `resumo.csv`: `run_id, teste`, todos os parâmetros de `&teste`, os observáveis acima, `tempo_parede_s, threads, data`. Ele é a entrada de `plot_compare.py conv`, `conv_local.py --analisar` e `make_manifest.py` (ajuste de tempo). Quem roda `bin/spline_omp` direto (§11) não obtém esse CSV.

## 7. Levas e varreduras

### 7.1 `tools/leva.py`: tudo a partir de um arquivo (recomendado)

Um arquivo de texto (formato TOML: `chave = valor`, `#` comenta) descreve a leva inteira: destino, threads, parâmetros fixos, grade de parâmetros, rodadas avulsas, gráficos por run, comparações, PNG e tabela. Um só comando roda o que falta e gera tudo.

```bash
cp levas/exemplo.toml levas/minha_leva.toml                # modelo comentado; editar
python3 tools/leva.py levas/minha_leva.toml --listar       # plano: rodadas (falta/ok/falhou) e comparações
python3 tools/leva.py levas/minha_leva.toml                # roda o que falta e gera tudo
python3 tools/leva.py levas/minha_leva.toml --so-graficos  # não roda nada; refaz gráficos, CSV e relatório
```

Opções de linha de comando: `--threads N` e `--simultaneos N` sobrepõem os valores do arquivo.

**No terminal**, cada rodada abre com `[#i/N] início <run_id> (hora)`. Em seguida aparece, ao vivo, o que os programas escrevem em `<run>/logs/spline.log` e `plot.log`: progresso da montagem com tempo decorrido e restante, tempos de montagem e do DGGEV, `Alfa = … (T s)`. A rodada fecha com `[#i/N] fim …`, uma linha com λ, α, ⟨ξ⟩ e o tempo total. Com `simultaneos > 1`, cada linha de log leva o prefixo `[#i]` da rodada correspondente.

**Estrutura do arquivo** (só `[leva] pasta` é obrigatório):

```toml
[leva]
nome  = "zeta_N"            # padrão: nome do arquivo
pasta = "Testes/zeta_N"     # destino dos runs: relativo à pasta "Versão Fable", ou caminho absoluto
threads = 8                 # por rodada
simultaneos = 1             # rodadas em paralelo (threads x simultaneos <= núcleos)
pos_processamento = true    # false = sem bin/plot (só λ, α; mais rápido)
tag = ""                    # sufixo nos run_id

[graficos]
por_run = true              # u, g, ψ, D⊥ de cada run em <run>/plots/
png = true                  # converte todos os .agr
ponto_decimal = false       # true = "0.5" nos PNG (LC_ALL=C); false = "0,5"

[base]                      # parâmetros fixos (qualquer chave de &teste ou &plot)
m2 = 1.5
Mtot = 1.9

[grade]                     # produto cartesiano das linhas
zeta = [0.0, 0.5, 1.0]
"NMG=NMZ" = [12, 16, 20]                         # mesmo valor nos dois (malha quadrada)
# "gamma_max,NMG" = [[3.0, 16], [4.5, 24]]       # pares que andam juntos

[[rodada]]                  # rodadas avulsas, somadas à grade (repetir o bloco)
m1 = 0.6043956043956044
m2 = 1.3956043956043955

[[convergencia]]            # observável x parâmetro (repetir o bloco para mais)
param = "NMZ"               # parâmetro de &teste ou "razao" (m2/m1)
y = ["lambda", "xi_medio"]  # lambda, alfa, xi_medio, norma_u
erro_log = true             # também |y - y_ref| em log (ref = último ponto, ou ref = valor)
eixo_x_log = false
filtro = { zeta = 1.0 }     # opcional: só essas rodadas

[[sobreposicao]]            # curvas de vários runs no mesmo gráfico
variar = "zeta"             # uma curva por valor
grandezas = ["u", "psi_z"]  # u, g_gamma, g_z, psi_gamma, psi_z, dperp
legenda = "zeta"            # padrão: o próprio "variar"
xi_medio = true             # retas verticais em ⟨ξ⟩ no gráfico de u
filtro = { NMG = 20 }

[tabela]
colunas = ["lambda", "alfa", "xi_medio", "norma_u", "tempo_s"]   # qualquer coluna do resumo.csv
```

**Como as comparações são agrupadas:** em `[[convergencia]]` e `[[sobreposicao]]`, as rodadas são separadas pelos valores dos *outros* parâmetros que variam na leva, e cada grupo com 2 ou mais rodadas vira um gráfico. No exemplo acima (ζ × N), `param = "NMZ"` dá um gráfico λ × NMZ para cada ζ, e `variar = "zeta"` dá uma sobreposição para cada malha (aqui restrita a `NMG = 20` pelo filtro). Parâmetros ligados (`"NMG=NMZ"`, `"gamma_max,NMG"`) contam como um só. O nome do arquivo leva só o que distingue os grupos, por exemplo `conv_NMZ_zeta-0.5_lambda.agr`.

**O que é gerado:**

```
<pasta>/
    resumo.csv, <run_id>/...                  como no run.py (§6); <run>/plots/ se por_run = true
    leva_<nome>/
        relatorio.md                          tabela (parâmetros que variam + colunas escolhidas), falhas, lista de gráficos
        rodadas.csv                           observáveis só das rodadas desta leva (+ coluna razao)
        csv/<comparação>.csv                  dados de cada gráfico de convergência
        conv_*.agr/.png, sobre_*.agr/.png     comparações
        graficos.log                          saída de todos os comandos de gráfico
        <nome>_<data>.toml                    cópia do arquivo usado (procedência)
```

**Comportamento:**

- **Validação antes de rodar**: seções, chaves, parâmetros (inclusive nos filtros), tipos (`NMG = 8.5` → erro), `y` e grandezas são conferidos antes de qualquer cálculo. Duas rodadas com o mesmo `run_id` e parâmetros diferentes também são recusadas; isso acontece ao variar uma chave de `&plot`, que não entra no `run_id`.
- **Relançável**: cada rodada passa pelo `run.py --skip-existing`. Rodar a mesma leva de novo (ou depois de interrompê-la) só calcula o que falta e regenera gráficos e relatório. Runs já existentes na pasta, de outras levas, são reaproveitados.
- **Rodadas que falharam** ficam listadas no fim e no `relatorio.md`, e o código de saída é 2. Para refazer uma delas, apague a pasta do run.
- **Gráfico sem dados** (ex.: |Δ⟨ξ⟩| = 0 com massas iguais, impossível em escala log) gera só um aviso; o resto continua.
- `pos_processamento = false` é incompatível com `por_run`, `[[sobreposicao]]` e `y = xi_medio/norma_u` (erro na validação).

Testado com uma leva de 7 rodadas (grade ζ × malha 6/8/10 + uma rodada avulsa do káon, 2 simultâneas): ~9 s ao todo, 18 gráficos de comparação, relançamento sem recalcular nada.

### 7.2 Laço simples com `run.py`

Cada chamada gera um `run_id` próprio e uma linha no `resumo.csv`:

```bash
for z in 0 0.5 1; do for n in 16 20 24; do
  python3 tools/run.py --teste varre --set zeta=$z NMG=$n NMZ=$n --no-plot --skip-existing
done; done
```

Com `--skip-existing`, o laço pode ser interrompido e relançado sem refazer nada.

### 7.3 `tools/conv_local.py`: varreduras de convergência prontas (T2)

Varre **um parâmetro por vez** em torno de dois casos-base (malha 16×16, Ng=Nz=32, Nv=Nu=24, γ_max=3), sempre em `Testes/convergencia_inicial/`:

- **E**, massas iguais: `m1=m2=1`, `Mtot=1.99`, `mu=0.5`, ζ=1
- **K**, káon: `m2/m1=2.309`, `Mtot=1.3671`, `mu=1.8` (unidades de m̄)

| varredura | valores |
|---|---|
| `NMG`, `NMZ` | 12, 16, 20, 24 |
| `Ng`, `Nz` | 16, 24, 32, 48, 64 |
| `Nv` | 12, 16, 24, 32, 48 |
| `Nu` (em ζ=0.5) | 12, 24, 36 |
| `gamma_max` | 2, 3, 4.5, 6, 9 |
| `gmax_dens` (NMG ∝ γ_max) | (3, 16), (4.5, 24), (6, 32) |
| `zeta` | 0, 0.5, 1, 2, 3 |
| `z_borda` | 0.9 … 0.9999999 |
| `e` (regulador) | 1e-12 … 1e-2 |
| `fis_razao`, `fis_Mtot`, `fis_mu` | estabilidade física (razão de massas, ligação, massa de troca) |

```bash
python3 tools/conv_local.py --listar               # lista as rodadas e marca as já feitas
python3 tools/conv_local.py --rodar --threads 8    # executa só o que falta
python3 tools/conv_local.py --analisar --png       # tabelas + gráficos
```

`--analisar` escreve em `Testes/convergencia_inicial/varreduras/`:

- `<varredura>.csv`: o subconjunto do `resumo.csv` daquela varredura, ordenado pelo parâmetro;
- `<varredura>_lambda.agr`, `<varredura>_xi_medio.agr`: observável × parâmetro (eixo x em log para `e`);
- `<varredura>_lambda_err.agr`: |λ_N − λ_ref| em escala log (só nas varreduras de malha);
- `tabelas.md`: uma tabela Markdown por varredura (λ, variação relativa, ⟨ξ⟩, norma, assimetria, ψ nas bordas, nº de complexos, tempo) e o ajuste de custo `t = a (NMG·NMZ)² Ng Nz Nv`;
- com `--png`, o PNG de cada `.agr`.

A análise dos resultados está em `orquestracao/relatorios/T2_convergencia.md`.

## 8. Gráficos Xmgrace

Os três scripts escrevem `.agr` que reproduzem o estilo de `tools/template_plot_long.agr`: fonte Times, ξ/ζ/λ/ψ/γ em Symbol, legenda sem caixa em (0.8, 0.8), ticks para dentro com 4 menores. A série 0 é azul tracejada, a 1 vermelho-tijolo contínua, a 2 cinza e as seguintes usam as cores padrão do Grace. O mundo (limites) e os ticks são calculados automaticamente com passos 1/2/5×10^k; o eixo y começa em 0 quando faz sentido, e os eixos log usam décadas inteiras e rótulos `10^k`. A precisão dos rótulos é ajustada para que ticks vizinhos não fiquem iguais, e rótulos longos do eixo y são rotacionados.

Robustez: NaN/Inf (e valores ≤ 0 em eixo log) são descartados com aviso, linhas ilegíveis são puladas e números Fortran (`1.0D-03`, `1.0-100`) são lidos. Avisos vão para o stderr com o prefixo `[agr]`.

### 8.1 `tools/plot_run.py`: todos os gráficos de um run

```bash
python3 tools/plot_run.py Testes/<teste>/<run_id> [--png]
```

| `.agr` gerado | fonte em `derived/` | eixos |
|---|---|---|
| `u_xi` | `dist_long.dat` | u × ξ, com reta vertical em ⟨ξ⟩ e legenda `⟨ξ⟩ = …` |
| `g_gamma` | `autovetoresG.dat` | g × γ |
| `g_z` | `autovetoresZ.dat` | g × z |
| `psi_gamma` | `psi.dat` | ψ × γ |
| `psi_z` | `psi_z.dat` | ψ × z |
| `dperp` | `dist_trans.dat` | D⊥(γ) × γ |

Saídas em `<run_dir>/plots/agr/` e, com `--png`, em `<run_dir>/plots/png/`. Um arquivo ausente só gera aviso e é pulado; o script dá erro apenas se nada puder ser gerado. O `run.py --plots` chama este script automaticamente.

### 8.2 `tools/plot_compare.py overlay`: curvas de vários runs no mesmo gráfico

```bash
python3 tools/plot_compare.py overlay --runs R1 R2 [R3 ...] \
        --quantity {u,g_gamma,g_z,psi_gamma,psi_z,dperp} \
        --legend-key {run_id,zeta,razao,<chave>} --out arq.agr [--vlines] [--png]
```

| opção | efeito |
|---|---|
| `--quantity` | qual grandeza sobrepor (mesmos nomes da tabela 8.1) |
| `--legend-key run_id` | legenda = `run_id` (padrão) |
| `--legend-key zeta` | legenda `ζ = …` |
| `--legend-key razao` | legenda `m2/m1 = …` |
| `--legend-key <chave>` | qualquer chave de `observaveis.txt`/`params.nml` (ex.: `NMZ`, `gamma_max`, `lambda`) |
| `--vlines` | só para `u`: reta vertical em ⟨ξ⟩ de cada run, na cor da série |
| `--png` | grava `arq.png` ao lado do `.agr` |

### 8.3 `tools/plot_compare.py conv`: observável × parâmetro, a partir do CSV

```bash
python3 tools/plot_compare.py conv --csv Testes/<teste>/resumo.csv --param NMZ \
        --y {lambda,alfa,xi_medio,norma_u} --out arq.agr \
        [--filter zeta=1.0 m2=1.0 ...] [--logerr [--ref V]] [--logx] [--png]
```

| opção | efeito |
|---|---|
| `--param` | coluna do CSV usada como eixo x (qualquer uma: `NMG`, `Nv`, `gamma_max`, `e`, `zeta`…) |
| `--y` | observável no eixo y |
| `--filter col=valor` | mantém só as linhas com esse valor (comparação numérica com tol. 1e-9, ou texto); pode repetir |
| `--logerr` | plota \|y − y_ref\| em escala log; a referência é o último ponto (maior x), ou `--ref V` |
| `--logx` | eixo x logarítmico (útil para o regulador `e`) |
| `--png` | grava o PNG ao lado do `.agr` |

Se o mesmo valor de `--param` aparecer mais de uma vez depois dos filtros, há um aviso (refine com `--filter`). Uma coluna inexistente dá erro listando as colunas disponíveis.

### 8.4 `tools/agr.py`: `.agr` de arquivos quaisquer

Serve para qualquer arquivo de colunas, inclusive as saídas do modo direto (§11), que não estão na árvore `Testes/`.

```bash
python3 tools/agr.py saida.agr dados1.dat [dados2.dat ...] \
        --xlabel '\f{Symbol}x' --ylabel '\f{Times-Italic}u' \
        --legend 'N = 16' --legend 'N = 20' [--cols 0,1] [--logx] [--logy] [--title T]
```

- `--cols x,y`: colunas a usar (base 0; padrão `0,1`);
- `--legend`: um por arquivo, na mesma ordem;
- os rótulos aceitam códigos do Grace: `\f{Symbol}g` (γ), `\f{Symbol}x` (ξ), `\f{Symbol}y` (ψ), `\f{Symbol}l` (λ), `\s…\N` (subscrito; ponha `\f{}` antes para sair da fonte Symbol).

Uso como biblioteca (é o que `plot_run.py` e `plot_compare.py` fazem):

```python
import sys; sys.path.insert(0, 'tools')
import agr
x, y = agr.read_columns('derived/dist_long.dat', (0, 1))
agr.write_agr('u.agr', [(x, y, 'N = 16', None)],          # estilo None = padrão por índice
              r'\f{Symbol}x', r'\f{Times-Italic}u',
              vlines=[(0.5, dict(color=agr.COR_CINZA, linestyle=3, linewidth=1.5), '')],
              world=(0, 0, 1, None),                          # None = automático naquele limite
              logy=False, legend_pos=(0.8, 0.8))
```

O estilo de cada série é um `dict` com as chaves opcionais `color`, `linestyle`, `linewidth`, `symbol` (`circle`, `square`, `diamond`, `triangle`, `plus`, `x`, `star`…), `symbol_size` e `comment`. As cores do template estão disponíveis como `agr.COR_AZUL`, `agr.COR_TIJOLO` e `agr.COR_CINZA`.

## 9. Conversão `.agr` → PNG / EPS / SVG

**Automática:** `png = true` no arquivo da leva, `--png` em `plot_run.py`, `plot_compare.py`, `conv_local.py --analisar` e `collect.py`, e `--plots` no `run.py`. Todos chamam:

```bash
gracebat -hardcopy -hdevice PNG -printfile saida.png entrada.agr
```

**Manual, um arquivo ou uma pasta inteira:**

```bash
gracebat -hardcopy -hdevice PNG -printfile u.png u.agr
gracebat -hardcopy -hdevice EPS -printfile u.eps u.agr       # para LaTeX
gracebat -hardcopy -hdevice SVG -printfile u.svg u.agr

for f in Testes/convergencia_inicial/varreduras/*.agr; do
  gracebat -hardcopy -hdevice PNG -printfile "${f%.agr}.png" "$f"
done
```

O `gracebat` aborta (`*** buffer overflow detected ***`) com caminhos de ~256 bytes ou mais. Por isso as ferramentas o chamam de dentro da pasta do PNG, com caminhos relativos. Na conversão manual, faça o mesmo: `cd` para a pasta e use nomes curtos.

Dispositivos testados com o Grace desta máquina: `PNG`, `EPS` e `SVG` funcionam; `PDF` **não** está disponível (para ter PDF, gere o EPS e converta com `epstopdf`). Para editar à mão, abra com `xmgrace arquivo.agr`.

Observação: a vírgula decimal nos rótulos do PNG vem do locale pt_BR do Grace. Para ponto decimal, use `LC_ALL=C gracebat ...` (ou `ponto_decimal = true` na leva).

## 10. Cluster

O plano completo, com custos, está em `cluster/plano.md` (110 jobs, ≈178 CPU·h). O fluxo:

```bash
# 1. gerar a lista de jobs (idx, teste, run_id, parametros, NMA, t_est_h, mem_est_GB, cpu_h)
python3 cluster/make_manifest.py --cpus 16       # -> cluster/manifest.csv + resumo de custo
python3 cluster/make_manifest.py --dryrun        # -> cluster/manifest_dryrun.csv (2 jobs, malha mínima)

# 2. submeter (depois de ajustar as linhas AJUSTAR do .slurm)
sbatch --array=1-110%20 cluster/job_array.slurm
MANIFEST=cluster/manifest_dryrun.csv sbatch --array=1-2 cluster/job_array.slurm   # teste

# 3. agregar e plotar
python3 cluster/collect.py Testes/cluster_convergencia --plots --png
```

- **`make_manifest.py`** monta quatro blocos sobre a malha fina (NMG=NMZ=24, Ng=Nz=Nv=48, Nu=12): A. convergência de malha (E e K); B. estabilidade em γ_max (com densidade de nós constante), `e` e `z_borda`; C. ζ × razão de massas; D. `Mtot` e `mu` em ζ=0 e 1. O tempo de cada job é estimado ajustando `t = a (NMG·NMZ)² Ng Nz Nv + b·(termo 1−ζ)` às rodadas locais de `Testes/convergencia_inicial/resumo.csv`; a memória, como ~6 matrizes NMA×NMA.
- **`job_array.slurm`** lê a linha `$SLURM_ARRAY_TASK_ID` do manifest, **confere** se o `run_id` do manifest é o que o `run.py --print-run-id` gera e roda `run.py --no-csv --skip-existing`. Cada tarefa escreve só no seu `run_id`, então não há arquivo compartilhado entre jobs. As variáveis `RAIZ`, `MANIFEST` e `BIN` podem ser trocadas pelo ambiente. O OpenMP fica com `PROC_BIND=close`, `PLACES=cores` e `STACKSIZE=512M`. Os logs vão para `cluster/logs/`.
- **`collect.py`** varre `Testes/<teste>/*/`, lê `observaveis.txt` e `run_info.txt` e reescreve o `resumo.csv`. Se já houver um CSV diferente, ele é guardado em `resumo.csv.bak-<data>`, e os runs falhos ou incompletos são listados. Com `--plots`, agrupa os runs que diferem em **um** parâmetro de malha (com ≥ 3 valores) e gera, em `comparacoes/`, `conv_<par>_NN.csv`, `conv_<par>_NN_lambda.agr` e `conv_<par>_NN_lambda_err.agr` (+ PNG com `--png`).

**Nada foi submetido ainda**: faltam a partição, o limite de tempo e os módulos do cluster (linhas `AJUSTAR` no `.slurm`).

## 11. Modo direto: vários testes num só `inputs.dat`

Sem o `run.py`, o programa principal lê **vários blocos `&teste`** de um arquivo e roda todos em sequência, como na Versão Original; os alfas ficam juntos em `alfas.dat` e `resumo.dat`.

```fortran
&saida pasta = 'Resultados/2026-11-17_conv_N' /
&plot Nz_plot = 60, Ng_plot = 60, Nv_plot = 60, z_fixo = 0.0, gamma_fixo = 0.0,
      gamma_den = 0.0, xi_den = 0.5, xi_num = 0.2, gamma_num = 2.0 /
&teste nome = 'N16', kernel_type = 'v', NMG = 16, NMZ = 16, Nz = 32, Ng = 32, Nv = 24, Nu = 24,
       Mtot = 1.99, m1 = 1.0, m2 = 1.0, mu = 0.5, zeta = 1.0,
       A_z = 1.0, A_g = 1.0, gamma_max = 3.0, z_borda = 0.9999999, e = 1.0d-6 /
&teste nome = 'N20', ... NMG = 20, NMZ = 20, ... /
```

```bash
source /opt/intel/oneapi/setvars.sh
export OMP_NUM_THREADS=10
"<caminho>/bin/spline_omp" inputs.dat                 # sem argumento, lê ./inputs.dat
"<caminho>/bin/plot" Resultados/2026-11-17_conv_N     # pós-processa TODOS os testes da pasta
```

- **Validação prévia**: todos os blocos são validados antes de qualquer cálculo; um erro de digitação ou dois testes com o mesmo `nome` param a execução logo no início.
- **`alfas.dat`**: uma linha por teste, `i NMG NMZ Nz Ng Nv Nu zeta Alfa dAlfa wr_max nome` (`dAlfa` = diferença para o teste anterior).
- **`resumo.dat`**: uma linha por teste com os parâmetros, `Alfa`, `wr_max`, tempo e norma (a norma é preenchida pelo `bin/plot`).
- **Por teste**: `<nome>_zeta<v>_coeficientes.dat`, `_autovalores.dat`, `_plot_dist_long.dat`, etc.
- **Sem sobrescrita**: se `<pasta>/resumo.dat` já existe, a execução é recusada.
- Os valores omitidos não têm padrão, exceto `Nu` (−1 → `Nv`), `usar_gmin` (ligado) e `DEBUG_TERMOS`/`SALVAR_MATRIZES`/`CHECAR_DENOM` (desligados). O limite é de 500 blocos por arquivo.

Nesse modo não há `observaveis.txt` nem `resumo.csv`. Para ter gráficos, use o `agr.py` diretamente:

```bash
cd Resultados/2026-11-17_conv_N
python3 "<caminho>/tools/agr.py" u_N.agr N16_zeta1_plot_dist_long.dat N20_zeta1_plot_dist_long.dat \
        --xlabel '\f{Symbol}x' --ylabel '\f{Times-Italic}u' --legend 'N = 16' --legend 'N = 20'
gracebat -hardcopy -hdevice PNG -printfile u_N.png u_N.agr
```

## 12. Parâmetros de `&teste` e `&plot`

| parâmetro | significado |
|---|---|
| `nome` | rótulo do teste (único); vazio → gerado automaticamente (o `run.py` usa `run`) |
| `kernel_type` | `'v'` vetorial; qualquer outro valor → escalar |
| `NMG`, `NMZ` | número de splines em γ e z (dimensão do problema: `NMG·NMZ`) |
| `Ng`, `Nz`, `Nv`, `Nu` | pontos de Gauss em γ', z', v e u (`Nu` só no termo `1−ζ`; −1 → `Nv`) |
| `Mtot`, `m1`, `m2`, `mu` | massa do estado ligado, massas dos constituintes e massa da troca |
| `zeta` | calibre (1 = Feynman, 0 = Landau); com ζ = 1 o termo P^(2) não é calculado |
| `A_g`, `A_z` | razão geométrica das malhas de colocação |
| `gamma_max` | corte superior em γ |
| `z_borda` | nós extremos em z: ±`z_borda` |
| `e` | regularização somada à diagonal de X |
| `usar_gmin` | `.TRUE.` (padrão): `γ' = γ̃ + γ_min(z')`; `.FALSE.`: `γ' ≥ 0` (comportamento antigo) |
| `CHECAR_DENOM` | grava `denominadores.dat` (mínimos dos denominadores nos pontos de quadratura) |
| `SALVAR_MATRIZES`, `DEBUG_TERMOS` | gravam `Z` e `X` antes do DGGEV / termos de diagnóstico (só no modo direto) |

`&plot`: `Nz_plot`, `Ng_plot`, `Nv_plot` (resolução das curvas), `z_fixo`, `gamma_fixo`, `gamma_den`, `xi_den`, `xi_num`, `gamma_num` (pontos fixos para as curvas de ψ e D⊥). Pelo `run.py`, todos são alteráveis com `--set`.

## 13. Resultados conhecidos e pontos em aberto

- A convergência local (T2) mostrou que λ **não** converge em `gamma_max`, em `e` e em `z_borda`, ou seja, o erro dominante não é de malha. Detalhes em `orquestracao/relatorios/T2_convergencia.md`.
- A convenção de sinal de γ_min(z) ainda precisa ser confirmada com o orientador (`orquestracao/relatorios/T1_gamma_min.md`).
- Limitações conhecidas dos gráficos (D⊥ muito concentrado perto de γ=0, legenda que pode cobrir a curva, etc.) estão em `orquestracao/relatorios/T4_agr.md`.
