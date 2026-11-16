# T4 - Ferramentas de gráficos Xmgrace (.agr)

## Arquivos criados
- `tools/agr.py` - gerador de .agr (API `write_agr`, `read_columns`, `nice_world`, `nice_tick`; CLI simples).
- `tools/plot_run.py` - gráficos de um run (`u_xi, g_gamma, g_z, psi_gamma, psi_z, dperp`), `--png`.
- `tools/plot_compare.py` - subcomandos `overlay` e `conv`.
- `Testes/_fixture_T4/convergencia/` - 8 runs (NMZ20/28/32/36 zeta1.0; NMZ28/32 zeta0.5 e zeta0.0) copiados de
  `Resultados/2026-11-10_zeta_Massas_Iguais` + `resumo.csv` sintético (lambda=col.18, alfa=col.16 de `resumo.dat`;
  xi_medio e norma_u por trapézio de `dist_long`) + `observaveis.txt`, `params.nml`, `run_info.txt`.
- `Testes/_fixture_T4/comparacoes/` - saídas de overlay/conv (agr+png); `Testes/_fixture_T4/regressao_estilo/`.
  Os fontes `.f90` não foram tocados.

## Uso
```
python3 tools/agr.py saida.agr d1.dat d2.dat --xlabel '\f{Symbol}x' --ylabel '\f{Times-Italic}u' --legend a --legend b [--logx --logy --cols 0,1 --title T]
python3 tools/plot_run.py <run_dir> [--png]
python3 tools/plot_compare.py overlay --runs R1 R2 ... --quantity u --legend-key zeta --out x.agr [--png] [--vlines]
python3 tools/plot_compare.py conv --csv resumo.csv --param NMZ --y lambda --out x.agr [--filter zeta=1.0] [--logerr] [--ref V] [--png]
```
Notas: `--png` em `overlay/conv` grava o PNG ao lado do .agr (mesmo nome). `plot_compare` importa `plot_run` e `agr`
(mesma pasta). Se `observaveis.txt` não tiver uma chave (ex.: m1, m2, zeta), há fallback leve para `input/params.nml`.

## Validações
1. **Fixture**: ok (8 runs, convenção T5; `resumo.csv` lido por `csv.DictReader`).
2. **plot_run --png / overlay / conv (normal e --logerr)**: ok. Todos os 20+ .agr da fixture passaram por
   `gracebat -hardcopy -hdevice PNG` sem nenhuma mensagem. PNGs inspecionados: `u_xi`, `dperp`, `overlay_u_zeta`,
   `overlay_psiz`, `conv_lambda`, `conv_lambda_err`, `conv_norma_err`. Observado: eixos, ticks (minor 4) e frame
   corretos; ξ, ζ, λ, ψ, γ em Symbol corretos; `⟨ξ⟩ = 0.5000` (símbolos `\#{e1}`/`\#{f1}` do Symbol) legível, dentro
   da área; legenda sem caixa em (0.8,0.8); série 0 azul tracejada, série 1 vermelho-tijolo, série 2 cinza.
   Problemas achados e corrigidos durante a validação:
   - rótulo y longo (`|λN - λref|`) era cortado (layout `perp` do template é horizontal) -> rótulos com >7
     caracteres visíveis passam a `layout para` (rotacionado); o rótulo de erro agora é `|Δy|` (referência na legenda);
   - subscritos `\s..\N` precisam de `\f{}` antes (senão herdam a fonte Symbol);
   - ticklabels com `prec 2` repetiam valores em faixas estreitas (0,077 0,077...) -> `prec` agora é calculado
     (nº de algarismos significativos necessários; mín. 2, igual ao template);
   - eixo log: formato `power` com `prec 0` (10^-5), em vez de `0,0001`/`1e-05` misturados.
3. **Regressão de estilo**: dados s0/s1 extraídos do template e regenerados por `agr.py` com os mesmos labels/legendas.
   `diff` das linhas `@` (e mesmo nº de linhas de dados, 806): world e ticks idênticos (0,0,1,2.5; 0.2; 0.5). Diferenças
   restantes (6 linhas):
   - `@timestamp def` - data de geração (esperado);
   - `s0 symbol color`, `s0 symbol fill color`, `s0 errorbar color` (1 no template, 16 no gerado) - no template a série 0
     tem esses itens em preto (sem efeito visual, símbolo é "none"); o gerador usa a cor da série;
   - `s0 comment` e `s1 comment` - o template guarda o caminho do arquivo de origem; o gerador deixa vazio (chave
     `comment` do estilo permite preencher).
4. **Bordas**: ok, sem traceback. NaN/Inf descartados com aviso (`série 'a': 2 ponto(s) NaN/Inf descartado(s)`);
   linhas ilegíveis puladas com aviso; série vazia ignorada, e se todas vazias -> erro claro, exit 1; valores <=0 em eixo
   log descartados com aviso; `run_dir` sem `derived/` ou inexistente -> mensagem clara, exit 1; arquivos ausentes em
   `plot_run` -> aviso e pula (exit 0 se ao menos um gerado); coluna inexistente no CSV lista as colunas disponíveis.
   Parser aceita `D`, `1.0E-003`, e `1.0-100` (Fortran sem letra do expoente).

## Divergências entre o resumo de estilo e o template
- Largura da série 1: o template tem `linewidth 1.5` (resumo dizia 2.0 para todas). Segui o template: série 0 = 2.0
  (tracejada, azul), série 1 = 1.5. Para séries >= 2 usei 2.0 (1.5 ficava fino demais em cinza).
- O `legend loctype view` e `view` conferem. Rótulo y no template é `layout perp` (horizontal): mantido para rótulos
  curtos; trocado automaticamente por `para` em rótulos longos.
- Ticklabel `prec 2` do template virou calculado (coincide com o template quando faixa é 0..1 ou 0..2.5).

## Limitações conhecidas
- Vírgula decimal nos PNG/tela vem do locale pt_BR do Grace (o PNG herda o locale do ambiente); não alterado.
- `dperp` é muito concentrado perto de gamma=0; o gráfico usa o world automático (0..gamma_max, pico em ~68): para
  ver a cauda use `--logy` via `agr.py` ou edite o .agr.
- World automático: y começa em 0 se ymin < 0.5 ymax e todos os dados >= 0, com 4% de folga no topo (nada garante
  que a legenda não cubra a curva em casos específicos; ajuste via parâmetro `world`/`legend_pos` na API).
- Eixo log: world em décadas inteiras (pode ser uma única década); minor ticks 8.
- `overlay` sem `--legend-key` conhecido na chave cai para `run_id` com aviso. Legenda `razao` com massas iguais dá
  `m2/m1 = 1` em todos (correto para esta fixture, mas indistinguível: use `--legend-key zeta`).
- `conv --filter` compara numericamente (tol. 1e-9) ou como texto; valores repetidos de `--param` geram aviso.
- Reta vertical em `plot_run`: azul tracejada; em `overlay`: pontilhada na cor da série (sem legenda própria).
