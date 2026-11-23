# Orquestração — IC: kernel vetorial, calibre arbitrário ζ, dois bósons de massas distintas

> **Como usar:** coloque este arquivo na raiz da pasta do código (`.../Spline 2D Vetorial/`) e abra o Claude Code com o Fable nessa pasta. Primeira mensagem:
> *"Leia `ORQUESTRACAO_IC_kernel_vetorial.md` e execute só a Fase 0. Pare no checkpoint."*

---

## 0. Papel e regras

Você (Fable) é o **orquestrador**. Você planeja, delega e revisa. Não implementa tudo sozinho e não lê tudo sozinho.

### Regras do usuário (inegociáveis)
1. **Plano antes de executar.** Toda fase começa com um plano curto (≤15 linhas) e **espera aprovação** do Victor.
2. **Nada é apagado, sobrescrito ou renomeado** sem mostrar antes a lista exata (origem → destino) e receber confirmação. Isso vale também para mover resultados antigos.
3. **Não mexer em nada fora da pasta de trabalho.** Não instalar pacotes de sistema sem perguntar.
4. **Não submeter jobs no cluster.** Só preparar scripts e fazer dry-run.

### Economia de tokens (obrigatório para você e para os subagentes)
- **Memória compartilhada em arquivo:** tudo que for aprendido vai para `orquestracao/`. Os subagentes leem `orquestracao/MAPA_CODIGO.md`, e não o código inteiro de novo.
- **Relatórios curtos:** cada subagente devolve ≤30 linhas e grava o relatório completo em `orquestracao/relatorios/T<n>_<nome>.md`. Você lê só o resumo, a menos que algo falhe.
- **Nunca despejar dados no contexto:** para `.dat`, use `head -3`, `wc -l` e resumos via `awk`/Python (mín/máx/norma). Nada de `cat` em arquivos de saída.
- **Compilação:** `make 2>&1 | grep -iE "error|warning" | head -20`. **Execução:** redirecionar para `logs/*.log` e ler só `tail -15`.
- **O PDF do TG (56 p.) não deve ser lido inteiro.** As equações necessárias estão na Seção 2. Se precisar de algo, use `pdftotext -f <p> -l <p>` em páginas específicas (o Cap. 3 está nas pp. ~21–37; a eq. 3.109 está na p. 36).
- **Modelos:** delegue tarefas mecânicas (mapear código, escrever o gerador de `.agr`, reorganizar saídas, scripts SLURM) a subagentes de modelo mais barato. A física do γ_min e a análise de convergência ficam com você.
- **Rodadas locais pequenas:** malhas grossas, cada execução local em menos de ~2 min. Malhas finas vão para o plano do cluster.

### Estado/progresso
Mantenha `orquestracao/STATUS.md` (tabela: tarefa | estado | branch/commit | próximo passo). Atualize ao fim de cada tarefa. Se a pasta não for um repositório git, **pergunte** antes de rodar `git init`. Use uma branch por tarefa e um commit por entrega.

---

## 1. Contexto do código (confirmar na Fase 0)

- Fortran 90. Arquivos conhecidos: `2D_Spline_Vetorial_Diferentes.f90` (monta as matrizes e resolve o autovalor generalizado via `DGGEV`), a variante `_omp.f90` (OpenMP + MKL) e `Plot_funcoes.f90` (reconstrói g, ψ na frente de luz, PDF u(ξ), TMD).
- Compilação MKL usada: `gfortran -O3 -fopenmp -m64 ... -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core` (após `source /opt/intel/oneapi/setvars.sh`).
- g(γ,z) é expandida em splines 2D (γ × z). As integrais em γ', z', v, u são feitas por quadratura de Gauss.
- O calibre ζ entra via o termo (1−ζ)P^(2) (TG eq. 3.109). **ζ = 1 deve reproduzir o código anterior (Feynman).**
- Problema conhecido: na versão káon (massas distintas), ψ reconstruída mostra um pico espúrio em z = −1, atribuído à falta de g(γ, z=±1)=0. **Hipótese a testar na T1:** a ausência de γ_min(z) pode contribuir para esse pico.

---

## 2. Física necessária (resumo, para não abrir o PDF)

Notação do TG (bósons m₁, m₂, massa ligada M, troca vetorial de massa μ):

- m̄ = (m₁+m₂)/2, Δ = (m₂−m₁)/2, κ² = m̄² − M²/4
- Equação (TG 3.109), forma atual com γ' a partir de **0**:
  ∫₀^∞ dγ' g(γ',z) / [γ + γ' + κ² + z²M²/4]² = (−λ²/32π²) ∫₋₁¹ dz' ∫₀^∞ dγ' [P⁽¹⁾ − (1−ζ)P⁽²⁾] g(γ',z')
  (note que κ² + z²M²/4 = (1−z²)κ² + z²m̄²)
- D₀(γ,z) = γ + (1−z²)κ² + (Δ − z m̄)²   (TG 3.87)

**Limitante inferior (referência: artigo bóson–férmion, eqs. 15–16, com Δ = (m_S − m_F)/2):**
  γ_min(z) = −2 z m̄ |Δ| + Δ²
Ele vem da exigência de ausência de cortes (sem propagação livre no estado ligado). Com ele, os dois limites inferiores das integrais em γ' passam a ser γ_min(z) no lado esquerdo e γ_min(z') no lado direito.

**Checagens analíticas que a implementação deve satisfazer:**
- γ_min(z) + z²m̄² = (z m̄ − |Δ|)² ≥ 0, logo o denominador do lado esquerdo fica estritamente positivo para κ² > 0.
- Δ = 0 ⇒ γ_min ≡ 0 ⇒ o resultado deve ser **idêntico** ao código atual (regressão).
- γ_min pode ser **negativo** (z > 0 com Δ ≠ 0). Isso é esperado: o mapeamento da quadratura em γ' não pode assumir γ' ≥ 0.
- ⚠️ **Convenção de sinal:** o artigo usa |Δ| com m_S > m_F. No nosso código, Δ = (m₂−m₁)/2 tem sinal, e a forma consistente com D₀ seria γ_min = −2 z m̄ Δ + Δ² (pois D₀ contém (Δ − z m̄)²). Antes de codificar, verifique qual convenção o código usa em D₀, D_u e D_d e **documente a decisão**. Teste de consistência: trocar m₁ ↔ m₂ deve deixar λ invariante e levar g(γ,z) → g(γ,−z) (se falhar, a convenção está errada). Em caso de dúvida, levante a questão para o Victor/orientador em vez de escolher em silêncio.

---

## 3. Grafo de tarefas

```
Fase 0 (mapa + baseline)
   ├── T1  γ_min(z) ───────────────┐
   └── T5  organização de saídas ──┤
              └── T4  .agr automático ─┤
                                       └── T2  convergência local ── T3  plano do cluster
```
T1 e T5 podem rodar em paralelo, com subagentes diferentes e branches diferentes. T4 depende da convenção de pastas da T5. T2 precisa de T1, T4 e T5 prontas.

---

## Fase 0 — Reconhecimento (subagente barato, somente leitura)

Produzir `orquestracao/MAPA_CODIGO.md` (≤150 linhas) com:
- módulos/subrotinas e o que cada um faz (1 linha cada);
- **onde** estão: o mapeamento da quadratura em γ' (ex.: γ' = f(t), t∈[0,1)); a base de splines em γ e z (nós, domínio, γ_max); os pontos de colocação; a montagem das matrizes; a chamada ao `DGGEV` e o critério de escolha do autovalor físico;
- **todos os parâmetros numéricos e físicos** (nome no código, valor atual e onde são lidos: hardcoded/namelist/stdin);
- **todos os arquivos de saída** (nome, quem escreve, colunas);
- como `Plot_funcoes.f90` lê os coeficientes e reconstrói g, ψ e u(ξ).

Depois, com aprovação: **rodada baseline** com massas iguais e ζ=1, malha pequena. Guardar λ e ⟨ξ⟩ em `orquestracao/baseline.md`. Essa é a referência de regressão.

**Checkpoint:** apresentar o mapa (resumo de 10 linhas) e os planos de T1 e T5.

---

## T1 — Limitante inferior γ_min(z)

**Abordagem recomendada (validar com o mapa):** usar a variável deslocada γ̃ = γ' − γ_min(z') ∈ [0, ∞) e expandir g(γ̃, z') nas splines atuais. Vantagens: a base, os nós e o mapeamento da quadratura ficam iguais, e γ_min só entra onde γ' aparece explicitamente:
- no lado esquerdo: o denominador vira γ̃ + γ + (1−z²)κ² + (z m̄ − Δ)²;
- no lado direito: em todo lugar onde o kernel (D_u, D_d, coeficientes B, P⁽¹⁾, P⁽²⁾, f₀, C₀, l_d...) contém γ', substituir γ' = γ̃ + γ_min(z').
- ⚠️ `Plot_funcoes.f90` **também** precisa da mesma transformação (a reconstrução de ψ, u(ξ) e da TMD integra em γ'). Esquecer isso gera um resultado inconsistente sem nenhum erro aparente.

Passos:
1. Implementar `gamma_min(z, mbar, Delta)` como função única, usada pelos dois programas.
2. Fazer as substituições, listando arquivo:linha de cada uma no relatório.
3. Testes (cada um com critério objetivo):
   - Δ=0: λ igual ao baseline (|Δλ/λ| < 1e-12);
   - massas distintas: todos os denominadores > 0 nos pontos de quadratura (assert temporário em modo debug);
   - troca m₁ ↔ m₂: λ invariante e g(γ,z) ↔ g(γ,−z);
   - ζ=1 com massas distintas: comparar λ antes e depois de γ_min e reportar a diferença (ela é esperada e deve ser explicada);
   - registrar se o pico espúrio em z=−1 de ψ muda.

**Checkpoint antes do merge:** mostrar o diff resumido e a tabela de testes.

---

## T5 — Organização dos arquivos gerados

Hoje cada execução despeja muitos arquivos soltos. Proposta de convenção (ajustar na Fase 0):

```
Testes/<nome_do_teste>/
  <run_id>/                      # ex.: m21-2.30_zeta-1.00_Ng40_Nz30_NG24
    input/params.nml             # cópia exata dos parâmetros usados
    input/run_info.txt           # data, hostname, commit git, flags de compilação, tempo de parede
    raw/                         # .dat escritos pelo Fortran (coeficientes, g, autovalores)
    derived/                     # ψ, u(ξ), TMD, ⟨ξ⟩ (saídas do Plot_funcoes)
    plots/agr/  plots/png/
    logs/
  resumo.csv                     # 1 linha por run_id: parâmetros + λ + ⟨ξ⟩ + tempo
```
Implementação mínima: os programas Fortran recebem `outdir` (namelist ou argumento) e prefixam todos os `open(...)`. Um wrapper `run.sh` (ou `run.py`) cria as pastas, grava `run_info.txt`, roda, cola a linha no `resumo.csv`. **Não mover resultados antigos** sem apresentar a lista e receber confirmação.

---

## T4 — Plots automáticos no Xmgrace

Escrever `tools/agr.py` (Python puro, sem dependências) que gera `.agr` no **estilo do template** `plot_long.agr`. Copie o arquivo para `tools/template_plot_long.agr` e reaproveite o cabeçalho (fontes/cores/página). Não releia o template inteiro: o estilo está resumido abaixo.

Estilo de referência (extraído do template):
- `@version 50125`, página 792×612; `view 0.18, 0.15, 1.15, 0.85`; fonte 0 (Times-Roman).
- Cores extras: 16 `azul` (0,0,205), 17 `vermelho-tijolo` (178,34,34), 18 `cinza` (100,100,100).
- Eixos: barra/ticks/frame com linewidth 1.5; ticks `in`, `place both`, minor ticks 4; ticklabel char size 1.30; label char size 1.75.
- Labels com códigos do Grace: `"\f{Symbol}x"` (ξ), `"\f{Times-Italic}u"`; legenda `"m\s2\N/m\s1\N = 2.3"`; legenda em `loctype view` (0.8, 0.8), sem caixa, char size 1.2.
- Séries: sem símbolo, linewidth 2.0; série 0 tracejada (linestyle 3) em azul; série 1 contínua em vermelho-tijolo; seguintes em cinza e cores padrão.
- `world` e `tick major` calculados dos dados (arredondados).

API sugerida: `write_agr(path, series=[(x, y, legenda, estilo)], xlabel, ylabel, world=None, vlines=[(x, estilo, legenda)])`. Use `vlines` para marcar ⟨ξ⟩ com retas verticais tracejadas.

Plots padrão gerados para cada `run_id` (script `tools/plot_run.py <run_dir>`):
- u(ξ); fatias g(γ, z_fixo) e g(γ_fixo, z); ψ; D⊥(γ) se existir.

Plots de comparação (script `tools/plot_compare.py`):
- vários run_id sobrepostos (ex.: m₂/m₁ = 1 vs 2.3; ζ = 0, 1, 3);
- convergência: λ vs N (por parâmetro), com eixo log opcional para |λ_N − λ_ref|.

Validação: se `gracebat` existir, `gracebat -hardcopy -hdevice PNG -printfile x.png x.agr` sem erros, e abrir o PNG para conferir. Se não existir, validar regenerando o próprio `plot_long.agr` a partir dos dois `.dat` e comparando (diff só nas linhas de dados/comentário).

---

## T2 — Testes de convergência iniciais (local)

**Observáveis:** λ (principal); ⟨ξ⟩; norma ∫u(ξ)dξ (deve ser 1); max|Δu(ξ)| entre refinamentos; com massas iguais, simetria u(ξ) = u(1−ξ).

**Parâmetros a varrer, um de cada vez em torno de um baseline** (nomes reais vêm do mapa):
- N de splines em γ; N de splines em z;
- pontos de Gauss em γ', z', v, u (separadamente, se o código permitir);
- escala/γ_max do mapeamento em γ;
- físicos (não são convergência, mas exigem estabilidade): ζ ∈ {0, 0.5, 1, 2, 3}, m₂/m₁ ∈ {1, 1.5, 2.3}, energia de ligação B, μ.

**Critérios:** convergido quando |λ_N − λ_{N−1}|/λ < 1e-4 e Δ⟨ξ⟩ < 1e-4. Rastrear o autovalor físico entre execuções (o `DGGEV` devolve muitos), por continuidade ou overlap com a solução anterior. **Reportar autovalores espúrios ou complexos.**

**Lembretes físicos:** na aproximação ladder, λ **não** é invariante de calibre. Variar ζ deve dar resultados contínuos e suaves, e ζ=1 deve reproduzir o código anterior. Não interpretar a dependência em ζ como bug.

**Entregas:** `Testes/convergencia_inicial/resumo.csv`; um `.agr` de λ×N por parâmetro (via T4); `orquestracao/relatorios/T2_convergencia.md` com uma tabela de recomendações de malha e o **tempo de parede vs N** (necessário para a T3).

---

## T3 — Plano de convergência em larga escala (cluster)

Produzir, **sem submeter**:
- `cluster/plano.md`: matriz de parâmetros (malhas finas + varredura física), número de jobs, estimativa de tempo e memória por job (extrapolada do ajuste tempo×N da T2: o `DGGEV` escala como O(n³) com n = N_γ·N_z) e custo total em CPU·h;
- `cluster/manifest.csv` (índice do array → parâmetros → run_id);
- `cluster/job_array.slurm`: `--array`, `--cpus-per-task`, `OMP_NUM_THREADS = MKL_NUM_THREADS = $SLURM_CPUS_PER_TASK`, módulos/ambiente MKL, escrita em `Testes/<teste>/<run_id>/` (convenção da T5), sem colisão entre jobs;
- `cluster/collect.py`: agrega `resumo.csv` e chama `plot_compare.py`;
- dry-run local de 2 entradas do manifest com malha mínima.

**Perguntar ao Victor:** qual cluster (partições, limites de tempo, módulos disponíveis) antes de fixar diretivas.

---

## 4. Definição de pronto (para o relatório final ao Victor)
- [ ] T1: γ_min implementado nos dois programas; os 4 testes passaram; convenção de sinal documentada.
- [ ] T5: uma execução gera exatamente a árvore acima; nada antigo apagado.
- [ ] T4: `plot_run.py` e `plot_compare.py` geram `.agr` que abrem no Xmgrace com o estilo do template.
- [ ] T2: tabela de convergência, malha recomendada, tempo×N.
- [ ] T3: plano, manifest, script SLURM e dry-run ok.
- [ ] `STATUS.md` atualizado. Relatório final ≤40 linhas: o que mudou, resultados-chave (λ, ⟨ξ⟩) e dúvidas para o orientador.
