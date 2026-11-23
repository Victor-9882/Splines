# T0 — Relatório de reconhecimento (somente leitura)
Fontes: `Versão Original/` — M=2D_Spline_Vetorial_Diferentes_zeta_omp.f90 (1626 l.), K=kernel.f90 (465 l.), P=Plot_funcoes_zeta.f90 (1008 l.), inputs.dat, RODAR_PARALELO.dat. Mapa condensado: `../MAPA_CODIGO.md`. Nada foi compilado nem executado; `git diff` (leitura) só para saber o que é não commitado.

## 1. Estado do repositório
`git status`: M (zeta_omp.f90), M (inputs.dat), RM (kernel_P2_zeta.f90 → kernel.f90). O diff de M contém as linhas adicionadas `gv(1) = 0.22d0` (M:363) e `CALL legauss(0.22d0,gamma_max,Ng,Y,dY,1.d-15)` (M:374), que substituem `gv(1) = 0.d0` e `legauss(0.d0,...)` do HEAD. Ou seja, γ_min=0.22 é um experimento local não commitado, e P (P:238, 334) continua com 0.

## 2. Mapeamentos e quadraturas (detalhe)
- γ': nenhuma transformação t→γ'. M:374 `CALL legauss(0.22d0,gamma_max,Ng,Y,dY,1.d-15)`; no builder `gp = Y(p)`, `dgp = dY(p)` (M:805-806; comentário M:804 "gamma variando de 0 a infinito" desatualizado). Peso de γ' entra em `peso = splg(k)*splz(l)*dzq*dgp*dv` (M:816) e no LHS `...*dgp` (M:859).
- z': X∈[0,1] (M:373). Ramo superior (θ(z'−z)): `dzq=(1.d0-z)*dX(q)`, `zq=(1.d0-z)*X(q)+z` (M:809-810). Ramo inferior (θ(z−z')): `dzq=(z+1.d0)*dx(q)`, `zq=(z+1.d0)*x(q)-1.d0` (M:833-834). Idêntico em BUILD_ZMATRIX_P2 (M:1127-1133) e TERMOS_* (M:971-972, 1008-1009, 1232, 1235).
- v: `legauss(0.d0,1.d0,Nv,W,dW,1.d-15)` (M:375); u=(1−v)t, t ∈ [0,1] com Nu_ef pontos (M:377; K:415); integração em u é numérica no P2_REG (K:414-424) e analítica no termo delta (K:439 comentário).
- `legauss` (M:1397-1457) é a rotina de Gauss-Legendre (Newton com deflação, tol ZZ=1e-15). `gauleg` (M:1459) e `setgaulag` (M:1502, usa dsyev) existem mas não são chamados.
- Para N=1, legauss usa X=0, DX=0.5 (M:1407-1409) antes do reescalonamento.

## 3. Splines e malhas
- SPLGR1/SPLMD1 (z) e SPLGR2/SPLMD2 (γ) são idênticas, mudando só o COMMON (FActz vs FActg). Constroem spline cúbica natural interpolante: `SPL(J)` = base cardinal (SPL(I)=1 no nó I), logo c(i,j) = g(γ_i,z_j) nos nós. Condições: Q(1,·)=0 (M:1289), C(N,·)=FAK(N,·)=0 (M:1305-1308) → derivada segunda zero nas pontas. Busca do intervalo: `IF (XA .GE. X(I+1) .AND. I .LT. N) GOTO 10`, `IF (I .EQ. 0) I=1` (M:1329-1330): abaixo de X(1) extrapola o cúbico do 1º intervalo; acima de X(N) usa I=N, onde FAK=0 → SPL(N)=1 (constante).
- Z: G1D(IW,-1.d0,N_intervalZ,A_z,1.d0,X) → COLLOC(NCOL=2) → `zv(i+1)=XG(i)`, i=1..2N, depois `zv(1)=-z_borda`, `zv(nmz)=z_borda` (M:346-354). Para NMZ par (=20): N=9, 18 interiores, nós 2..19, extremos nos 1 e 20. Para NMZ ímpar o último interior é sobrescrito pelo extremo (M:354). Limites ±z_borda=0.9999999 < 1.
- γ: G1D(0 → gamma_max) → COLLOC → `gv(i+1)=YG(i)`, `gv(1)=0.22d0`, `gv(nmg)=gamma_max` (M:356-364). **Verificação manual** com A_g=1, gamma_max=3, NMG=20: h=3/9, 1º ponto de colocação = 0.5(h(−0.57735)+h)=0.0704, o 2º=0.2630. Logo gv=(0.22, 0.0704, 0.2630,...): não monotônico, h(2)=gv(2)−gv(1)<0 em SPLGR2 (M:1346). Não confirmado numericamente (não executei), mas aritmética direta.
- Condição g(±1)=0: não há linha de código que a imponha (nem eliminação de funções de base nem linhas extras). Nós em z_borda e kernel com fatores (1±z) apenas; o resultado depende do autovetor.
- Sistema quadrado NMA=NMG·NMZ (M:254); `LWORK = 10*NMA` (M:255); `VL(NMA,2*NMA)` alocado (M:268) mas jobvl='N'.
- SPLMD lê `X(I+1)` com I=N → `zv(NMZ+1)`/`gv(NMG+1)`, alocados NMZ+1/NMG+1 (M:269) e nunca atribuídos; valor não usado graças ao `.AND. I .LT. N` (porém Fortran não garante curto-circuito).

## 4. Montagem
BUILD_ZMATRIX (M:702-894): laços i(γ colocação), j(z colocação) colapsados em OpenMP; dentro k,l (bases) → p,q,r (Gauss γ',z',v). O kernel não depende de k,l, mas é recalculado Nmg·Nmz vezes (M:792-860). Para kernel 'v': `zmatrix = zmatrix - contrib_C0_ku + contrib_f1_ku` (M:823-824), análogo no ramo inferior (M:844-845); escalar soma `contrib_esc*` (M:829, 850). SPLz é recalculado em z' dentro do laço e restaurado em z só no LHS (M:857). Progresso com `!$OMP ATOMIC CAPTURE` (M:866-869).
LHS: M:857-859. XMATRIX + e na diagonal (M:419-421). O RHS inclui `1/(32π²D0)` dentro dos KERNEL_* (K:133,148,172,187,203,218); o P2 vem sem esse fator e é multiplicado em M:411.
BUILD_ZMATRIX_P2 (M:1076-1186): para cada (i,j), q, ramo: P2_REG sobre γ' (p) devolve `valor`; acumula `tmp(:) += SG(:,p)*(valor*dY(p)*dzq)` e contrai com splz (M:1146-1154). Delta: `Z2DEL = tmp(k)*splz(l)` com splz em z_j (M:1159-1175). Observação: os comentários de laço `!z` em i e `!gamma` em j (M:1118-1119) estão trocados (i é γ, j é z).

## 5. Kernel — fórmulas literais
Prefácio do kernel (K:9): `-lambda^2/(32 pi^2) Int dz' Int dgamma' (P1 - (1-zeta) P2) g`.
- D0 K:43-44, Du K:55-60, Dd K:71-76, C0_U K:95-96, C0_D K:111-112, escalares K:133,148 e vetoriais K:172,187,203,218: ver mapa seção 2g. Variantes antigas de C0 estão comentadas (K:91-93, K:107-109): código morto.
- Du2: `DENOM_DU2 = (-4*mu**2*(1 + z)*(-1 + u + v - u*zeta) - Mtot**2*v*(1 + z)*(1 + zl)*(1 + (-1 + v)*z - v*zl) + v*(4*gp*(1 + z) + m1**2*(1 + z) + 2*m1*m2*(1 + z) - 4*g*(-1 + v)*(1 + zl) + m2**2*(1 - 3*z + 4*v*z + 4*zl - 4*v*zl)))/4.d0` (K:320-323).
- Dd2: `(4*mu**2*(-1 + z)*(-1 + u + v - u*zeta) + Mtot**2*v*(-1 + z)*(-1 + zl)*(-1 + (-1 + v)*z - v*zl) + v*(m2**2 - 2*m1*m2*(-1 + z) - 4*gp*z - m2**2*z + 4*(gp + g*(-1 + v)*(-1 + zl)) + m1**2*(1 + 3*z - 4*v*z + 4*(-1 + v)*zl)))/4.d0` (K:335-338).
- NUM_P2: K:358-370 (`d0t,d1,d2t,d3=1-v²zp²,d4=-4(v-1)²,d5=-1`; `NUM_P2 = 4*d0t - D2*(kpt*d1 + 4*d2t + c*M**2*d3 + k2t*d4 - 12*D2*d5)`). γ' só entra via D2.
- Consistência algébrica (conferida à mão): DENOM_D0 ≡ γ + (1−z²)κ² + (Δ−z m̄)² com Δ=(m2−m1)/2, m̄=(m1+m2)/2; ambos em K:43 e K:279. Du2 = Du − μ²(1−ζ)(1+z)u (K:314-315, e termo μ² em K:320 consistente).
- D2>0 é checado (K:421); P2_DELTA exige a<0 e b<0 (K:457).
- O kernel não usa Δ nem m̄ fora de D0_P2 (que serve de diagnóstico: D0_P2 é usada em P2_REG, K:408, ou seja, no cálculo real do P2, enquanto P1 usa DENOM_D0; os dois são iguais).

## 6. Convenção de sinal (crítico, sem correção)
| grandeza | código |
|---|---|
| m̄ | `m = (m1 + m2)/2` M:252; `mb = 0.5d0*(m1 + m2)` K:276 |
| Δ | `Dl = 0.5d0*(m2 - m1)` K:277 (m2−m1, não m1−m2) |
| κ² | `kappa = sqrt(m**2 - 0.25*Mtot**2)` M:253; `k2 = mb**2 - 0.25d0*Mtot**2` K:278 |
| z e Δ juntos | `(Dl - z*mb)**2` K:279 → (Δ − z m̄)² |
| M | `Mtot` (input 2.04) |
| B | NÃO ENCONTRADO como variável; B = m1+m2−Mtot só por inspeção (=0.06) |
| μ | `mu` (0.5) |
| LHS | `g + gp + (1-z**2)*kappa**2 + m**2*z**2` M:859, P:352,386,437 → sem Δ |
Associação de massas por ramo: ramo superior (z'≥z): COEF_C0_U usa m2 (K:95), polo `kut` usa m2 (K:262), Du tem m1²(z+1)+2m1m2(z+1)+m2²(...). Inferior: m1 (K:111, 263).
AMBÍGUO: o LHS (M:859) não contém Δ; para m1≠m2, κ²(1−z²)+m̄²z² ≠ (1−z²)κ²+(Δ−z m̄)² (diferença Δ²−2Δ z m̄). Não sei se o LHS correto deve ter (Δ−z m̄)² (análogo de D0 com γ→γ+γ'); a derivação está no .tex (não lido, fora do escopo).

## 7. DGGEV e autovalor
Chamada M:436-438 com A=ZMATRIX, B=XMATRIX, 'N','V'. λ=α/β (M:455-464, β limiar 1.D-16 → λ=1e16). Critério M:477-484: maior WR com WI==0 exato e 0<WR<1e15; sem candidato → wr(1) (M:485-493). Alfa=1/(16π wr_max) (M:509). INFO≠0: registra e segue (M:445-450); nesse caso nada é escrito em coeficientes. Autovetor gravado: vr(:,i_max) reordenado c(i,j)=vr(i+(j-1)*Nmg) (M:517-521). Não há normalização do autovetor. wr(1) e wr_max vão para resumo (M:534-535). Nenhuma filtragem por |WI|/|WR| tolerante.

## 8. Entrada
inputs.dat (7 linhas): L1 comentário; L2 `&saida pasta=...`; L3-4 `&plot Nz_plot=60,Ng_plot=60,Nv_plot=60,z_fixo=0.0,gamma_fixo=0.0,gamma_den=0.0,xi_den=0.5,xi_num=0.2,gamma_num=2.0`; L5-7 `&teste nome='NMZ=20', kernel_type='v', NMG=20, NMZ=20, Nz=40, Ng=40, Nv=30, Nu=30, Mtot=2.04, m1=1.0, m2=1.1, mu=0.5, zeta=1.0, A_z=1.0, A_g=1.0, gamma_max=3.0, z_borda=0.9999999, e=1.0d-6, DEBUG_TERMOS=.FALSE.`. Parâmetros omitidos num bloco herdam do bloco anterior (M:570-571); os do primeiro bloco são obrigatórios (NAO_DEF=-1e300, validação M:599-614). `pasta` existente com resumo.dat aborta (M:183-188); atual `Resultados/2026-11-10_zeta_Massas_Diferentes TESTE 2` (contém espaços; mkdir entre aspas M:190). Lê o arquivo duas vezes (REWIND M:148, 238). Cópia de inputs.dat é feita para a pasta (M:195).
Validações: NMG,NMZ ∈ [3,499]; κ² real; z_borda ∈ (0,1]; gamma_max>0 (M:617-629). Não valida gamma_max > 0.22 (com gamma_max≤0.22 a malha degenera).
RODAR_PARALELO.dat: notas soltas de shell; nenhum fonte o referencia. Linhas 3-4, 10, 17 compilam as versões `_omp` (sem zeta) ou zeta; só a linha 17-21 é a atual.

## 9. Saídas
Ver tabela do mapa seção 4. Detalhes: autovalores.dat (unid. 10) grava todos os NMA autovalores `(I4,2X,F20.12,2X,F20.12)` (M:469-471), depois autovetor do i_max com 9999ES16.8 (M:513). `erros.dat` fica só com diagnóstico G1D/COLLOC; abrir 14 (M:275) e IW=14 (M:313). `resumo.dat` formato M:534 largura 292 (confere com LARG_RESUMO P:118). Stdout: "Alfa = ..." (M:545), avisos; wr_max não é impresso em stdout (só em arquivos).

## 10. Plot_funcoes_zeta.f90
- Fluxo: lê `&plot` (P:146), para cada linha de resumo relê `_parametros.dat` e `_coeficientes.dat`, refaz as malhas com o mesmo código (P:217-242) e a mesma SPLGR/SPLMD.
- g(γ,z): `soma += c(i,j)*splg(i)*splz(j)`; `g1_00=1.d0` (P:259; o cálculo de referência em P:248-256 é descartado).
- Ψ(γ,ξ) = ∫dγ' g(γ',z)/[γ+γ'+m²z²+(1−z²)κ²]², z=1−2ξ; normalizada por Ψ(γ_den=0, ξ_den=0.5) (P:327-356). Gauss em γ': `CALL legauss(0.d0,gamma_max,Ng,Y,dY,1.d-15)` (P:334, 476), Ng=Ng_plot=60, linear, [0,gamma_max]; não há mapeamento nem γ_min.
- Norma: P:489-531 (Nz_plot pts em ξ, Nv_plot em γ externo). f1 em 3D: P:578-579 (malha 200×200). u(ξ): P:595-637 (400 pts em ξ, γ externo com Nv). D⊥(γ): P:651-694 (800 pts). Não calcula ⟨ξ⟩ (grep "media/xi_med" só encontra comentários): NÃO ENCONTRADO. Não calcula TMD além de f1(γ,ξ)=π|Ψ|².
- Clamps `z=±0.999999` (P:498-499 etc.) em vez de z_borda.
- Parâmetros duplicados/hardcoded: ver mapa seção 5. Os mais perigosos: `gv(1)=0.d0` (P:238) e γ'∈[0,γmax] (P:334) contra M:363, 374. O programa de plot também reinicializa `Nz,Ng,Nv` com os valores `*_plot` (P:200-202), mas os parametros.dat sobrescrevem NAMELIST teste (P:190) antes disso, ok.
- Ψ não inclui Δ no denominador (P:352), mesma forma do LHS (M:859).
- Atenção: unidades 16/17 são reutilizadas por P para arquivos diferentes (leitura de coeficientes em 16 → fechado P:215 → plot_psi.dat em 16, P:323).

## 11. Compilação
Linha 17-21 do RODAR_PARALELO.dat. M precisa estar no mesmo diretório que kernel.f90 (INCLUDE). `INCLUDE` coloca módulos antes de PROGRAM; dentro de K a ordem KERNEL_P1 → KERNEL_P2_ZETA (independentes). Os `.mod` presentes (kernel_p1.mod, kernel_p2_zeta.mod) são artefatos. MKL só aparece comentado.

## 12. Riscos relevantes para γ_min(z)
1. Y/dY são globais (COMMON /ALPHAINT/ no PROGRAM, M:57) e iguais para todos (i,j,q): γ_min(z') exigiria nós/pesos por (j,q,ramo) (ou mapeamento γ'=γ_min+(γmax−γ_min)t), dentro de laços já OpenMP; no ramo superior zq depende de q, logo γ_min(z') varia com q.
2. BUILD_ZMATRIX_P2 pré-tabela SG(k,p)=S_k(Y(p)) (M:1104-1107): deixaria de valer. LHS (M:859) também integra em γ' com Y.
3. A base em γ começa em gv(1)=0.22 (M:363); se γ_min(z')<0.22 a spline extrapola. P precisa do mesmo γ_min(z).
4. γ_min=0.22 vigente é experimento sujo (não commitado) e P está inconsistente com ele.
5. Nós não monotônicos para NMG=20,A_g=1,gamma_max=3 (seção 3).
6. LHS sem Δ (AMBÍGUO, seção 6).
7. Código morto: gauleg, setgaulag, C0 comentado (K:91-93,107-109), `det(1000)`, `XG/YG`, `NMAX_GAUSS` só limite.
Dúvidas pendentes: (a) se γ_min=0.22 é valor físico ou teste; (b) se LHS deve conter Δ; (c) ordem e papel de m1 vs m2 por ramo (só reportado).
