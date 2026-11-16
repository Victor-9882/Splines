# MAPA_CODIGO (Versão Original; só leitura)
Abreviações: **M**=2D_Spline_Vetorial_Diferentes_zeta_omp.f90, **K**=kernel.f90, **P**=Plot_funcoes_zeta.f90, **I**=inputs.dat. Detalhes em `relatorios/T0_mapa.md`.

## 1. Módulos, subrotinas, funções
| nome | local | função |
|---|---|---|
| PROGRAM SPLINE_VETORIAL | M:4-692 | lê NAMELIST, malhas, monta, DGGEV, grava |
| ZERAR_PARAMETROS / VALIDAR_E_NOMEAR / ERRO_TESTE / FORMATAR_ZETA | M:574 / 592 / 664 / 676 | internas: defaults NAO_DEF, validação, nome `_zeta<v>` |
| BUILD_ZMATRIX | M:702 | Z (P1) e X (LHS), OpenMP |
| BUILD_ZMATRIX_P2 | M:1076 | Z2REG, Z2DEL (termo 1-ζ), OpenMP |
| TERMOS_P1 / TERMOS_P2 | M:910 / 1194 | diagnóstico em termos.dat (DEBUG_TERMOS) |
| SPLGR1,SPLMD1 / SPLGR2,SPLMD2 | M:1280,1321 / 1339,1380 | spline cúbica natural cardinal em z (COMMON /FActz/) / em γ (/FActg/) |
| legauss | M:1397 | Gauss-Legendre em [XS,XL] (Newton, tol ZZ) |
| G1D / COLLOC | M:1548 / 1576 | malha geométrica razão A / pontos de colocação Gauss (NCOL=2 ou 3) |
| gauleg, setgaulag | M:1459, 1502 | CÓDIGO MORTO (nunca chamados; grep) |
| MODULE KERNEL_P1 | K:28-221 | DENOM_D0 K:39, DENOM_DU K:51, DENOM_DD K:67, COEF_C0_U K:87, COEF_C0_D K:103, KERNEL_{UPPER,LOWER}_ESC K:125/140, _C0 K:163/194, _F1 K:179/210 |
| MODULE KERNEL_P2_ZETA | K:246-465 | POLOS_KMENOS 257, D0_P2 271, KD_P2 288, LD_P2 300, DENOM_DU2 316, DENOM_DD2 331, NUM_P2 352, P2_REG 386, P2_DELTA 440 |
| P (programa sem PROGRAM) | P:1-723 | pós-processamento; cópias de SPLGR*/SPLMD*/legauss/G1D/COLLOC em P:729-997 |

Fluxo: `INCLUDE kernel.f90` (M:2) → ler `&saida` (M:140) → 1ª passada valida todos `&teste` (M:151-174) → cria pasta, resumo.dat/alfas.dat (M:183-222) → por teste (M:241): derivados (M:252-256), malhas G1D+COLLOC (M:346-364), SPLGR1/2 (M:367-368), legauss (M:373-377), BUILD_ZMATRIX (M:380), se P2_ATIVO BUILD_ZMATRIX_P2 e `ZMATRIX = ZMATRIX - (1.d0 - zeta)/(32.d0*PI**2) * (Z2REG + Z2DEL)` (M:411), `xmatrix(i,i)+e` (M:420), DGGEV (M:436), escolha de λ (M:477-509), c(i,j)=vr (M:517-521), grava. Depois `./plot` relê tudo de `<pasta>`.

## 2. Localizações
**a. Mapeamento γ'**: NÃO HÁ mapeamento γ'=f(t). Gauss linear direto: `CALL legauss(0.22d0,gamma_max,Ng,Y,dY,1.d-15)` (M:374); `gp = Y(p)`, `dgp = dY(p)` (M:805-806). Domínio γ' ∈ [0.22, gamma_max], gamma_max=3.0 (I:7) lido; **limite inferior 0.22 HARDCODED** (alteração não commitada; HEAD tinha 0.d0). γ_max é corte finito (não ∞). Não existe γ_min(z).
**b. Quadraturas** (todas Gauss-Legendre `legauss`, tol 1e-15, M:373-377): z': `legauss(0,1,Nz,X,dX)`, mapa ramo sup `dzq=(1.d0-z)*dX(q)`, `zq=(1.d0-z)*X(q)+z` (M:809-810); inf `dzq=(z+1.d0)*dx(q)`, `zq=(z+1.d0)*x(q)-1.d0` (M:833-834); mesmos em M:1128-1132. v: `legauss(0,1,Nv,W,dW)` (M:375), v=W(r) (M:802). u: t em `legauss(0,1,Nu_ef,T,dT)` (M:377), `u = (1.d0 - v)*T(s)` (K:415), peso `dW(r)*v**2*(1.d0 - v)*soma_u` (K:425). Pontos atuais Nz=40, Ng=40, Nv=30, Nu=30. P2_DELTA: só em v (K:453-459). Nu<0 → Nu=Nv (M:259).
**c. Splines**: cúbica NATURAL CARDINAL (2ª derivada nula nas pontas: `Q(1,I)=0`, `C(N,J)=0`, M:1289,1305); SPL(J)=valor no nó J; g(γ,z)=Σ c(i,j) S_i(γ) S_j(z), c = valores nodais. Nós = pontos de colocação. z: `zv(1)=-z_borda`, `zv(nmz)=z_borda` (M:353-354), interiores = COLLOC 2 pts por intervalo de G1D(-1→1, A_z) (M:346-351), z_borda=0.9999999. γ: `gv(1) = 0.22d0`, `gv(nmg) = gamma_max` (M:363-364), interiores COLLOC de G1D(0→gamma_max, A_g) (M:356-361). N_intervalZ=(NMZ-1)/2 (M:314). Nº funções = NMG*NMZ=400 (NMG=NMZ=20). **g(γ,z=±1)=0 NÃO é imposto** (nenhum vínculo; nós extremos em ±z_borda livres). Fora de [zv(1),zv(NMZ)] SPLMD1 extrapola constante (FAK(N,.)=0, M:1305-1308).
**d. Colocação**: pontos = nós; (γ_i,z_j)=(gv(i),zv(j)), `g=gv(i)`, `z=zv(j)` (M:789-790); 2 pts de Gauss (±0.577350269189626) por intervalo (M:1591-1594); extremos incluídos à mão. `index1=(j-1)*Nmg+i` (M:791).
**e. Montagem**: M:778-887 `!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(...) COLLAPSE(2) SCHEDULE(DYNAMIC)` sobre (i,j); dentro: k,l (splines) → p(γ'),q(z'),r(v). RHS: `zmatrix(index1,index2) - contrib_C0_ku + contrib_f1_ku` (M:823-824) e idem `kd` (M:844-845), escalar `+ contrib_escu/escd`. peso=`splg(k)*splz(l)*dzq*dgp*dv` (M:816). LHS (dentro do laço p): `xmatrix(index1,index2) + 1.0d0 / ((g +gp + ((1-z**2)*kappa**2) + m**2*z**2)**2)*splg(k)*splz(l)*dgp` (M:858-859). P2: M:1114-1179 PARALLEL DO COLLAPSE(2) REDUCTION(.AND.:ok), CRITICAL M:1142,1165; tabela SG(k,p)=S_k(Y(p)) fixa (M:1104-1107).
**f. DGGEV**: `CALL DGGEV('N','V',NMA,ZMATRIX,...,XMATRIX,...,ALPHAR,ALPHAI,BETA,VL,...,VR,...)` (M:436-438): Z v = λ X v. λ=α/β se |β|>1e-16 senão 1e16 (M:455-464). Critério: `IF (WI(I) .EQ. 0.d0 .AND. WR(I) .GT. wr_max .AND. WR(I) .LT. 1.D15)` (M:480) → maior real positivo i_max; se nenhum, usa wr(1) (M:485-493). `Alfa = 1.0d0 / (wr_max * 16.0d0 * PI)` (M:509). Partes imaginárias: só exige WI==0 exato.
**g. Denominadores/coeficientes** (K; PI=π; γ'=gp, γ=g, z'=zq/zp/zl):
- D0 `DENOM_D0 = 0.25d0*(4.d0*g + Mtot**2*(z**2 - 1.d0) + 2.d0*m1**2*(z + 1.d0) - 2.d0*m2**2*(z - 1.d0))` K:43-44; forma derivada `D0_P2 = g + (1.d0 - z**2)*k2 + (Dl - z*mb)**2` K:279.
- Du K:55-60: `0.25d0*( Mtot**2*(-v)*(z+1)*(zq+1)*((v-1)*z - v*zq + 1) + v*( m1**2*(z+1) + 2*m1*m2*(z+1) + m2**2*(4*v*z - 4*v*zq - 3*z + 4*zq + 1) - 4*g*(v-1)*(zq+1) + 4*gp*(z+1) ) - 4*mu**2*(v-1)*(z+1) )`. γ' aparece só em `4.0d0 * gp * (z + 1.0d0)` (K:59).
- Dd K:71-76: `0.25d0*( Mtot**2*v*(z-1)*(zq-1)*((v-1)*z - v*zq - 1) + v*( m1**2*(-4*v*z + 4*(v-1)*zq + 3*z + 1) - 2*m1*m2*(z-1) - m2**2*z + m2**2 + 4*(gp + g*(v-1)*(zq-1)) - 4*gp*z ) + 4*mu**2*(v-1)*(z-1) )`. γ' em `4.0d0*(gp + ...) - 4.0d0*gp*z` (K:75).
- C0_U K:95-96 `(-2 + v)*(1 + z)*(m2**2*((-2 + v)*z - v*zq) - g*(-2 + v + v*zq)) - (Mtot**2*(1 + z)**2*(4 + ((-2 + v)*z - v*zq)*(-2 + v + v*zq)))/4.`; C0_D K:111-112 `(-2 + v)*(-1 + z)*(g*(-2 + v - v*zq) + m1**2*((-2 + v)*z - v*zq)) - (Mtot**2*(-1 + z)**2*(4 + (2 + v*(-1 + zq))*((-2 + v)*z - v*zq)))/4.` (sem γ').
- Kernels: `1.d0 / (32*PI**2*D0) * (v**2 / (Du**2)) * Co_ku` K:172 (C0 sup), `2.d0 / (32*PI**2*D0) * (v**2 / Du) * (z + 1.d0)` K:187 (f1 sup); inf K:203 (`Dd**2`,`Co_kd`), K:218 (`(1.d0 - z)`); escalar K:133/148 `(1.d0 + z)**2` / `(1.d0 - z)**2`. Builder soma −C0+f1 (M:824,845).
- P2: `kut = (-4*g - 4*m2**2 + Mtot**2*(1 + z))/(2.*Mtot)` K:262; `kdt = (4*g + 4*m1**2 + Mtot**2*(-1 + z))/(2.*Mtot)` K:263. `KD_P2 = 0.5d0*Mtot*v*(1.d0 - v)*(zp - z)` K:292. `LD_P2 = -v*(1-v)*(g + 0.25*Mtot**2*z*zp) - v*(k2 + 0.25*Mtot**2*v*zp**2 + gp) - (1-v)*mu**2 + mu**2*(1-zeta)*u` K:306-308 (γ' em K:307). DENOM_DU2 K:320-323 (γ': `4*gp*(1 + z)`), DENOM_DD2 K:335-338 (γ': `- 4*gp*z ... + 4*(gp + g*(-1 + v)*(-1 + zl))`). NUM_P2 K:358-370 (sem γ' explícito; usa kmt). P2_REG `valor = -0.5d0*soma/D0` K:428; integrando `dT(s)*c*N/D2**3` K:423. P2_DELTA: `a = -v*(1.d0 - v)*g - v*(k2 + gp + 0.25d0*z**2*Mtot**2) - (1.d0 - v)*mu**2` K:455; `b = a + (1.d0 - zeta)*(1.d0 - v)*mu**2` K:456; `valor = 2.d0*soma`, soma `dW(r)*v**3/(a*b)` K:458,461.
- LHS: `(g +gp + ((1-z**2)*kappa**2) + m**2*z**2)**2` M:859 (sem Δ!).
**h. Convenção Δ/m̄/κ²**: M:252 `m = (m1 + m2)/2`; M:253 `kappa = sqrt(m**2 - 0.25*Mtot**2)` (κ²=m̄²−M²/4); K:276 `mb = 0.5d0*(m1 + m2)`; K:277 `Dl = 0.5d0*(m2 - m1)` (**Δ=(m2−m1)/2**); K:278 `k2 = mb**2 - 0.25d0*Mtot**2`; K:279 `(Dl - z*mb)**2` → **(Δ − z m̄)²** (confere algebricamente com DENOM_D0). M=`Mtot` (input 2.04); energia de ligação B NÃO é variável (B=m1+m2−Mtot=0.06 só derivado); μ=`mu` (0.5). Δ/m̄ só aparecem em D0_P2 (K:276-279); resto usa m1,m2 explícitos. Sinal físico de z: ramo sup (z'≥z) usa m2 (K:87,262), inf usa m1 (K:103,263).
**i. ζ**: `P2_ATIVO = (kernel_type .EQ. 'v') .AND. (DABS(1.d0 - zeta) .GE. TOL_ZETA)` M:260 (TOL_ZETA=1.D-14 M:67); fator `-(1-ζ)/(32π²)` M:411; ζ dentro de LD_P2 K:308, DENOM_DU2 K:320 `(-1 + u + v - u*zeta)`, DENOM_DD2 K:335, P2_DELTA K:456. Com ζ=1 (I:6) P2 nunca é calculado. Kernel escalar ignora ζ (M:414).

## 3. Parâmetros (todos de `inputs.dat`, NAMELIST; ler 1 `&saida`, 1 `&plot` (só P), 1+ `&teste`; blocos herdam valores)
Formato I: L1 comentário `!`; L2 `&saida pasta='Resultados/2026-11-10_zeta_Massas_Diferentes TESTE 2' /`; L3-4 `&plot ... /`; L5-7 `&teste ... /`. Arg 1 da linha de comando substitui `inputs.dat` (M:130-131).
| nome | significado | valor | origem |
|---|---|---|---|
| nome | prefixo; recebe `_zeta1.0` | 'NMZ=20' | I:5 |
| kernel_type | 'v' vetorial, outro escalar | 'v' | I:5 |
| NMG, NMZ | nº splines γ, z | 20, 20 | I:5 |
| Nz, Ng, Nv, Nu | pts Gauss z', γ', v, u (Nu=-1→Nv) | 40, 40, 30, 30 | I:5 |
| Mtot, m1, m2, mu | M, massas, μ | 2.04, 1.0, 1.1, 0.5 | I:6 |
| zeta | calibre ζ | 1.0 | I:6 |
| A_z, A_g | razão geométrica das malhas | 1.0, 1.0 | I:7 |
| gamma_max | corte superior γ,γ' | 3.0 | I:7 |
| z_borda | nós extremos ±z_borda | 0.9999999 | I:7 |
| e | regulador diag. de X | 1.0d-6 | I:7 |
| DEBUG_TERMOS / SALVAR_MATRIZES | flags | .FALSE. / (ausente → .FALSE., M:584-585) | I:7 |
| plot: Nz_plot,Ng_plot,Nv_plot | Gauss do plot | 60,60,60 | I:3 |
| plot: z_fixo, gamma_fixo, gamma_den, xi_den, xi_num, gamma_num | cortes/normaliz. | 0.0, 0.0, 0.0, 0.5, 0.2, 2.0 | I:3-4 |
| γ' limite inf / gv(1) | 0.22 | HARDCODED M:363,374 |
| NCOL=2, iw=14, i_dbg..r_dbg=10, NMAX_SPL=499, TOL_ZETA=1e-14 | | HARDCODED M:226,313,229-235,24,67 |
| tol legauss | 1.d-15 | HARDCODED M:373-377 |
`RODAR_PARALELO.dat`: apenas notas de comandos (compilar/rodar, MKL, OMP_NUM_THREADS=10); NINGUÉM o lê (grep em *.f90 sem ocorrência). Linhas 16-21 = "Testes Zeta".

## 4. Saídas (`<pasta>/<nome>_zeta<v>_...`, prefixo M:272)
| arquivo | un | open / write | conteúdo |
|---|---|---|---|
| resumo.dat | 30 | M:197 / M:199-210, 534 | 1 linha/teste: nome,NMG,NMZ,Nz,Ng,Nv,Nu,k,Mtot,m1,m2,mu,zeta,A_z,A_g,Alfa,wr(1),wr_max,tempo,Norma(preenchida por P, P:707) |
| alfas.dat | 31 | M:213 / M:539 | i,NMG,NMZ,Nz,Ng,Nv,Nu,zeta,Alfa,dAlfa,wr_max,nome |
| `_autovalores.dat` | 10 | M:273 / M:470, 506, 510, 513 | cabeçalho, i, Re λ, Im λ (todos NMA), "Maior autovalor real positivo", Alfa, autovetor i_max |
| `_alfa.dat` | 12 | M:274 / M:304-311, 507, 511 | cabeçalho, λ_max, Alfa |
| `_erros.dat` | 14 | M:275 / só G1D/COLLOC via IW=14 (M:1565-1570, 1620-1622) | só diagnóstico de malha; nenhum WRITE(14) explícito |
| `_coeficientes.dat` | 16 | M:276 / M:525 | matriz c(i=γ linhas, j=z colunas), ES16.8 |
| `_termos.dat` | 18 | M:278 / M:319-343 + TERMOS_* | só se DEBUG_TERMOS |
| `_parametros.dat` | 19 | M:282 / M:287 (NML teste) | reutilizável; P o relê |
| `_matrizes.bin` | 17 | M:425 / M:427 | só se SALVAR_MATRIZES (stream) |
| stdout | 6 | M:545 | "Alfa = ..." e avisos; λ=wr_max só nos arquivos; sem fort.NN implícito (legauss escreve unid. 5/2 só se N≤0, M:1402) |
P grava: `plot_parametros.dat`(19, P:155), `_autovetoresG.dat`(11, P:244: γ, g(γ,z_fixo)), `_autovetoresZ.dat`(15, P:245: z, g(γ_fixo,z)), `_plot_psi.dat`(16, P:323), `_plot_psi_z.dat`(17, P:403), `_plot_psi2_3d.dat`(21, P:540: γ,ξ,f1), `_plot_dist_long.dat`(22, P:592: ξ,u), `_plot_dist_trans.dat`(23, P:648: γ,D⊥), reescreve resumo.dat (25, P:717).

## 5. Plot_funcoes_zeta.f90
Lê: `inputs.dat` (`&saida`,`&plot`) ou `<pasta>/inputs.dat` (arg1) (P:128-152); `resumo.dat` (nomes, P:165-181); `<nome>_parametros.dat` (NML teste, P:185-191); `<nome>_coeficientes.dat` (P:211-215). g=`Σ c(i,j)*splg(i)*splz(j)` (P:254, 347). Ψ(γ,ξ)=Σ_p `g_val*dgp/(D**2)` com `D_num = gamma_ext + gammap + (m**2)*(z_ext**2) + (1.d0 - z_ext**2)*(kappa**2)` (P:521-524), ξ=(1−z)/2, `z = 1 - 2*xi` (P:497), Ψ dividido por psi_den=Ψ(γ_den,ξ_den) (P:392). Norma=π∫dξ∫dγ|Ψ|² com `legauss(0,gamma_max,Nv)`, `legauss(0,1,Nz)` (P:489-490,529). `u(ξ)` = `dist_long + PI*(psi_norm**2)*dgamma_ext` / Norma (P:634,637); D⊥ P:691-694. ⟨ξ⟩: NÃO ENCONTRADO (não calcula). TMD/f1: `PI*(psi_norm**2)/Norma` (P:579). Integral em γ': `legauss(0.d0,gamma_max,Ng,Y,dY,...)` (P:334,476), linear, [0,gamma_max]. Relê do parametros.dat: NMG,NMZ,Mtot,m1,m2,A_z,A_g,gamma_max,z_borda (consistente). **Duplicados/hardcoded divergentes**: `gv(1) = 0.d0` (P:238) vs 0.22 (M:363); γ' em [0,γmax] (P:334) vs [0.22,γmax] (M:374); clamp |z|≤0.999999 (P:498) vs z_borda; N_PLOT fixos 1000/1500/200/400/800; g1_00=1 (P:259, normalização desativada).

## 6. Compilação (RODAR_PARALELO.dat:17-21)
`gfortran -O3 -fopenmp -g -fbacktrace 2D_Spline_Vetorial_Diferentes_zeta_omp.f90 -o spline_omp -llapack -lblas`; `export OMP_NUM_THREADS=10`; `gfortran Plot_funcoes_zeta.f90 -o plot`. kernel.f90 entra por `INCLUDE` (M:2) → compilar M no diretório do kernel; módulos KERNEL_P1 antes de KERNEL_P2_ZETA no próprio K; P não usa módulos nem OpenMP. MKL: só comentado (linhas 25-36, para versão antiga `_omp`). Dependência: LAPACK DGGEV (+ BLAS), OpenMP (gfortran).

## 7. Riscos
1. **γ_min=0.22 hardcoded só em M** (M:363,374; não commitado) e `gv(1)=0.d0` em P (P:238): coeficientes reinterpretados em malha diferente no plot.
2. Com A_g=1, NMG=20, γmax=3: `gv(2)`≈0.0704 < `gv(1)`=0.22 → vetor de nós NÃO monotônico (spline com h<0); deve ser checado antes de qualquer γ_min.
3. LHS (M:859) usa m̄²z² sem Δ, enquanto D0 usa (Δ−z m̄)²: AMBÍGUO se intencional (coincide só com m1=m2); P usa o mesmo LHS.
4. γ_min(z) exigirá: Y/dY (M:374) dependentes de (z ou z'), SG(k,p) tabelado independente de z (M:1104-1107) e X/Z no mesmo laço p; spline γ definida em [gv(1),γmax].
5. g(z=±1)=0 não imposto; fatores (1±z) vêm do kernel; SPLMD acessa X(N+1) não inicializado (zv(NMZ+1)) sem efeito.
6. Laço (k,l) externo recalcula kernel idêntico NMG·NMZ vezes (M:792-860): custo desnecessário.
7. ζ: zeta<0 só aviso (M:261); P2 com ok=.FALSE. aborta (M:397-403).

## 8. Revisão do orquestrador (16/11/2026; conferido por amostragem em M:346-376, 470-512, 852-859; K:39-44, 272-279; P:232-240, 334, 352, 476)
- Localizações e expressões literais dos itens 2a, 2c, 2e, 2f, 2g(D0, LHS), 2h conferidas: corretas.
- Risco 3 **resolvido, não é bug**: o LHS `g+gp+(1-z²)κ²+m̄²z²` é exatamente o da TG 3.109 (γ+γ'+κ²+z²M²/4). A diferença para D0 é D0 − [γ+(1−z²)κ²+z²m̄²] = Δ² − 2 z m̄ Δ, que é o próprio γ_min(z) na convenção com sinal. Com γ̃ = γ' − γ_min(z) o LHS vira γ̃ + D0(γ,z).
- **Convenção do código:** Δ = (m2−m1)/2 com sinal (K:277), D0 contém (Δ − z m̄)²; z=+1 ↔ m1, z=−1 ↔ m2 (DENOM_D0 = γ − (1−z²)M²/4 + m1²(1+z)/2 + m2²(1−z)/2). Forma candidata consistente: γ_min(z) = −2 z m̄ Δ + Δ² (sem módulo). A validar na T1 pelo teste m1↔m2.
- Risco 2 confirmado por aritmética: N_intervalG=(NMG−1)/2=9 (M:315), h=1/3, gv(2)=0.0704 < gv(1)=0.22. O estado não commitado (0.22 em M:363 e M:374) dá nós não monotônicos e é inconsistente com P (0.d0). O HEAD usa 0.d0 nos dois lugares.
- ⟨ξ⟩ não é calculado por nenhum programa: obter de `_plot_dist_long.dat` (ξ, u; 401 pontos) por integração externa.
- Custo: rodada de 10/11 (m1=m2=1, Mtot=1.99, ζ=1, NMG=NMZ=20, Nz=Ng=40, Nv=30) levou 85 s; escala ≈ (NMG·NMZ)²·Ng·Nz·Nv.
