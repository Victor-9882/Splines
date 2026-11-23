# T2 — Convergência local (16/11/2026)

91 rodadas, todas concluídas, em `Testes/convergencia_inicial/` (`resumo.csv`; tabelas completas e 64 `.agr`/`.png` em `varreduras/`; driver `tools/conv_local.py`).
Um parâmetro de cada vez em torno da malha-base NMG=NMZ=16, Ng=Nz=32, Nv=Nu=24, γ_max=3, z_borda=0.9999999, e=1e-6, com γ_min ligado.
Casos: **E** (m1=m2=1, Mtot=1.99, μ=0.5, ζ=1) e **K** (káon: m2/m1=2.309, Mtot=1.3671, μ=1.8, ζ=1).

## 1. Parâmetros de discretização (critério: |Δλ|/λ < 1e-4 e Δ⟨ξ⟩ < 1e-4 no último refinamento)

| parâmetro | valores | E: \|Δλ\|/λ último passo | K: \|Δλ\|/λ último passo | K: Δ⟨ξ⟩ último passo | veredito |
|---|---|---|---|---|---|
| NMG | 12,16,20,24 | 1.2e-7 (16→20: 1.6e-4) | 3.5e-5 | 1.7e-2 | λ ok com NMG ≥ 20 |
| NMZ | 12,16,20,24 | 4.0e-4 | 3.4e-4 | 4.8e-3 | **não convergiu** (deriva constante) |
| Ng | 16…64 | 5.7e-8 (32→48: 4.5e-5) | 6.2e-8 | 1.3e-3 | λ ok com Ng ≥ 48 |
| Nz | 16…64 | 3.4e-6 | 1.2e-4 | 2.0e-4 | ok em E; no limite em K |
| Nv | 12…48 | 1.2e-3 | 2.1e-3 (piorando) | 2.6e-3 | **não convergiu** (lento em E, diverge em K) |
| Nu (ζ=0.5) | 12,24,36 | 1e-13 | 1e-14 | 5e-13 | Nu = 12 basta |

- ∫u dξ = 1 dentro de 7e-5 em todas; simetria u(ξ)=u(1−ξ) com massas iguais ≤ 5e-10 (≤ 4e-8 com ligação forte).
- No káon ⟨ξ⟩ oscila na 2ª casa (0.272–0.297) mesmo quando λ está parado em 1e-7 (varredura em Ng): u(ξ) não está convergida; os picos perto das bordas mudam de peso com a malha.
- O valor de ψ nas bordas z=±1 **cresce** com Nv, Nz e NMZ (E: 0.0011 → 0.0173 de Nv=12 a 48; K em z=+1: 0.0065 → 0.097).

## 2. Parâmetros que NÃO convergem (resultado principal)

| parâmetro | E: λ | K: λ |
|---|---|---|
| γ_max = 2, 3, 4.5, 6, 9 (NMG=16) | 0.07608, 0.07717, 0.07838, 0.07930, 0.08074 | 0.02698, 0.03039, 0.03390, 0.03640, 0.03989 |
| γ_max = 3, 4.5, 6 com NMG = 16, 24, 32 (mesma densidade de nós) | 0.07717, 0.07835, 0.07925 | 0.03039, 0.03390, 0.03641 |
| e = 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2 | 0.07812, 0.07891, 0.07802, 0.07717, 0.07604, 0.07360 | 0.03656, 0.03495, 0.03295, 0.03039, 0.02675, 0.02002 |
| z_borda = 0.9, 0.99, 0.999, 0.9999, 0.99999, 0.9999999 | 0.07680, 0.07663, 0.07603, 0.07082*, 0.07683, 0.07717 | 0.02920, 0.02773, 0.03980, 0.08278, 0.05947, 0.03039 |

(*) único caso em que o maior autovalor real positivo não foi o de índice 1 (índice 3).
- **γ_max:** λ cresce sem saturar, e o efeito não vem da resolução (mesmos valores com a densidade de nós mantida). No káon dλ/d ln γ_max ≈ 0.0086 constante de 2 a 9, isto é, crescimento logarítmico: 7–11% por passo. Em E, 1–2% por passo.
- **Regulador e (somado à diagonal de X):** não há patamar. Káon: ~8% por duas décadas; E: ~1%. Para e ≤ 1e-10 o número de autovalores complexos dobra e ⟨ξ⟩ do káon fica errático (0.62, 0.21).
- **z_borda:** λ varia de forma errática; no káon entre 0.028 e 0.083.
- Interpretação (minha, a confirmar): (i) γ_max age como corte ultravioleta e o kernel vetorial na escada pode ter dependência logarítmica do corte; (ii) X é mal condicionada (inversão do tipo Stieltjes) e λ depende de direções quase nulas de X; (iii) nas linhas de colocação dos nós extremos z=±z_borda os integrandos em v e z' têm picos de largura ~(1−z_borda) que a quadratura de Gauss não resolve, o que explicaria também a não convergência em Nv/Nz/NMZ e os picos espúrios perto de z=±1.

## 3. Parâmetros físicos (estabilidade; malha-base)
- ζ = 0, 0.5, 1, 2, 3: λ = 0.0547, 0.0663, 0.0772, 0.0987, 0.1195 (E) e 0.0178, 0.0246, 0.0304, 0.0403, 0.0486 (K). Contínuo e monótono; ζ=1 reproduz o código anterior (T1). Na escada λ não é invariante de calibre: não é bug.
- m2/m1 = 1, 1.5, 2.309 (Mtot=1.3671, μ=1.8): λ = 0.03102, 0.03085, 0.03039; ⟨ξ⟩ = 0.500, 0.378, 0.284.
- Mtot = 1.0, 1.367, 1.6, 1.9, 1.99: λ = 0.0385, 0.0419, 0.0457, 0.0582, 0.0772. μ = 0.15, 0.5, 1.0, 1.8: λ = 0.1495, 0.0772, 0.0553, 0.0413. Tudo suave.
- Autovalores: 40–55% do espectro é complexo (80–214 de 192–384), espúrio. O físico foi rastreado como o maior real positivo; ficou no índice 1 em 90 de 91 rodadas, e o segundo real positivo é ~0.3–0.4 λ (sem cruzamento).

## 4. Recomendações de malha e tempo
- Malha de trabalho para λ a ~1e-3 na discretização: **NMG=24, NMZ=24, Ng=48, Nz=48, Nv=48, Nu=12** (~7 min com 8 threads). Não vale refinar além disso antes de tratar γ_max, e e z_borda: o erro deles é 10 a 100 vezes maior.
- **Tempo de parede:** t ≈ a·(NMG·NMZ)²·Ng·Nz·Nv com a = 1.19e-8 s (8 threads, 79 rodadas com ζ=1; t/ajuste entre 0.81 e 1.99, máquina compartilhada). O termo (1−ζ) e o DGGEV (0.08 s para n=256) são desprezíveis.
- O custo é dominado pelo laço (k,l) externo em `BUILD_ZMATRIX`, que recalcula o mesmo kernel NMG·NMZ vezes. Tirar o kernel desse laço (como já é feito em `BUILD_ZMATRIX_P2`) reduziria o tempo em cerca de duas ordens de grandeza. Não foi feito (fora do escopo; muda a ordem das somas).
