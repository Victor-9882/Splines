#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gera cluster/manifest.csv (índice do array -> parâmetros -> run_id) e as estimativas de custo.

    python3 cluster/make_manifest.py                 # escreve cluster/manifest.csv e imprime o resumo
    python3 cluster/make_manifest.py --dryrun        # escreve cluster/manifest_dryrun.csv (2 entradas, malha mínima)
    python3 cluster/make_manifest.py --cpus 16       # threads por job usadas na estimativa de tempo

Tempo: t = a (NMG NMZ)^2 Ng Nz Nv * (8/threads) [+ termo (1-zeta)], com `a` ajustado às rodadas
locais da T2 (Testes/convergencia_inicial/resumo.csv). Memória: ~6 matrizes NMA x NMA em dupla precisão.
Só usa a biblioteca padrão.
"""
import argparse
import csv
import os
import sys

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RAIZ, 'tools'))
import run as R  # noqa: E402

KAON = dict(Mtot=1.3670631868131868, m1=0.6043956043956044, m2=1.3956043956043955, mu=1.8)
IGUAIS = dict(Mtot=1.99, m1=1.0, m2=1.0, mu=0.5)

# Malha de referência das rodadas do cluster (definida a partir da T2; ver cluster/plano.md)
FINA = dict(NMG=24, NMZ=24, Ng=48, Nz=48, Nv=48, Nu=12, gamma_max=3.0)

VARRE_FINA = {
    'NMG': [20, 24, 32, 40],
    'NMZ': [20, 24, 32, 40],
    'Ng': [32, 48, 64, 96],
    'Nz': [32, 48, 64, 96],
    'Nv': [32, 48, 64, 96, 128],
}
# Estabilidade (os três parâmetros que NÃO convergiram na T2)
GMAX_NMG = [(3.0, 24), (4.5, 36), (6.0, 48), (9.0, 72), (12.0, 96)]     # densidade de nós em gamma constante
REGULADOR = [1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2]
Z_BORDA = [0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999, 0.9999999]


def razao(r, Mtot):
    m1 = 2.0 / (1.0 + r)
    return dict(m1=m1, m2=2.0 - m1, Mtot=Mtot)


def entradas():
    """Lista de (teste, overrides)."""
    es = []
    for caso, fis in (('E', IGUAIS), ('K', KAON)):
        base = dict(fis, zeta=1.0, **FINA)
        # A. convergência: um parâmetro de malha de cada vez
        for par, valores in VARRE_FINA.items():
            for v in valores:
                es.append(('cluster_convergencia', dict(base, **{par: v})))
        # B. estabilidade: corte em gamma, regulador e nós extremos em z
        for gm, nmg in GMAX_NMG:
            es.append(('cluster_estabilidade', dict(base, gamma_max=gm, NMG=nmg)))
        for e in REGULADOR:
            es.append(('cluster_estabilidade', dict(base, e=e)))
        for zb in Z_BORDA:
            es.append(('cluster_estabilidade', dict(base, z_borda=zb)))
    # C. física: calibre x razão de massas (ligação do káon)
    for r in (1.0, 1.5, 2.309090909090909):
        for zeta in (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
            es.append(('cluster_fisica', dict(razao(r, KAON['Mtot']), mu=KAON['mu'], zeta=zeta, **FINA)))
    # D. física: energia de ligação e massa de troca (massas iguais), zeta = 1 e 0
    for zeta in (1.0, 0.0):
        for M in (1.0, 1.367, 1.6, 1.8, 1.9, 1.99):
            es.append(('cluster_fisica', dict(IGUAIS, Mtot=M, zeta=zeta, **FINA)))
        for mu in (0.15, 0.5, 1.0, 1.8):
            es.append(('cluster_fisica', dict(IGUAIS, mu=mu, zeta=zeta, **FINA)))
    return es


def ajusta_tempo():
    """(a, b): t[8 threads] = a (NMG NMZ)^2 Ng Nz Nv + b NMG NMZ Ng Nz Nv Nu [zeta != 1]."""
    caminho = os.path.join(RAIZ, 'Testes', 'convergencia_inicial', 'resumo.csv')
    a, b, n = 1.0e-8, 2.0e-8, 0
    if os.path.isfile(caminho):
        with open(caminho, newline='', encoding='utf-8') as f:
            linhas = list(csv.DictReader(f))

        def c1(l):
            return (float(l['NMG']) * float(l['NMZ'])) ** 2 * float(l['Ng']) * float(l['Nz']) * float(l['Nv'])

        def c2(l):
            return float(l['NMG']) * float(l['NMZ']) * float(l['Ng']) * float(l['Nz']) * float(l['Nv']) * float(l['Nu'])
        z1 = [l for l in linhas if float(l['zeta']) == 1.0 and l['tempo_s']]
        if z1:
            a = sum(c1(l) * float(l['tempo_s']) for l in z1) / sum(c1(l) ** 2 for l in z1)
            n = len(z1)
        zn = [l for l in linhas if float(l['zeta']) != 1.0 and l['tempo_s']]
        if zn:
            b = max(0.0, sum(c2(l) * (float(l['tempo_s']) - a * c1(l)) for l in zn) / sum(c2(l) ** 2 for l in zn))
    return a, b, n


def custo(t, a, b, cpus):
    nma = t['NMG'] * t['NMZ']
    seg = a * nma ** 2 * t['Ng'] * t['Nz'] * t['Nv']
    if abs(1.0 - t['zeta']) > 1e-14:
        seg += b * nma * t['Ng'] * t['Nz'] * t['Nv'] * t['Nu']
    seg *= 8.0 / cpus                      # ajuste feito com 8 threads; escala ~linear na montagem
    seg += 2.0e-9 * nma ** 3               # DGGEV, O(n^3) (desprezível frente à montagem)
    mem = 6 * nma ** 2 * 8 / 1024 ** 3 + 0.2
    return seg, mem


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dryrun', action='store_true', help='manifest de 2 entradas com malha mínima')
    ap.add_argument('--cpus', type=int, default=16, help='threads por job (estimativa de tempo)')
    a = ap.parse_args()

    ta, tb, n = ajusta_tempo()
    if a.dryrun:
        minima = dict(NMG=8, NMZ=8, Ng=12, Nz=12, Nv=8, Nu=8)
        es = [('cluster_dryrun', dict(IGUAIS, zeta=1.0, **minima)),
              ('cluster_dryrun', dict(KAON, zeta=0.5, **minima))]
        destino = os.path.join(RAIZ, 'cluster', 'manifest_dryrun.csv')
    else:
        es = entradas()
        destino = os.path.join(RAIZ, 'cluster', 'manifest.csv')

    vistos, linhas = set(), []
    for teste, over in es:
        t = dict(R.TESTE_PADRAO)
        t.update(over)
        rid = R.monta_run_id(t)
        if (teste, rid) in vistos:
            continue
        vistos.add((teste, rid))
        seg, mem = custo(t, ta, tb, a.cpus)
        params = ' '.join('%s=%r' % (k, v) for k, v in over.items())
        linhas.append(dict(idx=len(linhas) + 1, teste=teste, run_id=rid, parametros=params,
                           NMA=t['NMG'] * t['NMZ'], t_est_h='%.3f' % (seg / 3600.0),
                           mem_est_GB='%.2f' % mem, cpu_h='%.2f' % (seg / 3600.0 * a.cpus)))
    with open(destino, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['idx', 'teste', 'run_id', 'parametros', 'NMA', 't_est_h', 'mem_est_GB', 'cpu_h'])
        w.writeheader()
        w.writerows(linhas)

    tot = sum(float(l['cpu_h']) for l in linhas)
    print('manifest: %s' % os.path.relpath(destino, RAIZ))
    print('jobs: %d   cpus/job: %d   ajuste: a = %.3e s, b = %.3e s (%d rodadas locais, 8 threads)' % (len(linhas), a.cpus, ta, tb, n))
    for teste in sorted({l['teste'] for l in linhas}):
        sub = [l for l in linhas if l['teste'] == teste]
        ts = sorted(float(l['t_est_h']) for l in sub)
        print('  %-22s %3d jobs   t/job: mín %.2f h, mediana %.2f h, máx %.2f h   CPU·h: %.0f   mem máx: %s GB' % (
            teste, len(sub), ts[0], ts[len(ts) // 2], ts[-1], sum(float(l['cpu_h']) for l in sub),
            max(l['mem_est_GB'] for l in sub)))
    print('total: %.0f CPU·h' % tot)


if __name__ == '__main__':
    main()
