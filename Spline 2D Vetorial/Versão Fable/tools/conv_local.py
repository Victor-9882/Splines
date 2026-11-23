#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Varreduras de convergência locais (T2): um parâmetro de cada vez em torno de um caso-base.

    python3 tools/conv_local.py --listar            # mostra as rodadas, não executa
    python3 tools/conv_local.py --rodar             # executa o que falta (pula run_id existentes)
    python3 tools/conv_local.py --analisar          # tabelas + .agr a partir do resumo.csv

Casos-base (malha 16x16, Ng=Nz=32, Nv=Nu=24, gamma_max=3):
    E  massas iguais   m1=m2=1, Mtot=1.99, mu=0.5, zeta=1
    K  káon            m2/m1=2.309, Mtot=1.3671, mu=1.8, zeta=1 (unidades de mbar)
Só usa a biblioteca padrão.
"""
import argparse
import csv
import math
import os
import subprocess
import sys

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RAIZ, 'tools'))
import run as R  # noqa: E402

TESTE = 'convergencia_inicial'
KAON = dict(Mtot=1.3670631868131868, m1=0.6043956043956044, m2=1.3956043956043955, mu=1.8)
CASOS = {'E': {}, 'K': dict(KAON)}

# parâmetro -> valores (o valor-base está incluído em cada lista)
NUMERICOS = {
    'NMG': [12, 16, 20, 24],
    'NMZ': [12, 16, 20, 24],
    'Ng': [16, 24, 32, 48, 64],
    'Nz': [16, 24, 32, 48, 64],
    'Nv': [12, 16, 24, 32, 48],
    'gamma_max': [2.0, 3.0, 4.5, 6.0, 9.0],
}


def razao(r, Mtot):
    """m1, m2 com mbar = 1 e m2/m1 = r."""
    m1 = 2.0 / (1.0 + r)
    return dict(m1=m1, m2=2.0 - m1, Mtot=Mtot)


def varreduras():
    """Lista de (nome_da_varredura, parametro_x, [dict de overrides])."""
    vs = []
    for caso, base in CASOS.items():
        for par, valores in NUMERICOS.items():
            vs.append(('%s_%s' % (caso, par), par, [dict(base, **{par: v}) for v in valores]))
        # Nu só existe no termo (1-zeta): varre em zeta = 0.5
        vs.append(('%s_Nu_zeta0.5' % caso, 'Nu', [dict(base, zeta=0.5, Nu=v) for v in (12, 24, 36)]))
        # calibre
        vs.append(('%s_zeta' % caso, 'zeta', [dict(base, zeta=z) for z in (0.0, 0.5, 1.0, 2.0, 3.0)]))
        # gamma_max com a densidade de nós em gamma mantida (NMG proporcional a gamma_max)
        vs.append(('%s_gmax_dens' % caso, 'gamma_max',
                   [dict(base, gamma_max=gm, NMG=n) for gm, n in ((3.0, 16), (4.5, 24), (6.0, 32))]))
        # posição dos nós extremos em z e regulador da diagonal de X
        vs.append(('%s_z_borda' % caso, 'z_borda', [dict(base, z_borda=zb) for zb in (0.9, 0.99, 0.999, 0.9999, 0.99999, 0.9999999)]))
        vs.append(('%s_e' % caso, 'e', [dict(base, e=e) for e in (1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)]))
    # físicos (estabilidade): razão de massas com a ligação do káon; ligação e mu com massas iguais
    vs.append(('fis_razao', 'razao', [dict(razao(r, KAON['Mtot']), mu=KAON['mu']) for r in (1.0, 1.5, 2.309090909090909)]))
    vs.append(('fis_Mtot', 'Mtot', [dict(Mtot=M) for M in (1.0, 1.367, 1.6, 1.9, 1.99)]))
    vs.append(('fis_mu', 'mu', [dict(mu=mu) for mu in (0.15, 0.5, 1.0, 1.8)]))
    return vs


def params(over):
    t = dict(R.TESTE_PADRAO)
    t.update(over)
    return t


def run_id(over):
    return R.monta_run_id(params(over))


def le_csv():
    caminho = os.path.join(RAIZ, 'Testes', TESTE, 'resumo.csv')
    if not os.path.isfile(caminho):
        return {}
    with open(caminho, newline='', encoding='utf-8') as f:
        return {l['run_id']: l for l in csv.DictReader(f)}


def num(l, k):
    try:
        return float(l[k])
    except (KeyError, ValueError, TypeError):
        return float('nan')


def rodar(threads, so_listar):
    vistos, fila = set(), []
    for nome, par, lista in varreduras():
        for over in lista:
            rid = run_id(over)
            if rid not in vistos:
                vistos.add(rid)
                fila.append((nome, over, rid))
    feitos = le_csv()
    print('%d rodadas distintas (%d já feitas)' % (len(fila), sum(1 for _, _, r in fila if r in feitos)))
    for i, (nome, over, rid) in enumerate(fila, 1):
        if so_listar:
            print('%3d %-16s %s%s' % (i, nome, rid, '  [feito]' if rid in feitos else ''))
            continue
        if rid in feitos:
            continue
        cmd = [sys.executable, os.path.join(RAIZ, 'tools', 'run.py'), '--teste', TESTE,
               '--threads', str(threads), '--skip-existing', '--set'] + ['%s=%r' % kv for kv in over.items()]
        saida = subprocess.run(cmd, capture_output=True, text=True)
        print('%3d/%d %s' % (i, len(fila), (saida.stdout.strip() or saida.stderr.strip())[:200]), flush=True)


def analisar(png):
    feitos = le_csv()
    pasta = os.path.join(RAIZ, 'Testes', TESTE, 'varreduras')
    os.makedirs(pasta, exist_ok=True)
    linhas_md = []
    for nome, par, lista in varreduras():
        rows = []
        for over in lista:
            l = feitos.get(run_id(over))
            if l is None:
                continue
            l = dict(l)
            l['razao'] = repr(num(l, 'm2') / num(l, 'm1'))
            rows.append(l)
        if len(rows) < 2:
            continue
        rows.sort(key=lambda l: num(l, par))
        caminho = os.path.join(pasta, nome + '.csv')
        with open(caminho, 'w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        linhas_md.append('\n### %s (x = %s)\n' % (nome, par))
        linhas_md.append('| %s | lambda | dlam/lam (vs anterior) | <xi> | d<xi> | norma_u | assim_u | psi_z(-1) | psi_z(+1) | n_compl | i_max | t (s) |' % par)
        linhas_md.append('|---|---|---|---|---|---|---|---|---|---|---|---|')
        ant = None
        for l in rows:
            lam, xi = num(l, 'lambda'), num(l, 'xi_medio')
            dl = abs(lam - ant[0]) / abs(lam) if ant else float('nan')
            dx = abs(xi - ant[1]) if ant else float('nan')
            linhas_md.append('| %s | %.10f | %s | %.6f | %s | %.6f | %.1e | %+.4f | %+.4f | %s | %s | %.1f |' % (
                ('%.8g' % num(l, par)), lam, '—' if ant is None else '%.1e' % dl, xi,
                '—' if ant is None else '%.1e' % dx, num(l, 'norma_u'), num(l, 'assim_u'),
                num(l, 'psi_z_borda_menos'), num(l, 'psi_z_borda_mais'), l.get('n_complexos', ''),
                l.get('i_max', ''), num(l, 'tempo_s')))
            ant = (lam, xi)
        for y in ('lambda', 'xi_medio'):
            cmd = [sys.executable, os.path.join(RAIZ, 'tools', 'plot_compare.py'), 'conv', '--csv', caminho,
                   '--param', par, '--y', y, '--out', os.path.join(pasta, '%s_%s.agr' % (nome, y))]
            if par == 'e':
                cmd.append('--logx')
            subprocess.run(cmd + (['--png'] if png else []), capture_output=True, text=True)
        if not nome.startswith('fis_') and not nome.endswith(('_zeta', '_z_borda', '_e', '_gmax_dens')):
            cmd = [sys.executable, os.path.join(RAIZ, 'tools', 'plot_compare.py'), 'conv', '--csv', caminho,
                   '--param', par, '--y', 'lambda', '--logerr',
                   '--out', os.path.join(pasta, '%s_lambda_err.agr' % nome)]
            subprocess.run(cmd + (['--png'] if png else []), capture_output=True, text=True)

    # tempo x custo: t ~ a * (NMG*NMZ)^2 * Ng*Nz*Nv (zeta = 1)
    pts = []
    for l in feitos.values():
        if num(l, 'zeta') == 1.0 and num(l, 'tempo_s') > 0:
            custo = (num(l, 'NMG') * num(l, 'NMZ')) ** 2 * num(l, 'Ng') * num(l, 'Nz') * num(l, 'Nv')
            pts.append((custo, num(l, 'tempo_s'), num(l, 'threads')))
    if pts:
        a = sum(c * t for c, t, _ in pts) / sum(c * c for c, t, _ in pts)
        desvios = sorted(t / (a * c) for c, t, _ in pts)
        linhas_md.append('\n### Tempo de parede x custo (zeta = 1)\n')
        linhas_md.append('Ajuste por mínimos quadrados de t = a (NMG NMZ)^2 Ng Nz Nv com %d rodadas, %d threads: '
                         'a = %.3e s; razão t/ajuste entre %.2f e %.2f (mediana %.2f).' % (
                             len(pts), int(pts[0][2]), a, desvios[0], desvios[-1], desvios[len(desvios) // 2]))
    saida = os.path.join(pasta, 'tabelas.md')
    with open(saida, 'w', encoding='utf-8') as f:
        f.write('# Varreduras de convergência (geradas por tools/conv_local.py --analisar)\n')
        f.write('\n'.join(linhas_md) + '\n')
    print('tabelas em', saida)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--listar', action='store_true')
    ap.add_argument('--rodar', action='store_true')
    ap.add_argument('--analisar', action='store_true')
    ap.add_argument('--png', action='store_true')
    ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()
    if a.listar or a.rodar:
        rodar(a.threads, so_listar=a.listar and not a.rodar)
    if a.analisar:
        analisar(a.png)
    if not (a.listar or a.rodar or a.analisar):
        ap.print_help()


if __name__ == '__main__':
    main()
