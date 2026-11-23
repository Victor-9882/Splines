#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Agrega os run_id de um teste em resumo.csv e gera os gráficos de comparação.

    python3 cluster/collect.py Testes/<teste> [--plots] [--png]

Os jobs do cluster rodam com --no-csv (sem escrita concorrente); este script
reconstrói o resumo.csv a partir de <run_id>/derived/observaveis.txt e
<run_id>/input/run_info.txt. Se já existir um resumo.csv diferente, ele é
preservado como resumo.csv.bak-<data>. Só usa a biblioteca padrão.
"""
import argparse
import csv
import datetime
import os
import subprocess
import sys

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RAIZ, 'tools'))
import run as R  # noqa: E402


def le_chave_valor(caminho):
    d = {}
    if os.path.isfile(caminho):
        with open(caminho, encoding='utf-8', errors='replace') as f:
            for linha in f:
                if '=' in linha and not linha.startswith('#'):
                    k, v = linha.split('=', 1)
                    d[k.strip()] = v.strip()
    return d


def coleta(pasta_teste):
    linhas, incompletos = [], []
    for nome in sorted(os.listdir(pasta_teste)):
        run_dir = os.path.join(pasta_teste, nome)
        if not os.path.isdir(os.path.join(run_dir, 'input')):
            continue
        obs = le_chave_valor(os.path.join(run_dir, 'derived', 'observaveis.txt'))
        info = le_chave_valor(os.path.join(run_dir, 'input', 'run_info.txt'))
        if 'lambda' not in obs or info.get('rc_spline') != '0':
            incompletos.append(nome)
            continue
        linha = {k: obs.get(k, '') for k in R.COLUNAS_CSV}
        linha.update(run_id=nome, teste=os.path.basename(os.path.normpath(pasta_teste)),
                     threads=info.get('threads', ''), tempo_parede_s=info.get('tempo_parede_s', ''),
                     data=info.get('data', ''))
        linhas.append(linha)
    return linhas, incompletos


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('pasta_teste')
    ap.add_argument('--plots', action='store_true', help='gera λ × N por parâmetro com tools/plot_compare.py')
    ap.add_argument('--png', action='store_true')
    a = ap.parse_args()
    if not os.path.isdir(a.pasta_teste):
        print('ERRO: pasta não encontrada: ' + a.pasta_teste, file=sys.stderr)
        return 1
    linhas, incompletos = coleta(a.pasta_teste)
    destino = os.path.join(a.pasta_teste, 'resumo.csv')
    tmp = destino + '.tmp'
    with open(tmp, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=R.COLUNAS_CSV)
        w.writeheader()
        w.writerows(linhas)
    if os.path.isfile(destino):
        with open(destino, encoding='utf-8') as f1, open(tmp, encoding='utf-8') as f2:
            igual = f1.read() == f2.read()
        if not igual:
            bak = destino + '.bak-' + datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
            os.rename(destino, bak)
            print('resumo.csv anterior preservado em ' + bak)
    os.replace(tmp, destino)
    print('%d run_id agregados em %s' % (len(linhas), destino))
    if incompletos:
        print('%d run_id incompletos/falhos (fora do CSV):' % len(incompletos))
        for nome in incompletos:
            print('   ' + nome)

    if a.plots and linhas:
        pasta = os.path.join(a.pasta_teste, 'comparacoes')
        os.makedirs(pasta, exist_ok=True)
        # um gráfico λ × N para cada parâmetro de malha que varia, por grupo de parâmetros restantes iguais
        malha = ['NMG', 'NMZ', 'Ng', 'Nz', 'Nv', 'Nu', 'gamma_max']
        fixos = ['m1', 'm2', 'Mtot', 'mu', 'zeta', 'usar_gmin'] + malha
        feitos = 0
        for par in malha:
            grupos = {}
            for l in linhas:
                chave = tuple((k, l[k]) for k in fixos if k != par)
                grupos.setdefault(chave, []).append(l)
            for i, (chave, grupo) in enumerate(sorted(grupos.items())):
                if len({l[par] for l in grupo}) < 3:
                    continue
                sub = os.path.join(pasta, 'conv_%s_%02d.csv' % (par, i))
                with open(sub, 'w', newline='', encoding='utf-8') as f:
                    w = csv.DictWriter(f, fieldnames=R.COLUNAS_CSV)
                    w.writeheader()
                    w.writerows(grupo)
                for extra, suf in (([], ''), (['--logerr'], '_err')):
                    cmd = [sys.executable, os.path.join(RAIZ, 'tools', 'plot_compare.py'), 'conv', '--csv', sub,
                           '--param', par, '--y', 'lambda', '--out', sub[:-4] + '_lambda%s.agr' % suf] + extra
                    subprocess.run(cmd + (['--png'] if a.png else []), capture_output=True, text=True)
                feitos += 1
        print('%d grupos de convergência plotados em %s' % (feitos, pasta))
    return 0


if __name__ == '__main__':
    sys.exit(main())
