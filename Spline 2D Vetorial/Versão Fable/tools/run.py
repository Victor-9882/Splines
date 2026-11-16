#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Executa UMA rodada (um run_id) dos programas Fortran e organiza as saídas.

Árvore gerada (convenção da T5):

    Testes/<teste>/<run_id>/
        input/params.nml      namelist exata usada (&saida, &plot, &teste)
        input/run_info.txt    data, host, commit, flags, threads, tempo de parede
        raw/                  saídas do programa principal (coeficientes, autovalores, ...)
        derived/              saídas do programa de plot (u(xi), psi, ...) + observaveis.txt
        plots/agr  plots/png  gráficos (tools/plot_run.py)
        logs/                 stdout/stderr dos dois programas
    Testes/<teste>/resumo.csv 1 linha por run_id

Uso:
    python3 tools/run.py --teste baseline
    python3 tools/run.py --teste conv --set NMG=20 NMZ=20 zeta=0.5
    python3 tools/run.py --teste conv --set m2=1.5 --dry-run

Nunca sobrescreve: se o run_id já existir, a rodada é recusada (use --tag para
criar outro run_id, ou --skip-existing para apenas pular).
Só usa a biblioteca padrão do Python.
"""
import argparse
import csv
import datetime
import fcntl
import math
import os
import re
import shutil
import socket
import subprocess
import sys
import time

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ----------------------------------------------------------------------------
# Parâmetros padrão = baseline B0 (massas iguais, zeta = 1, malha pequena)
# ----------------------------------------------------------------------------
TESTE_PADRAO = dict(
    kernel_type='v',
    NMG=16, NMZ=16, Nz=32, Ng=32, Nv=24, Nu=24,
    Mtot=1.99, m1=1.0, m2=1.0, mu=0.5, zeta=1.0,
    A_z=1.0, A_g=1.0, gamma_max=3.0, z_borda=0.9999999, e=1.0e-6,
    usar_gmin=None,          # None -> não escreve (vale o padrão do Fortran)
    CHECAR_DENOM=None,       # idem
)
PLOT_PADRAO = dict(
    Nz_plot=60, Ng_plot=60, Nv_plot=60, z_fixo=0.0, gamma_fixo=0.0,
    gamma_den=0.0, xi_den=0.5, xi_num=0.2, gamma_num=2.0,
)
INTEIROS = {'NMG', 'NMZ', 'Nz', 'Ng', 'Nv', 'Nu', 'Nz_plot', 'Ng_plot', 'Nv_plot'}
LOGICOS = {'usar_gmin', 'CHECAR_DENOM'}
TEXTOS = {'kernel_type'}

RAW = ['coeficientes.dat', 'autovalores.dat', 'alfa.dat', 'parametros.dat',
       'erros.dat', 'termos.dat', 'matrizes.bin', 'denominadores.dat']
DERIVED = {'plot_dist_long.dat': 'dist_long.dat', 'plot_dist_trans.dat': 'dist_trans.dat',
           'plot_psi.dat': 'psi.dat', 'plot_psi_z.dat': 'psi_z.dat',
           'plot_psi2_3d.dat': 'psi2_3d.dat',
           'autovetoresG.dat': 'autovetoresG.dat', 'autovetoresZ.dat': 'autovetoresZ.dat'}
SOLTOS_RAW = ['resumo.dat', 'alfas.dat', 'inputs.dat', 'plot_parametros.dat']

COLUNAS_CSV = ['run_id', 'teste', 'kernel_type', 'NMG', 'NMZ', 'Nz', 'Ng', 'Nv', 'Nu',
               'Mtot', 'm1', 'm2', 'mu', 'zeta', 'A_z', 'A_g', 'gamma_max', 'z_borda', 'e',
               'usar_gmin', 'lambda', 'alfa', 'xi_medio', 'norma_u', 'assim_u',
               'psi_z_borda_menos', 'psi_z_borda_mais', 'i_max', 'n_complexos',
               'lambda_2', 'norma_psi', 'tempo_s', 'tempo_parede_s', 'threads', 'data']


def erro(msg):
    print('ERRO: ' + msg, file=sys.stderr)
    sys.exit(1)


def converte(chave, texto):
    if chave in INTEIROS:
        return int(texto)
    if chave in LOGICOS:
        t = texto.strip().lower().strip('.')
        if t in ('t', 'true', '1', 'sim'):
            return True
        if t in ('f', 'false', '0', 'nao', 'não'):
            return False
        erro('valor lógico inválido para %s: %s' % (chave, texto))
    if chave in TEXTOS:
        return texto
    return float(texto.lower().replace('d', 'e'))


def fmt_nml(chave, v):
    if chave in LOGICOS:
        return '.TRUE.' if v else '.FALSE.'
    if chave in TEXTOS:
        return "'%s'" % v
    if chave in INTEIROS:
        return str(int(v))
    return repr(float(v)).replace('e', 'd') if 'e' in repr(float(v)) else repr(float(v)) + 'd0'


def g(v, n=4):
    """Número curto para compor o run_id."""
    return ('%.*g' % (n, v))


def monta_run_id(t, tag=None):
    partes = ['m1-%s' % g(t['m1'], 5), 'm2-%s' % g(t['m2'], 5), 'M-%s' % g(t['Mtot'], 6),
              'mu-%s' % g(t['mu']), 'zeta-%s' % g(t['zeta']),
              'G%dx%d' % (t['NMG'], t['NMZ']),
              'q%d-%d-%d-%d' % (t['Ng'], t['Nz'], t['Nv'], t['Nu']),
              'gmax-%s' % g(t['gamma_max'])]
    for chave, rot in (('A_g', 'Ag'), ('A_z', 'Az'), ('z_borda', 'zb'), ('e', 'e')):
        if t[chave] != TESTE_PADRAO[chave]:
            partes.append('%s-%s' % (rot, g(t[chave], 8)))
    if t['kernel_type'] != 'v':
        partes.append('kern-%s' % t['kernel_type'])
    if t.get('usar_gmin') is False:
        partes.append('gmin0')
    if tag:
        partes.append(tag)
    return '_'.join(partes)


def escreve_nml(caminho, t, p, comentario):
    with open(caminho, 'w', encoding='utf-8') as f:
        f.write('! %s\n' % comentario)
        f.write("&saida pasta = 'work' /\n")
        f.write('&plot ' + ', '.join('%s = %s' % (k, fmt_nml(k, v)) for k, v in p.items()) + ' /\n')
        itens = ["nome = 'run'"] + ['%s = %s' % (k, fmt_nml(k, v)) for k, v in t.items() if v is not None]
        f.write('&teste ' + ',\n       '.join(', '.join(itens[i:i + 6]) for i in range(0, len(itens), 6)) + ' /\n')


def ambiente_mkl(env):
    """Se o binário precisa do MKL e o ambiente não está carregado, carrega o setvars."""
    if env.get('MKLROOT'):
        return env
    setvars = '/opt/intel/oneapi/setvars.sh'
    if not os.path.isfile(setvars):
        return env
    try:
        saida = subprocess.run(['bash', '-c', 'source %s >/dev/null 2>&1; env -0' % setvars],
                               capture_output=True, check=True).stdout.decode('utf-8', 'replace')
    except (OSError, subprocess.CalledProcessError):
        return env
    novo = dict(env)
    for item in saida.split('\0'):
        if '=' in item:
            k, v = item.split('=', 1)
            novo[k] = v
    return novo


def git_info():
    def roda(*args):
        try:
            return subprocess.run(['git'] + list(args), cwd=RAIZ, capture_output=True,
                                  text=True, check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return ''
    commit = roda('rev-parse', '--short', 'HEAD') or 'sem-git'
    ramo = roda('rev-parse', '--abbrev-ref', 'HEAD')
    sujo = roda('status', '--porcelain', '--', 'src')
    return commit, ramo, bool(sujo)


def le_colunas(caminho):
    xs, ys = [], []
    with open(caminho, encoding='utf-8', errors='replace') as f:
        for linha in f:
            c = linha.split()
            if len(c) < 2 or c[0][0] in '#!':
                continue
            try:
                x, y = float(c[0].replace('D', 'E')), float(c[1].replace('D', 'E'))
            except ValueError:
                continue
            xs.append(x)
            ys.append(y)
    return xs, ys


def trapz(xs, ys):
    return sum(0.5 * (ys[i] + ys[i + 1]) * (xs[i + 1] - xs[i]) for i in range(len(xs) - 1))


def interp(xs, ys, x):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    w = (x - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] * (1 - w) + ys[hi] * w


def observaveis(run_dir):
    """Lê as saídas organizadas e devolve um dict com lambda, alfa, <xi>, etc."""
    obs = {}
    raw, der = os.path.join(run_dir, 'raw'), os.path.join(run_dir, 'derived')

    # resumo.dat: nome,NMG,NMZ,Nz,Ng,Nv,Nu,k,Mtot,m1,m2,mu,zeta,A_z,A_g,Alfa,wr(1),wr_max,tempo,Norma
    arq = os.path.join(raw, 'resumo.dat')
    if os.path.isfile(arq):
        with open(arq, encoding='utf-8', errors='replace') as f:
            for linha in f:
                c = linha.split()
                if len(c) >= 19 and not linha.startswith('#'):
                    obs['alfa'] = float(c[15])
                    obs['lambda'] = float(c[17])
                    obs['tempo_s'] = float(c[18])
                    if len(c) >= 20:
                        try:
                            obs['norma_psi'] = float(c[19])
                        except ValueError:
                            pass

    # autovalores.dat: i, Re, Im  -> complexos, índice usado e 2º maior real positivo
    arq = os.path.join(raw, 'autovalores.dat')
    if os.path.isfile(arq):
        reais, ncomp = [], 0
        with open(arq, encoding='utf-8', errors='replace') as f:
            for linha in f:
                c = linha.split()
                m = re.search(r'Maior autovalor real positivo:\s*(\S+)\s+indice\s+(\d+)', linha)
                if m:
                    obs['i_max'] = int(m.group(2))
                    obs.setdefault('lambda', float(m.group(1)))
                    continue
                if len(c) == 3:
                    try:
                        int(c[0])
                        wr, wi = float(c[1]), float(c[2])
                    except ValueError:
                        continue
                    if wi != 0.0:
                        ncomp += 1
                    elif 0.0 < wr < 1e15:
                        reais.append(wr)
        obs['n_complexos'] = ncomp
        reais.sort(reverse=True)
        if len(reais) > 1:
            obs['lambda_2'] = reais[1]

    # parametros.dat (namelist regravada pelo Fortran): usar_gmin
    arq = os.path.join(raw, 'parametros.dat')
    obs['usar_gmin'] = 'n/a'
    if os.path.isfile(arq):
        with open(arq, encoding='utf-8', errors='replace') as f:
            m = re.search(r'USAR_GMIN\s*=\s*([TF])', f.read(), re.I)
        if m:
            obs['usar_gmin'] = m.group(1).upper()

    # u(xi): norma, <xi>, assimetria u(xi) - u(1-xi)
    arq = os.path.join(der, 'dist_long.dat')
    if os.path.isfile(arq):
        xs, us = le_colunas(arq)
        if len(xs) > 2 and all(math.isfinite(u) for u in us):
            norma = trapz(xs, us)
            obs['norma_u'] = norma
            if norma != 0.0:
                obs['xi_medio'] = trapz(xs, [x * u for x, u in zip(xs, us)]) / norma
            umax = max(abs(u) for u in us)
            if umax > 0.0:
                obs['assim_u'] = max(abs(u - interp(xs, us, 1.0 - x)) for x, u in zip(xs, us)) / umax

    # psi(gamma_fixo, z): valor nas bordas z = -1 e z = +1 relativo ao máximo (pico espúrio)
    arq = os.path.join(der, 'psi_z.dat')
    if os.path.isfile(arq):
        zs, ps = le_colunas(arq)
        if len(zs) > 2 and all(math.isfinite(p) for p in ps):
            pmax = max(abs(p) for p in ps)
            if pmax > 0.0:
                i_menos = min(range(len(zs)), key=lambda i: zs[i])
                i_mais = max(range(len(zs)), key=lambda i: zs[i])
                obs['psi_z_borda_menos'] = ps[i_menos] / pmax
                obs['psi_z_borda_mais'] = ps[i_mais] / pmax
    return obs


def organiza(run_dir):
    """Move as saídas de work/ (geradas agora) para raw/ e derived/ com nomes canônicos."""
    work = os.path.join(run_dir, 'work')
    movidos = []
    if not os.path.isdir(work):
        return movidos
    prefixo = None
    for nome in sorted(os.listdir(work)):
        if nome.endswith('_parametros.dat') and nome != 'plot_parametros.dat':
            prefixo = nome[:-len('parametros.dat')]
    for nome in sorted(os.listdir(work)):
        origem = os.path.join(work, nome)
        destino = None
        if nome in SOLTOS_RAW:
            destino = os.path.join(run_dir, 'raw', nome)
        elif prefixo and nome.startswith(prefixo):
            resto = nome[len(prefixo):]
            if resto in RAW:
                destino = os.path.join(run_dir, 'raw', resto)
            elif resto in DERIVED:
                destino = os.path.join(run_dir, 'derived', DERIVED[resto])
        if destino and not os.path.exists(destino):
            shutil.move(origem, destino)
            movidos.append(os.path.relpath(destino, run_dir))
    if not os.listdir(work):
        os.rmdir(work)
    return movidos


def acrescenta_csv(caminho, linha):
    """Acrescenta uma linha ao resumo.csv com trava (seguro para jobs simultâneos)."""
    with open(caminho, 'a+', newline='', encoding='utf-8') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.seek(0, os.SEEK_END)
        w = csv.DictWriter(f, fieldnames=COLUNAS_CSV, extrasaction='ignore')
        if f.tell() == 0:
            w.writeheader()
        w.writerow(linha)
        f.flush()
        fcntl.flock(f, fcntl.LOCK_UN)


def linha_csv(run_id, teste, t, obs, threads, parede, data):
    linha = {k: '' for k in COLUNAS_CSV}
    linha.update({k: t[k] for k in t if k in COLUNAS_CSV and t[k] is not None})
    for k, v in obs.items():
        if k in COLUNAS_CSV:
            linha[k] = repr(v) if isinstance(v, float) else v
    linha.update(run_id=run_id, teste=teste, threads=threads,
                 tempo_parede_s='%.2f' % parede, data=data)
    return linha


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--teste', required=True, help='nome do teste (subpasta de Testes/)')
    ap.add_argument('--set', nargs='*', default=[], metavar='CHAVE=VALOR',
                    help='parâmetros de &teste ou &plot a alterar em relação ao padrão')
    ap.add_argument('--tag', help='sufixo extra do run_id')
    ap.add_argument('--threads', type=int, default=int(os.environ.get('OMP_NUM_THREADS', '8')))
    ap.add_argument('--bin', default=os.path.join(RAIZ, 'bin'), help='pasta dos executáveis')
    ap.add_argument('--raiz-testes', default=os.path.join(RAIZ, 'Testes'))
    ap.add_argument('--no-plot', action='store_true', help='não roda o programa de plot')
    ap.add_argument('--plots', action='store_true', help='gera .agr/.png com tools/plot_run.py')
    ap.add_argument('--no-csv', action='store_true', help='não acrescenta ao resumo.csv (use collect.py)')
    ap.add_argument('--skip-existing', action='store_true', help='se o run_id já existir, pula sem erro')
    ap.add_argument('--dry-run', action='store_true', help='só mostra run_id e parâmetros; não cria nada')
    ap.add_argument('--print-run-id', action='store_true', help='imprime só o run_id e sai')
    a = ap.parse_args()

    t, p = dict(TESTE_PADRAO), dict(PLOT_PADRAO)
    for item in a.set:
        if '=' not in item:
            erro('use CHAVE=VALOR em --set (recebido: %s)' % item)
        k, v = item.split('=', 1)
        alvo = t if k in t else p if k in p else None
        if alvo is None:
            erro('parâmetro desconhecido: %s (válidos: %s)' % (k, ', '.join(list(t) + list(p))))
        alvo[k] = converte(k, v)
    if t['Nu'] < 0:
        t['Nu'] = t['Nv']

    run_id = monta_run_id(t, a.tag)
    run_dir = os.path.join(a.raiz_testes, a.teste, run_id)
    if a.print_run_id:
        print(run_id)
        return 0
    if a.dry_run:
        print('run_id  :', run_id)
        print('run_dir :', run_dir)
        print('&teste  :', ', '.join('%s=%s' % kv for kv in t.items() if kv[1] is not None))
        print('&plot   :', ', '.join('%s=%s' % kv for kv in p.items()))
        print('(dry-run: nada foi criado)')
        return 0

    if os.path.exists(run_dir):
        if a.skip_existing:
            print('já existe, pulando: ' + run_dir)
            return 0
        erro('run_id já existe (nada foi sobrescrito): %s\n      use --tag para outro run_id ou --skip-existing' % run_dir)
    spline, plot = os.path.join(a.bin, 'spline_omp'), os.path.join(a.bin, 'plot')
    if not os.path.isfile(spline):
        erro('executável não encontrado: %s (rode make)' % spline)

    for sub in ('input', 'raw', 'derived', 'plots/agr', 'plots/png', 'logs'):
        os.makedirs(os.path.join(run_dir, sub))
    agora = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    escreve_nml(os.path.join(run_dir, 'input', 'params.nml'), t, p,
                'gerado por tools/run.py em %s; teste=%s; run_id=%s' % (agora, a.teste, run_id))

    env = ambiente_mkl(dict(os.environ))
    env['OMP_NUM_THREADS'] = env['MKL_NUM_THREADS'] = str(a.threads)
    commit, ramo, sujo = git_info()
    flags = ''
    if os.path.isfile(os.path.join(a.bin, 'flags.txt')):
        with open(os.path.join(a.bin, 'flags.txt'), encoding='utf-8') as f:
            flags = f.read().strip()

    t0 = time.time()
    with open(os.path.join(run_dir, 'logs', 'spline.log'), 'w') as log:
        rc = subprocess.run([spline, 'input/params.nml'], cwd=run_dir, env=env,
                            stdout=log, stderr=subprocess.STDOUT).returncode
    rc_plot = None
    if rc == 0 and not a.no_plot and os.path.isfile(plot):
        with open(os.path.join(run_dir, 'logs', 'plot.log'), 'w') as log:
            rc_plot = subprocess.run([plot, 'work'], cwd=run_dir, env=env,
                                     stdout=log, stderr=subprocess.STDOUT).returncode
    parede = time.time() - t0

    movidos = organiza(run_dir)
    obs = observaveis(run_dir)
    with open(os.path.join(run_dir, 'derived', 'observaveis.txt'), 'w', encoding='utf-8') as f:
        f.write('# observáveis do run_id %s (tools/run.py)\n' % run_id)
        for k, v in t.items():
            if v is not None:
                f.write('%s = %s\n' % (k, v))
        for k in sorted(obs):
            f.write('%s = %s\n' % (k, repr(obs[k]) if isinstance(obs[k], float) else obs[k]))

    with open(os.path.join(run_dir, 'input', 'run_info.txt'), 'w', encoding='utf-8') as f:
        f.write('data             = %s\n' % agora)
        f.write('hostname         = %s\n' % socket.gethostname())
        f.write('teste            = %s\nrun_id           = %s\n' % (a.teste, run_id))
        f.write('git_commit       = %s%s\ngit_branch       = %s\n' % (commit, ' (src modificado)' if sujo else '', ramo))
        f.write('flags_compilacao = %s\n' % flags)
        f.write('threads          = %d\n' % a.threads)
        f.write('tempo_parede_s   = %.2f\n' % parede)
        f.write('rc_spline        = %s\nrc_plot          = %s\n' % (rc, rc_plot))
        f.write('comando          = %s\n' % ' '.join(sys.argv))
        f.write('slurm_job        = %s\n' % os.environ.get('SLURM_JOB_ID', ''))

    if not a.no_csv:
        acrescenta_csv(os.path.join(a.raiz_testes, a.teste, 'resumo.csv'),
                       linha_csv(run_id, a.teste, t, obs, a.threads, parede, agora))

    if a.plots and rc == 0:
        plot_run = os.path.join(RAIZ, 'tools', 'plot_run.py')
        if os.path.isfile(plot_run):
            with open(os.path.join(run_dir, 'logs', 'plot_run.log'), 'w') as log:
                subprocess.run([sys.executable, plot_run, run_dir, '--png'],
                               stdout=log, stderr=subprocess.STDOUT)

    estado = 'ok' if rc == 0 and rc_plot in (0, None) else 'FALHOU (rc_spline=%s, rc_plot=%s)' % (rc, rc_plot)
    print('%s | %s | lambda=%s alfa=%s <xi>=%s | %.1fs | %d arquivos' % (
        estado, run_id, obs.get('lambda', 'nan'), obs.get('alfa', 'nan'),
        obs.get('xi_medio', 'nan'), parede, len(movidos)))
    return 0 if rc == 0 and rc_plot in (0, None) else 2


if __name__ == '__main__':
    sys.exit(main())
