#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""leva.py -- roda uma leva de testes descrita num arquivo TOML e gera tudo.

    python3 tools/leva.py levas/minha_leva.toml              # roda o que falta e gera tudo
    python3 tools/leva.py levas/minha_leva.toml --listar     # só mostra o plano (nada é criado)
    python3 tools/leva.py levas/minha_leva.toml --so-graficos  # não roda; refaz gráficos/tabela

O arquivo define: pasta de destino, threads, parâmetros fixos ([base]), grade de
parâmetros ([grade], produto cartesiano), rodadas avulsas ([[rodada]]), gráficos por
run, comparações ([[convergencia]], [[sobreposicao]]) e a tabela de resultados.
Modelo comentado: levas/exemplo.toml.

Cada rodada é feita por tools/run.py (com --skip-existing): relançar a mesma leva só
roda o que falta. Saídas da leva (tabela, comparações, cópia do .toml) ficam em
<pasta>/leva_<nome>/. Só usa a biblioteca padrão (Python >= 3.11, por causa do tomllib).
"""
import argparse
import csv
import datetime
import itertools
import os
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    import tomllib
except ImportError:  # Python < 3.11
    print('[leva] ERRO: precisa de Python >= 3.11 (módulo tomllib)', file=sys.stderr)
    sys.exit(1)

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FERR = os.path.join(RAIZ, 'tools')
sys.path.insert(0, FERR)
import run as R            # noqa: E402
import plot_run as PR      # noqa: E402
import plot_compare as PC  # noqa: E402

PARAMS = list(R.TESTE_PADRAO) + list(R.PLOT_PADRAO)
CHAVES = {
    'leva': {'nome', 'pasta', 'threads', 'simultaneos', 'pos_processamento', 'tag'},
    'graficos': {'por_run', 'png', 'ponto_decimal'},
    'convergencia': {'nome', 'param', 'y', 'erro_log', 'ref', 'eixo_x_log', 'filtro'},
    'sobreposicao': {'nome', 'variar', 'grandezas', 'legenda', 'xi_medio', 'filtro'},
    'tabela': {'colunas'},
}
SECOES = set(CHAVES) | {'base', 'grade', 'rodada'}
COLUNAS_TABELA = ['lambda', 'alfa', 'xi_medio', 'norma_u', 'tempo_s']


def erro(msg):
    print('[leva] ERRO: ' + msg, file=sys.stderr)
    sys.exit(1)


def aviso(msg):
    print('[leva] AVISO: ' + msg, file=sys.stderr)


def fv(v):
    """Valor curto para nomes de arquivo e tabelas."""
    if isinstance(v, bool):
        return 'T' if v else 'F'
    if isinstance(v, (int, float)):
        return '%g' % v
    return str(v)


# ----------------------------------------------------------------------------
# Leitura e validação do arquivo
# ----------------------------------------------------------------------------
def valida_chaves(onde, d, permitidas):
    for k in d:
        if k not in permitidas:
            erro('%s: chave desconhecida "%s" (válidas: %s)' % (onde, k, ', '.join(sorted(permitidas))))


def valida_param(onde, k):
    if k not in PARAMS:
        erro('%s: parâmetro desconhecido "%s"\n      válidos: %s' % (onde, k, ', '.join(PARAMS)))


def fmt_set(k, v):
    """Valor do TOML -> texto de --set do run.py (com checagem de tipo)."""
    if k in R.LOGICOS:
        if not isinstance(v, bool):
            erro('%s deve ser true/false (recebido %r)' % (k, v))
        return 'T' if v else 'F'
    if k in R.TEXTOS:
        if not isinstance(v, str):
            erro('%s deve ser texto entre aspas (recebido %r)' % (k, v))
        return v
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        erro('%s deve ser número (recebido %r)' % (k, v))
    if k in R.INTEIROS:
        if float(v) != int(v):
            erro('%s deve ser inteiro (recebido %r)' % (k, v))
        return str(int(v))
    return repr(float(v))


def le_config(caminho):
    try:
        with open(caminho, 'rb') as f:
            cfg = tomllib.load(f)
    except OSError as e:
        erro(str(e))
    except tomllib.TOMLDecodeError as e:
        erro('%s: erro de sintaxe TOML: %s' % (caminho, e))
    for s in cfg:
        if s not in SECOES:
            erro('seção desconhecida [%s] (válidas: %s)' % (s, ', '.join(sorted(SECOES))))
    for s in ('leva', 'graficos', 'tabela'):
        valida_chaves('[%s]' % s, cfg.get(s, {}), CHAVES[s])
    for s in ('convergencia', 'sobreposicao'):
        if s in cfg and not isinstance(cfg[s], list):
            erro('use [[%s]] (colchete duplo) para cada comparação' % s)
        for i, c in enumerate(cfg.get(s, []), 1):
            valida_chaves('[[%s]] nº %d' % (s, i), c, CHAVES[s])
    if 'rodada' in cfg and not isinstance(cfg['rodada'], list):
        erro('use [[rodada]] (colchete duplo) para cada rodada avulsa')
    for k in cfg.get('base', {}):
        valida_param('[base]', k)

    lv = cfg.setdefault('leva', {})
    lv.setdefault('nome', os.path.splitext(os.path.basename(caminho))[0])
    if 'pasta' not in lv:
        erro('[leva] precisa de "pasta" (destino dos runs, ex.: "Testes/minha_leva")')
    pasta = lv['pasta']
    lv['pasta'] = os.path.normpath(pasta if os.path.isabs(pasta) else os.path.join(RAIZ, pasta))
    lv.setdefault('threads', int(os.environ.get('OMP_NUM_THREADS', '8')))
    lv.setdefault('simultaneos', 1)
    lv.setdefault('pos_processamento', True)
    lv.setdefault('tag', '')
    gr = cfg.setdefault('graficos', {})
    gr.setdefault('por_run', False)
    gr.setdefault('png', True)
    gr.setdefault('ponto_decimal', False)
    tb = cfg.setdefault('tabela', {})
    tb.setdefault('colunas', COLUNAS_TABELA)
    for c in tb['colunas']:
        if c not in R.COLUNAS_CSV:
            erro('[tabela] coluna desconhecida "%s" (válidas: %s)' % (c, ', '.join(R.COLUNAS_CSV)))

    for i, c in enumerate(cfg.get('convergencia', []), 1):
        onde = '[[convergencia]] nº %d' % i
        p = c.get('param')
        if p is None:
            erro('%s: falta "param"' % onde)
        if p != 'razao' and p not in R.COLUNAS_CSV:
            erro('%s: param "%s" não é coluna do resumo (use um parâmetro de &teste ou "razao")' % (onde, p))
        for y in _lista(c.get('y', ['lambda'])):
            if y not in PC.Y_ROTULOS:
                erro('%s: y "%s" inválido (válidos: %s)' % (onde, y, ', '.join(PC.Y_ROTULOS)))
        valida_filtro(onde, c.get('filtro'))
    for i, c in enumerate(cfg.get('sobreposicao', []), 1):
        onde = '[[sobreposicao]] nº %d' % i
        p = c.get('variar')
        if p is None:
            erro('%s: falta "variar"' % onde)
        if p != 'razao':
            valida_param(onde, p)
        for q in _lista(c.get('grandezas', ['u'])):
            if q not in PR.QUANTIDADES:
                erro('%s: grandeza "%s" inválida (válidas: %s)' % (onde, q, ', '.join(PR.QUANTIDADES)))
        valida_filtro(onde, c.get('filtro'))

    precisa_derived = gr['por_run'] or cfg.get('sobreposicao') or any(
        y in ('xi_medio', 'norma_u') for c in cfg.get('convergencia', []) for y in _lista(c.get('y', [])))
    if precisa_derived and not lv['pos_processamento']:
        erro('gráficos por run, [[sobreposicao]] e y = xi_medio/norma_u precisam de pos_processamento = true')
    return cfg


def _lista(v):
    return v if isinstance(v, list) else [v]


def le_grade(grade):
    """[grade] -> lista de dimensões (chaves, [tuplas de valores]).
    "A=B" = [..]      mesmo valor para A e B
    "A,B" = [[a,b]..] pares (A, B) que andam juntos"""
    dims, vistas = [], set()
    for chave, valores in grade.items():
        if not isinstance(valores, list) or not valores:
            erro('[grade] %s: use uma lista não vazia, ex.: %s = [1, 2]' % (chave, chave))
        if '=' in chave:
            ks = [k.strip() for k in chave.split('=')]
            tuplas = [tuple([v] * len(ks)) for v in valores]
        elif ',' in chave:
            ks = [k.strip() for k in chave.split(',')]
            tuplas = []
            for v in valores:
                if not isinstance(v, list) or len(v) != len(ks):
                    erro('[grade] "%s": cada item deve ser uma lista com %d valores, ex.: [3.0, 16]'
                         % (chave, len(ks)))
                tuplas.append(tuple(v))
        else:
            ks, tuplas = [chave], [(v,) for v in valores]
        for k in ks:
            valida_param('[grade]', k)
            if k in vistas:
                erro('[grade]: "%s" aparece em mais de uma linha' % k)
            vistas.add(k)
        dims.append((ks, tuplas))
    return dims


def monta_rodadas(cfg):
    base = cfg.get('base', {})
    dims = le_grade(cfg.get('grade', {}))
    lista = []
    if dims:
        for combo in itertools.product(*[d[1] for d in dims]):
            o = dict(base)
            for (ks, _), vals in zip(dims, combo):
                o.update(zip(ks, vals))
            lista.append(o)
    for i, r in enumerate(cfg.get('rodada', []), 1):
        for k in r:
            valida_param('[[rodada]] nº %d' % i, k)
        o = dict(base)
        o.update(r)
        lista.append(o)
    if not lista:
        lista = [dict(base)]

    rodadas, por_id = [], {}
    for o in lista:
        sets = ['%s=%s' % (k, fmt_set(k, v)) for k, v in o.items()]
        t, p = dict(R.TESTE_PADRAO), dict(R.PLOT_PADRAO)
        for s in sets:
            k, v = s.split('=', 1)
            (t if k in t else p)[k] = R.converte(k, v)
        if t['Nu'] < 0:
            t['Nu'] = t['Nv']
        rid = R.monta_run_id(t, cfg['leva']['tag'] or None)
        efetivo = dict(t)
        efetivo.update(p)
        if rid in por_id:
            if por_id[rid]['params'] != efetivo:
                erro('duas rodadas com o mesmo run_id e parâmetros diferentes:\n      %s\n'
                     '      (parâmetros de &plot não entram no run_id: não os varie na mesma leva)' % rid)
            continue
        r = dict(sets=sets, run_id=rid, params=efetivo,
                 run_dir=os.path.join(cfg['leva']['pasta'], rid))
        por_id[rid] = r
        rodadas.append(r)
    return rodadas, dims


# ----------------------------------------------------------------------------
# Execução
# ----------------------------------------------------------------------------
def le_chave_valor(caminho):
    d = {}
    if os.path.isfile(caminho):
        with open(caminho, encoding='utf-8', errors='replace') as f:
            for linha in f:
                if '=' in linha and not linha.startswith('#'):
                    k, v = linha.split('=', 1)
                    d[k.strip()] = v.strip()
    return d


def estado(r):
    info = le_chave_valor(os.path.join(r['run_dir'], 'input', 'run_info.txt'))
    if not info:
        return 'falta' if not os.path.exists(r['run_dir']) else 'falhou'
    if info.get('rc_spline') == '0' and info.get('rc_plot') in ('0', 'None'):
        return 'ok'
    return 'falhou'


_SAIDA = threading.Lock()


def diz(msg):
    """print protegido: com simultaneos > 1 as linhas de rodadas diferentes não se misturam."""
    with _SAIDA:
        print(msg, flush=True)


class Acompanha:
    """Repete no terminal as linhas novas de um log (spline.log, plot.log) enquanto ele cresce.
    Ignora um log antigo (de uma tentativa anterior) até ele ser reescrito pela rodada atual."""

    def __init__(self, caminho, prefixo, t0):
        self.caminho, self.prefixo, self.t0 = caminho, prefixo, t0
        self.pos, self.buf = 0, ''

    def le(self, final=False):
        try:
            if os.path.getmtime(self.caminho) < self.t0 - 1:
                return
            with open(self.caminho, 'rb') as f:
                f.seek(self.pos)
                dados = f.read()
        except OSError:
            return
        self.pos += len(dados)
        self.buf += dados.decode('utf-8', 'replace')
        linhas = self.buf.split('\n')
        self.buf = linhas.pop()
        if final and self.buf:
            linhas.append(self.buf)
            self.buf = ''
        for l in linhas:
            if l.strip():
                diz(self.prefixo + l.rstrip())


def roda_uma(r, lv, prefixo=''):
    cmd = [sys.executable, os.path.join(FERR, 'run.py'),
           '--teste', os.path.basename(lv['pasta']), '--raiz-testes', os.path.dirname(lv['pasta']),
           '--threads', str(lv['threads']), '--skip-existing']
    if not lv['pos_processamento']:
        cmd.append('--no-plot')
    if lv['tag']:
        cmd += ['--tag', lv['tag']]
    cmd += ['--set'] + r['sets']
    t0 = time.time()
    logs = [Acompanha(os.path.join(r['run_dir'], 'logs', n), prefixo, t0) for n in ('spline.log', 'plot.log')]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    while True:
        try:
            out, err = proc.communicate(timeout=1.0)
            break
        except subprocess.TimeoutExpired:
            for lg in logs:
                lg.le()
    for lg in logs:
        lg.le(final=True)
    return (out.strip() or err.strip()).splitlines()[-1:] or ['(sem saída)']


def roda_tudo(rodadas, lv):
    falta = [r for r in rodadas if estado(r) == 'falta']
    print('%d rodadas na leva, %d a rodar (%d threads cada, %d simultânea(s))'
          % (len(rodadas), len(falta), lv['threads'], lv['simultaneos']), flush=True)
    if lv['threads'] * lv['simultaneos'] > (os.cpu_count() or 1):
        aviso('threads x simultaneos = %d > %d núcleos' % (lv['threads'] * lv['simultaneos'], os.cpu_count()))
    n = len(falta)

    def tarefa(i, r):
        # com várias rodadas simultâneas, cada linha do log leva o número da rodada
        prefixo = '   [#%d] ' % i if lv['simultaneos'] > 1 else '   '
        diz('[#%d/%d] início %s  (%s)' % (i, n, r['run_id'], datetime.datetime.now().strftime('%H:%M:%S')))
        return roda_uma(r, lv, prefixo)

    feitas = 0
    with ThreadPoolExecutor(max_workers=max(1, lv['simultaneos'])) as ex:
        futuros = {ex.submit(tarefa, i, r): (i, r) for i, r in enumerate(falta, 1)}
        for fut in as_completed(futuros):
            feitas += 1
            diz('[#%d/%d] fim (%d de %d concluídas): %s'
                % (futuros[fut][0], n, feitas, n, fut.result()[0][:220]))


def env_graficos(gr):
    env = dict(os.environ)
    if gr['ponto_decimal']:
        env['LC_ALL'] = 'C'
    return env


def graficos_por_run(rodadas, gr, log):
    env = env_graficos(gr)
    feitos = 0
    for r in rodadas:
        if estado(r) != 'ok' or not os.path.isdir(os.path.join(r['run_dir'], 'derived')):
            continue
        d_agr, d_png = (os.path.join(r['run_dir'], 'plots', s) for s in ('agr', 'png'))
        tem = lambda d, ext: os.path.isdir(d) and any(n.endswith(ext) for n in os.listdir(d))  # noqa: E731
        if tem(d_agr, '.agr') and (not gr['png'] or tem(d_png, '.png')):
            continue
        cmd = [sys.executable, os.path.join(FERR, 'plot_run.py'), r['run_dir']] + (['--png'] if gr['png'] else [])
        s = subprocess.run(cmd, capture_output=True, text=True, env=env)
        log.write('$ %s\n%s%s\n' % (' '.join(cmd), s.stdout, s.stderr))
        feitos += 1
    return feitos


# ----------------------------------------------------------------------------
# Resultados, comparações e tabela
# ----------------------------------------------------------------------------
def linha_resultado(r):
    obs = le_chave_valor(os.path.join(r['run_dir'], 'derived', 'observaveis.txt'))
    info = le_chave_valor(os.path.join(r['run_dir'], 'input', 'run_info.txt'))
    linha = {k: obs.get(k, '') for k in R.COLUNAS_CSV}
    linha.update(run_id=r['run_id'], teste=os.path.basename(os.path.dirname(r['run_dir'])),
                 threads=info.get('threads', ''), tempo_parede_s=info.get('tempo_parede_s', ''),
                 data=info.get('data', ''))
    try:
        linha['razao'] = repr(float(linha['m2']) / float(linha['m1']))
    except (ValueError, ZeroDivisionError):
        linha['razao'] = ''
    return linha


COLS_SAIDA = R.COLUNAS_CSV + ['razao']


def escreve_csv(caminho, linhas):
    with open(caminho, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=COLS_SAIDA, extrasaction='ignore')
        w.writeheader()
        w.writerows(linhas)


def valor_param(r, p):
    if p == 'razao':
        return r['params']['m2'] / r['params']['m1']
    return r['params'][p]


def chaves_de(p, dims):
    """Chaves que andam junto com p (dimensão ligada da grade)."""
    if p == 'razao':
        return {'m1', 'm2', 'razao'}
    for ks, _ in dims:
        if p in ks:
            return set(ks)
    return {p}


def valida_filtro(onde, filtro):
    if filtro is None:
        return
    if not isinstance(filtro, dict):
        erro('%s: filtro deve ser uma tabela, ex.: filtro = { zeta = 1.0 }' % onde)
    for k in filtro:
        valida_param(onde + ' filtro', k)


def filtra(rodadas, filtro, onde):
    if not filtro:
        return rodadas

    def igual(a, b):
        if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(b, bool):
            return abs(a - b) <= 1e-9 * max(1.0, abs(b))
        return a == b
    return [r for r in rodadas if all(igual(r['params'][k], v) for k, v in filtro.items())]


def agrupa(rodadas, p, dims):
    """Separa as rodadas por valores dos OUTROS parâmetros que variam; ordena cada grupo por p."""
    if not rodadas:
        return []
    ligadas = chaves_de(p, dims)
    variam = [k for k in rodadas[0]['params']
              if len({repr(r['params'][k]) for r in rodadas}) > 1 and k not in ligadas]
    grupos = {}
    for r in rodadas:
        chave = tuple((k, r['params'][k]) for k in variam)
        grupos.setdefault(chave, []).append(r)
    # só grupos com >= 2 rodadas viram gráfico; o nome leva só o que os distingue entre si
    uteis = [(chave, g) for chave, g in grupos.items() if len(g) >= 2]
    distintos = {k for i, k in enumerate(variam) if len({repr(c[i][1]) for c, _ in uteis}) > 1}
    saida = []
    for chave, g in uteis:
        g.sort(key=lambda r: valor_param(r, p))
        sufixo = '_'.join('%s-%s' % (k, fv(v)) for k, v in chave if k in distintos)
        saida.append((sufixo, g))
    return saida


def chama(cmd, env, log, gerados):
    s = subprocess.run(cmd, capture_output=True, text=True, env=env)
    log.write('$ %s\n%s%s\n' % (' '.join(cmd), s.stdout, s.stderr))
    novos = [l for l in s.stdout.splitlines() if l.endswith(('.agr', '.png'))]
    gerados.extend(novos)
    if s.returncode != 0:
        alvo = cmd[cmd.index('--out') + 1] if '--out' in cmd else ' '.join(cmd[1:3])
        msg = (s.stderr.strip().splitlines() or ['?'])[-1]
        if 'nenhuma série com dados válidos' in msg:
            msg = 'sem dados para plotar (ex.: |Δy| = 0 em todos os pontos no gráfico de erro em log)'
        aviso('%s não gerado: %s' % (os.path.basename(alvo), msg))
    return s.returncode


def comparacoes(cfg, rodadas, dims, dir_leva, linhas, log):
    gr = cfg['graficos']
    env = env_graficos(gr)
    png = ['--png'] if gr['png'] else []
    gerados = []
    ok = [r for r in rodadas if r['run_id'] in linhas]
    if not ok:
        return gerados
    pasta_csv = os.path.join(dir_leva, 'csv')

    for i, c in enumerate(cfg.get('convergencia', []), 1):
        onde = '[[convergencia]] nº %d' % i
        p = c.get('param')
        if p is None:
            erro('%s: falta "param"' % onde)
        if p != 'razao' and p not in R.COLUNAS_CSV:
            erro('%s: param "%s" não é coluna do resumo (use um parâmetro de &teste ou "razao")' % (onde, p))
        ys = _lista(c.get('y', ['lambda']))
        for y in ys:
            if y not in PC.Y_ROTULOS:
                erro('%s: y "%s" inválido (válidos: %s)' % (onde, y, ', '.join(PC.Y_ROTULOS)))
        nome = c.get('nome', 'conv_%s' % p)
        sel = filtra(ok, c.get('filtro'), onde)
        for sufixo, g in agrupa(sel, p, dims) if sel else []:
            if len(g) < 2:
                continue
            base = nome + ('_' + sufixo if sufixo else '')
            os.makedirs(pasta_csv, exist_ok=True)
            sub = os.path.join(pasta_csv, base + '.csv')
            escreve_csv(sub, [linhas[r['run_id']] for r in g])
            for y in ys:
                cmd = [sys.executable, os.path.join(FERR, 'plot_compare.py'), 'conv', '--csv', sub,
                       '--param', p, '--y', y, '--out', os.path.join(dir_leva, '%s_%s.agr' % (base, y))]
                if c.get('eixo_x_log'):
                    cmd.append('--logx')
                chama(cmd + png, env, log, gerados)
                if c.get('erro_log'):
                    cmd = [sys.executable, os.path.join(FERR, 'plot_compare.py'), 'conv', '--csv', sub,
                           '--param', p, '--y', y, '--logerr',
                           '--out', os.path.join(dir_leva, '%s_%s_err.agr' % (base, y))]
                    if c.get('ref') is not None:
                        cmd += ['--ref', repr(float(c['ref']))]
                    if c.get('eixo_x_log'):
                        cmd.append('--logx')
                    chama(cmd + png, env, log, gerados)

    for i, c in enumerate(cfg.get('sobreposicao', []), 1):
        onde = '[[sobreposicao]] nº %d' % i
        p = c.get('variar')
        if p is None:
            erro('%s: falta "variar"' % onde)
        if p != 'razao':
            valida_param(onde, p)
        qs = _lista(c.get('grandezas', ['u']))
        for q in qs:
            if q not in PR.QUANTIDADES:
                erro('%s: grandeza "%s" inválida (válidas: %s)' % (onde, q, ', '.join(PR.QUANTIDADES)))
        legenda = c.get('legenda', p)
        nome = c.get('nome', 'sobre_%s' % p)
        sel = filtra(ok, c.get('filtro'), onde)
        for sufixo, g in agrupa(sel, p, dims) if sel else []:
            if len(g) < 2:
                continue
            base = nome + ('_' + sufixo if sufixo else '')
            for q in qs:
                cmd = [sys.executable, os.path.join(FERR, 'plot_compare.py'), 'overlay',
                       '--runs'] + [r['run_dir'] for r in g] + [
                       '--quantity', q, '--legend-key', legenda,
                       '--out', os.path.join(dir_leva, '%s_%s.agr' % (base, q))]
                if q == 'u' and c.get('xi_medio', True):
                    cmd.append('--vlines')
                chama(cmd + png, env, log, gerados)
    return gerados


def relatorio(cfg, caminho_cfg, rodadas, dims, linhas, gerados, dir_leva):
    cols = cfg['tabela']['colunas']
    variam = [k for k in rodadas[0]['params'] if len({repr(r['params'][k]) for r in rodadas}) > 1]
    L = ['# Leva `%s`' % cfg['leva']['nome'], '',
         '- gerado por `tools/leva.py` em %s' % datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
         '- configuração: `%s` (cópia nesta pasta)' % (
             os.path.relpath(os.path.abspath(caminho_cfg), RAIZ)
             if os.path.abspath(caminho_cfg).startswith(RAIZ + os.sep) else os.path.abspath(caminho_cfg)),
         '- runs em: `%s`' % cfg['leva']['pasta'],
         '- rodadas: %d (%d ok)' % (len(rodadas), len(linhas)), '']
    falhas = [r['run_id'] for r in rodadas if r['run_id'] not in linhas]
    if falhas:
        L += ['## Rodadas que falharam ou faltam', ''] + ['- `%s`' % f for f in falhas] + ['']
    L += ['## Resultados', '',
          '| ' + ' | '.join(variam + cols + ['run_id']) + ' |',
          '|' + '---|' * (len(variam) + len(cols) + 1)]
    for r in rodadas:
        l = linhas.get(r['run_id'])
        if l is None:
            continue
        cel = [fv(r['params'][k]) for k in variam]
        for c in cols:
            try:
                cel.append('%.10g' % float(l[c]))
            except (ValueError, TypeError):
                cel.append(l.get(c, '') or '—')
        L.append('| ' + ' | '.join(cel + ['`%s`' % r['run_id']]) + ' |')
    if gerados:
        L += ['', '## Gráficos de comparação', ''] + sorted(
            '- `%s`' % os.path.relpath(g, dir_leva) for g in gerados if os.path.dirname(os.path.abspath(g)) == dir_leva)
    with open(os.path.join(dir_leva, 'relatorio.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('config', help='arquivo .toml da leva')
    ap.add_argument('--listar', action='store_true', help='só mostra o plano; não cria nada')
    ap.add_argument('--so-graficos', action='store_true', help='não roda nada; refaz gráficos, CSV e relatório')
    ap.add_argument('--threads', type=int, help='sobrepõe [leva] threads')
    ap.add_argument('--simultaneos', type=int, help='sobrepõe [leva] simultaneos')
    a = ap.parse_args()

    cfg = le_config(a.config)
    lv, gr = cfg['leva'], cfg['graficos']
    if a.threads:
        lv['threads'] = a.threads
    if a.simultaneos:
        lv['simultaneos'] = a.simultaneos
    rodadas, dims = monta_rodadas(cfg)
    dir_leva = os.path.join(lv['pasta'], 'leva_' + lv['nome'])

    if a.listar:
        print('leva    : %s' % lv['nome'])
        print('pasta   : %s' % lv['pasta'])
        print('saídas  : %s' % dir_leva)
        print('rodadas : %d   (threads %d, simultâneas %d, pós-processamento %s)'
              % (len(rodadas), lv['threads'], lv['simultaneos'], 'sim' if lv['pos_processamento'] else 'não'))
        for i, r in enumerate(rodadas, 1):
            print('%4d  %-7s %s' % (i, estado(r), r['run_id']))
        print('gráficos por run: %s   png: %s   ponto decimal: %s'
              % (gr['por_run'], gr['png'], gr['ponto_decimal']))
        for s, campo in (('convergencia', 'param'), ('sobreposicao', 'variar')):
            for c in cfg.get(s, []):
                p = c.get(campo)
                ng = len([g for _, g in agrupa(filtra(rodadas, c.get('filtro'), s), p, dims) if len(g) >= 2]) \
                    if p in PARAMS or p == 'razao' else 0
                print('[[%s]] %s = %s: %d grupo(s) de >= 2 rodadas' % (s, campo, p, ng))
        print('(listar: nada foi criado)')
        return 0

    if not a.so_graficos:
        os.makedirs(lv['pasta'], exist_ok=True)
        roda_tudo(rodadas, lv)

    os.makedirs(dir_leva, exist_ok=True)
    copia = os.path.join(dir_leva, '%s_%s.toml' % (lv['nome'], datetime.datetime.now().strftime('%Y%m%d-%H%M%S')))
    shutil.copy2(a.config, copia)

    linhas = {r['run_id']: linha_resultado(r) for r in rodadas if estado(r) == 'ok'}
    escreve_csv(os.path.join(dir_leva, 'rodadas.csv'), [linhas[r['run_id']] for r in rodadas if r['run_id'] in linhas])
    with open(os.path.join(dir_leva, 'graficos.log'), 'w', encoding='utf-8') as log:
        if gr['por_run']:
            n = graficos_por_run(rodadas, gr, log)
            print('gráficos por run: %d run(s) processado(s)' % n)
        gerados = comparacoes(cfg, rodadas, dims, dir_leva, linhas, log)
    relatorio(cfg, a.config, rodadas, dims, linhas, gerados, dir_leva)

    falhas = [r for r in rodadas if r['run_id'] not in linhas]
    print('%d/%d rodadas ok; %d arquivo(s) de comparação; relatório: %s'
          % (len(linhas), len(rodadas), len(gerados), os.path.join(dir_leva, 'relatorio.md')))
    for r in falhas:
        print('   FALHOU/FALTA: %s (ver %s/logs/)' % (r['run_id'], r['run_dir']))
    return 2 if falhas else 0


if __name__ == '__main__':
    sys.exit(main())
