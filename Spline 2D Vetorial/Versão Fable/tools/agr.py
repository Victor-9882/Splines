#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""agr.py -- gerador de arquivos Xmgrace (.agr) no estilo do template
`template_plot_long.agr`.  Python 3 puro (apenas biblioteca padrão).

API principal:
    write_agr(path, series, xlabel, ylabel, world=None, vlines=None, title=None,
              logy=False, logx=False, legend_pos=(0.8, 0.8))

    series : lista de (x, y, legenda, estilo)
             estilo = None (padrão por índice) ou dict com chaves opcionais
             color, linestyle, linewidth, symbol, symbol_size, comment
    vlines : lista de (x, estilo, legenda) -- retas verticais (2 pontos)

Auxiliares: read_columns, nice_world, nice_tick.

Uso em linha de comando:
    python3 agr.py saida.agr dados1.dat [dados2.dat ...] --xlabel ... --ylabel ... \
        --legend ... --legend ...
"""
import argparse
import math
import re
import sys
import time

# --------------------------------------------------------------------------
# Constantes de estilo (copiadas do cabeçalho do template)
# --------------------------------------------------------------------------
FONTES = [
    (8, "Courier"), (10, "Courier-Bold"), (11, "Courier-BoldOblique"),
    (9, "Courier-Oblique"), (4, "Helvetica"), (6, "Helvetica-Bold"),
    (7, "Helvetica-BoldOblique"), (5, "Helvetica-Oblique"),
    (14, "NimbusMonoPS-BoldItalic"), (15, "NimbusMonoPS-Italic"),
    (16, "NimbusMonoPS-Regular"), (17, "NimbusRoman-Bold"),
    (18, "NimbusRoman-BoldItalic"), (19, "NimbusRoman-Italic"),
    (20, "NimbusRoman-Regular"), (21, "NimbusSans-Bold"),
    (22, "NimbusSans-BoldItalic"), (23, "NimbusSans-Italic"),
    (24, "NimbusSans-Regular"), (25, "NimbusSansNarrow-Bold"),
    (26, "NimbusSansNarrow-BoldOblique"), (27, "NimbusSansNarrow-Oblique"),
    (28, "NimbusSansNarrow-Regular"), (29, "P052-Bold"),
    (30, "P052-BoldItalic"), (31, "P052-Italic"), (32, "P052-Roman"),
    (33, "StandardSymbolsPS-Regular"), (12, "Symbol"), (2, "Times-Bold"),
    (3, "Times-BoldItalic"), (1, "Times-Italic"), (0, "Times-Roman"),
    (39, "URWBookman-Demi"), (40, "URWBookman-DemiItalic"),
    (41, "URWBookman-Light"), (42, "URWBookman-LightItalic"),
    (43, "URWGothic-Book"), (44, "URWGothic-BookOblique"),
    (45, "URWGothic-Demi"), (46, "URWGothic-DemiOblique"),
    (47, "Z003-MediumItalic"), (13, "ZapfDingbats"),
]

CORES = [
    (0, (255, 255, 255), "white"), (1, (0, 0, 0), "black"),
    (2, (255, 0, 0), "red"), (3, (0, 255, 0), "green"),
    (4, (0, 0, 255), "blue"), (5, (255, 255, 0), "yellow"),
    (6, (188, 143, 143), "brown"), (7, (220, 220, 220), "grey"),
    (8, (148, 0, 211), "violet"), (9, (0, 255, 255), "cyan"),
    (10, (255, 0, 255), "magenta"), (11, (255, 165, 0), "orange"),
    (12, (114, 33, 188), "indigo"), (13, (103, 7, 72), "maroon"),
    (14, (64, 224, 208), "turquoise"), (15, (0, 139, 0), "green4"),
    (16, (0, 0, 205), "azul"), (17, (178, 34, 34), "vermelho-tijolo"),
    (18, (100, 100, 100), "cinza"),
]

COR_AZUL, COR_TIJOLO, COR_CINZA = 16, 17, 18
# cores "padrão do Grace" usadas a partir da série 3
CORES_EXTRAS = [4, 15, 11, 12, 13, 14, 10, 2, 8, 9]

SIMBOLOS = {"none": 0, "circle": 1, "square": 2, "diamond": 3,
            "triangle": 4, "triangle_up": 4, "triangle_left": 5,
            "triangle_down": 6, "triangle_right": 7, "plus": 8, "x": 9,
            "star": 10}


def estilo_padrao(i):
    """Estilo padrão da série de índice i (template: 0 azul tracejada,
    1 vermelho-tijolo contínua, 2 cinza, depois cores padrão do Grace)."""
    if i == 0:
        return dict(color=COR_AZUL, linestyle=3, linewidth=2.0)
    if i == 1:
        return dict(color=COR_TIJOLO, linestyle=1, linewidth=1.5)
    if i == 2:
        return dict(color=COR_CINZA, linestyle=1, linewidth=2.0)
    return dict(color=CORES_EXTRAS[(i - 3) % len(CORES_EXTRAS)],
                linestyle=1, linewidth=2.0)


def _aviso(msg):
    print("[agr] AVISO: " + msg, file=sys.stderr)


# --------------------------------------------------------------------------
# Leitura de dados
# --------------------------------------------------------------------------
_RE_FORTRAN_SEM_E = re.compile(r"^([+-]?\d*\.?\d+)([+-]\d{3})$")


def _para_float(tok):
    """Converte token numérico Fortran (D, E, expoente de 3 dígitos, ou
    expoente sem letra '1.0-100') em float.  Levanta ValueError."""
    t = tok.strip().replace("D", "E").replace("d", "e")
    try:
        return float(t)
    except ValueError:
        m = _RE_FORTRAN_SEM_E.match(t)
        if m:
            return float(m.group(1) + "E" + m.group(2))
        raise


def read_columns(path, cols=(0, 1)):
    """Lê colunas numéricas de um arquivo texto.
    Ignora linhas vazias e comentários (# ou !).  Linhas ilegíveis ou
    curtas são puladas com aviso.  Retorna lista de listas (uma por coluna).
    NaN/Inf são mantidos aqui (write_agr os descarta com aviso)."""
    cols = tuple(cols)
    saida = [[] for _ in cols]
    ruins = 0
    with open(path, "r", errors="replace") as f:
        for n, linha in enumerate(f, 1):
            s = linha.strip()
            if not s or s[0] in "#!":
                continue
            toks = s.replace(",", " ").split()
            try:
                if max(cols) >= len(toks):
                    raise ValueError("colunas insuficientes")
                vals = [_para_float(toks[c]) for c in cols]
            except ValueError:
                ruins += 1
                if ruins <= 3:
                    _aviso("%s:%d linha ilegível ignorada: %r" % (path, n, s[:60]))
                continue
            for k, v in enumerate(vals):
                saida[k].append(v)
    if ruins > 3:
        _aviso("%s: %d linhas ilegíveis ignoradas no total" % (path, ruins))
    return saida


# --------------------------------------------------------------------------
# Mundo e ticks "bonitos"
# --------------------------------------------------------------------------
def nice_tick(span, n=5):
    """Espaçamento 'bonito' (1, 2, 5 x 10^k) para ~n intervalos em span."""
    if not (span > 0) or math.isinf(span) or math.isnan(span):
        return 1.0
    bruto = span / float(n)
    mag = 10.0 ** math.floor(math.log10(bruto))
    f = bruto / mag
    if f < 1.5:
        b = 1.0
    elif f < 3.5:
        b = 2.0
    elif f < 7.5:
        b = 5.0
    else:
        b = 10.0
    return b * mag


def _arred(v, tick):
    """Remove ruído de ponto flutuante de múltiplos de tick."""
    return float("%.12g" % v)


def _eixo_linear(vmin, vmax, base_zero=False, folga=0.0):
    if vmax == vmin:
        d = abs(vmax) * 0.1 if vmax != 0 else 1.0
        vmin, vmax = vmin - d, vmax + d
    span = vmax - vmin
    if base_zero and vmin >= 0 and vmin < 0.5 * vmax:
        vmin = 0.0
        span = vmax - vmin
    vmax = vmax + folga * span
    tick = nice_tick(vmax - vmin)
    lo = math.floor(vmin / tick + 1e-9) * tick
    hi = math.ceil(vmax / tick - 1e-9) * tick
    tick = nice_tick(hi - lo)
    # recalcula limites com o tick final (garante múltiplos)
    lo = math.floor(vmin / tick + 1e-9) * tick
    hi = math.ceil(vmax / tick - 1e-9) * tick
    return _arred(lo, tick), _arred(hi, tick), _arred(tick, tick)


def _eixo_log(vmin, vmax):
    lo = math.floor(math.log10(vmin) + 1e-9)
    hi = math.ceil(math.log10(vmax) - 1e-9)
    if hi <= lo:
        hi = lo + 1
    return 10.0 ** lo, 10.0 ** hi


def nice_world(xs, ys, logx=False, logy=False):
    """Calcula (xmin, ymin, xmax, ymax, xtick, ytick) a partir dos dados.
    xs, ys: listas planas de números (ou listas de listas).  Para eixos
    log os ticks retornados são 10 (década)."""
    def achata(a):
        out = []
        for e in a:
            if isinstance(e, (list, tuple)):
                out.extend(e)
            else:
                out.append(e)
        return [v for v in out if math.isfinite(v)]
    xs, ys = achata(xs), achata(ys)
    if not xs or not ys:
        raise ValueError("sem dados finitos para calcular o world")
    if logx:
        x0, x1 = _eixo_log(min(xs), max(xs))
        xt = 10.0
    else:
        x0, x1, xt = _eixo_linear(min(xs), max(xs))
    if logy:
        y0, y1 = _eixo_log(min(ys), max(ys))
        yt = 10.0
    else:
        y0, y1, yt = _eixo_linear(min(ys), max(ys), base_zero=True, folga=0.04)
    return x0, y0, x1, y1, xt, yt


# --------------------------------------------------------------------------
# Escrita do .agr
# --------------------------------------------------------------------------
def _fmt(v):
    return "%.10g" % v


def _q(s):
    """Aspas para strings do Grace (barra invertida é código do Grace)."""
    return '"' + str(s).replace('"', "'") + '"'


def _cabecalho_fixo():
    L = ["# Grace project file", "#",
         "@version 50125",
         "@page size 792, 612", "@page scroll 5%", "@page inout 5%",
         "@link page off"]
    for n, nome in FONTES:
        L.append('@map font %d to "%s", "%s"' % (n, nome, nome))
    for n, rgb, nome in CORES:
        L.append('@map color %d to (%d, %d, %d), "%s"' % ((n,) + rgb + (nome,)))
    L += ["@reference date 0", "@date wrap off", "@date wrap year 1950",
          "@default linewidth 1.5", "@default linestyle 1", "@default color 1",
          "@default pattern 1", "@default font 0",
          "@default char size 1.000000", "@default symbol size 1.000000",
          '@default sformat "%.8g"', "@background color 0",
          "@page background fill on", "@timestamp off",
          "@timestamp 0.03, 0.03", "@timestamp color 1", "@timestamp rot 0",
          "@timestamp font 0", "@timestamp char size 1.000000",
          '@timestamp def "%s"' % time.strftime("%a %b %e %H:%M:%S %Y")]
    for r in range(5):
        L += ["@r%d off" % r, "@link r%d to g0" % r, "@r%d type above" % r,
              "@r%d linestyle 1" % r, "@r%d linewidth 1.0" % r,
              "@r%d color 1" % r, "@r%d line 0, 0, 0, 0" % r]
    L += ["@g0 on", "@g0 hidden false", "@g0 type XY", "@g0 stacked false",
          "@g0 bar hgap 0.000000", "@g0 fixedpoint off",
          "@g0 fixedpoint type 0", "@g0 fixedpoint xy 0.000000, 0.000000",
          "@g0 fixedpoint format general general", "@g0 fixedpoint prec 6, 6"]
    return L


def _comprimento_visivel(label):
    """Nº aproximado de caracteres visíveis de um rótulo (sem códigos do Grace)."""
    t = re.sub(r"\\f\{[^}]*\}|\\#\{[0-9a-fA-F]+\}|\\[sSNcCxXbB]", "", label)
    return len(t)


def _prec_sig(lo, hi, tick):
    """Algarismos significativos para o ticklabel (formato general) distinguir
    os ticks vizinhos; mínimo 2 (como no template)."""
    m = max(abs(lo), abs(hi))
    if m <= 0 or tick <= 0:
        return 2
    sig = int(math.floor(math.log10(m)) - math.floor(math.log10(tick)) + 1)
    return max(2, min(sig, 8))


def _eixo(nome, label, major, minor_ticks, label_place, perp, log, prec=2):
    # rótulo longo no eixo y: layout 'perp' (horizontal) estouraria a margem
    if perp and _comprimento_visivel(label) > 7:
        perp = False
    p = "@    %saxis  " % nome
    L = [p + "on", p + "type zero false", p + "offset 0.000000 , 0.000000",
         p + "bar on", p + "bar color 1", p + "bar linestyle 1",
         p + "bar linewidth 1.5",
         p + "label " + _q(label),
         p + "label layout " + ("perp" if perp else "para"),
         p + "label place spec",
         p + "label place %s" % label_place,
         p + "label char size 1.750000", p + "label font 0",
         p + "label color 1", p + "label place normal",
         p + "tick on", p + "tick major " + _fmt(major),
         p + "tick minor ticks %d" % minor_ticks,
         p + "tick default 6", p + "tick place rounded true", p + "tick in",
         p + "tick major size 1.100000", p + "tick major color 1",
         p + "tick major linewidth 1.5", p + "tick major linestyle 1",
         p + "tick major grid off", p + "tick minor color 1",
         p + "tick minor linewidth 1.0", p + "tick minor linestyle 1",
         p + "tick minor grid off", p + "tick minor size 0.650000",
         p + "ticklabel on",
         p + "ticklabel format " + ("power" if log else "general"),
         p + "ticklabel prec %d" % (0 if log else prec), p + 'ticklabel formula ""',
         p + 'ticklabel append ""', p + 'ticklabel prepend ""',
         p + "ticklabel angle 0", p + "ticklabel skip 0",
         p + "ticklabel stagger 0", p + "ticklabel place normal",
         p + "ticklabel offset auto",
         p + "ticklabel offset 0.000000 , 0.010000",
         p + "ticklabel start type auto", p + "ticklabel start 0.000000",
         p + "ticklabel stop type auto", p + "ticklabel stop 0.000000",
         p + "ticklabel char size 1.300000", p + "ticklabel font 0",
         p + "ticklabel color 1", p + "tick place both",
         p + "tick spec type none"]
    return L


def _bloco_serie(i, est, legenda):
    c = est["color"]
    sim = est.get("symbol", 0)
    if isinstance(sim, str):
        sim = SIMBOLOS.get(sim, 0)
    p = "@    s%d " % i
    return [
        p + "hidden false", p + "type xy", p + "symbol %d" % sim,
        p + "symbol size %s" % ("%.6f" % est.get("symbol_size", 1.0)),
        p + "symbol color %d" % c, p + "symbol pattern 1",
        p + "symbol fill color %d" % c,
        p + "symbol fill pattern %d" % (1 if sim else 0),
        p + "symbol linewidth 1.5", p + "symbol linestyle 1",
        p + "symbol char 65", p + "symbol char font 0", p + "symbol skip 0",
        p + "line type 1", p + "line linestyle %d" % est["linestyle"],
        p + "line linewidth %s" % ("%.1f" % est["linewidth"]),
        p + "line color %d" % c, p + "line pattern 1",
        p + "baseline type 0", p + "baseline off", p + "dropline off",
        p + "fill type 0", p + "fill rule 0", p + "fill color 1",
        p + "fill pattern 1", p + "avalue off", p + "avalue type 2",
        p + "avalue char size 1.000000", p + "avalue font 0",
        p + "avalue color 1", p + "avalue rot 0",
        p + "avalue format general", p + "avalue prec 3",
        p + 'avalue prepend ""', p + 'avalue append ""',
        p + "avalue offset 0.000000 , 0.000000", p + "errorbar on",
        p + "errorbar place both", p + "errorbar color %d" % c,
        p + "errorbar pattern 1", p + "errorbar size 1.000000",
        p + "errorbar linewidth 1.5", p + "errorbar linestyle 1",
        p + "errorbar riser linewidth 1.5", p + "errorbar riser linestyle 1",
        p + "errorbar riser clip off", p + "errorbar riser clip length 0.100000",
        p + "comment " + _q(est.get("comment", "")),
        p + "legend  " + _q(legenda or ""),
    ]


def _limpa_serie(x, y, nome, logx, logy):
    """Descarta NaN/Inf (e <=0 em eixo log), com aviso em stderr."""
    if len(x) != len(y):
        n = min(len(x), len(y))
        _aviso("série '%s': x e y com tamanhos diferentes (%d, %d); usando %d"
               % (nome, len(x), len(y), n))
        x, y = x[:n], y[:n]
    xo, yo = [], []
    nnan = nlog = 0
    for a, b in zip(x, y):
        if not (math.isfinite(a) and math.isfinite(b)):
            nnan += 1
            continue
        if (logx and a <= 0) or (logy and b <= 0):
            nlog += 1
            continue
        xo.append(a)
        yo.append(b)
    if nnan:
        _aviso("série '%s': %d ponto(s) NaN/Inf descartado(s)" % (nome, nnan))
    if nlog:
        _aviso("série '%s': %d ponto(s) <= 0 descartado(s) (eixo log)" % (nome, nlog))
    return xo, yo


def write_agr(path, series, xlabel, ylabel, world=None, vlines=None,
              title=None, logy=False, logx=False, legend_pos=(0.8, 0.8)):
    """Escreve o .agr.  world = (xmin, ymin, xmax, ymax) ou None (automático;
    ticks sempre calculados).  Levanta ValueError se não houver nenhum dado
    válido."""
    vlines = vlines or []
    limpas = []
    for i, s in enumerate(series):
        x, y, leg, est = (tuple(s) + (None, None))[:4] if len(s) < 4 else s
        nome = leg or "série %d" % i
        x, y = _limpa_serie(list(x), list(y), nome, logx, logy)
        if not x:
            _aviso("série '%s' vazia: ignorada" % nome)
            continue
        e = estilo_padrao(len(limpas))
        if est:
            e.update(est)
        limpas.append((x, y, leg, e))
    if not limpas:
        raise ValueError("nenhuma série com dados válidos para %s" % path)

    xs = [v for s in limpas for v in s[0]]
    ys = [v for s in limpas for v in s[1]]
    xs_w = xs + [v[0] for v in vlines if math.isfinite(v[0])
                 and not (logx and v[0] <= 0)]
    x0, y0, x1, y1, xt, yt = nice_world(xs_w, ys, logx, logy)
    if world is not None:
        wx0, wy0, wx1, wy1 = world
        if wx0 is not None: x0 = wx0
        if wy0 is not None: y0 = wy0
        if wx1 is not None: x1 = wx1
        if wy1 is not None: y1 = wy1
        if not logx:
            xt = nice_tick(x1 - x0)
        if not logy:
            yt = nice_tick(y1 - y0)

    # retas verticais: 2 pontos cobrindo o world em y
    todas = list(limpas)
    for (xv, est, leg) in vlines:
        if not math.isfinite(xv) or (logx and xv <= 0):
            _aviso("reta vertical em x=%r descartada" % (xv,))
            continue
        e = dict(color=COR_CINZA, linestyle=3, linewidth=1.5)
        if est:
            e.update(est)
        todas.append(([xv, xv], [y0, y1], leg, e))

    L = _cabecalho_fixo()
    L += ["@with g0",
          "@    world %s, %s, %s, %s" % (_fmt(x0), _fmt(y0), _fmt(x1), _fmt(y1)),
          "@    stack world 0, 0, 0, 0", "@    znorm 1",
          "@    view 0.180000, 0.150000, 1.150000, 0.850000",
          "@    title " + _q(title or ""), "@    title font 0",
          "@    title size 1.500000", "@    title color 1",
          '@    subtitle ""', "@    subtitle font 0",
          "@    subtitle size 1.000000", "@    subtitle color 1",
          "@    xaxes scale " + ("Logarithmic" if logx else "Normal"),
          "@    yaxes scale " + ("Logarithmic" if logy else "Normal"),
          "@    xaxes invert off", "@    yaxes invert off"]
    L += _eixo("x", xlabel, xt, 8 if logx else 4, "0.000000, 0.080000", False, logx,
                 _prec_sig(x0, x1, xt))
    L += _eixo("y", ylabel, yt, 8 if logy else 4, "0.000000, 0.110000", True, logy,
                 _prec_sig(y0, y1, yt))
    L += ["@    altxaxis  off", "@    altyaxis  off", "@    legend on",
          "@    legend loctype view",
          "@    legend %s, %s" % (_fmt(legend_pos[0]), _fmt(legend_pos[1])),
          "@    legend box color 1", "@    legend box pattern 1",
          "@    legend box linewidth 1.0", "@    legend box linestyle 0",
          "@    legend box fill color 0", "@    legend box fill pattern 0",
          "@    legend font 0", "@    legend char size 1.200000",
          "@    legend color 1", "@    legend length 4", "@    legend vgap 1",
          "@    legend hgap 1", "@    legend invert false",
          "@    frame type 0", "@    frame linestyle 1",
          "@    frame linewidth 1.5", "@    frame color 1",
          "@    frame pattern 1", "@    frame background color 0",
          "@    frame background pattern 0"]
    for i, (x, y, leg, e) in enumerate(todas):
        L += _bloco_serie(i, e, leg)
    for i, (x, y, leg, e) in enumerate(todas):
        L += ["@target G0.S%d" % i, "@type xy"]
        L += ["%s %s" % (_fmt(a), _fmt(b)) for a, b in zip(x, y)]
        L.append("&")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
    return path


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Gera um .agr (Xmgrace) a partir de arquivos de dados de 2 colunas.")
    ap.add_argument("saida")
    ap.add_argument("dados", nargs="+")
    ap.add_argument("--xlabel", default="x")
    ap.add_argument("--ylabel", default="y")
    ap.add_argument("--legend", action="append", default=[])
    ap.add_argument("--title", default=None)
    ap.add_argument("--logx", action="store_true")
    ap.add_argument("--logy", action="store_true")
    ap.add_argument("--cols", default="0,1", help="colunas x,y (base 0)")
    a = ap.parse_args(argv)
    try:
        cols = tuple(int(c) for c in a.cols.split(","))
        series = []
        for i, d in enumerate(a.dados):
            x, y = read_columns(d, cols)
            leg = a.legend[i] if i < len(a.legend) else ""
            series.append((x, y, leg, None))
        write_agr(a.saida, series, a.xlabel, a.ylabel, title=a.title,
                  logx=a.logx, logy=a.logy)
    except (OSError, ValueError) as e:
        print("[agr] ERRO: %s" % e, file=sys.stderr)
        return 1
    print(a.saida)
    return 0


if __name__ == "__main__":
    sys.exit(main())
