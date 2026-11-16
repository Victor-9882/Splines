#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""plot_compare.py -- comparações entre runs (gráficos Xmgrace .agr).

Subcomandos:
  overlay --runs R1 R2 ... --quantity {u,g_gamma,g_z,psi_gamma,psi_z,dperp}
          --legend-key {razao,zeta,run_id,<chave de observaveis.txt>}
          --out arq.agr [--png] [--vlines]
  conv    --csv resumo.csv --param NMZ --y {lambda,alfa,xi_medio,norma_u}
          --out arq.agr [--filter col=valor ...] [--logerr] [--ref V] [--png]
"""
import argparse
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import agr  # noqa: E402
import plot_run as pr  # noqa: E402

Y_ROTULOS = {
    "lambda":   r"\f{Symbol}l",
    "alfa":     r"\f{Symbol}a",
    "xi_medio": r"\f{Symbol}\#{e1}x\#{f1}\f{}",
    "norma_u":  r"\f{Times-Italic}N\su\N\f{}",
}


def _g(v):
    return "%g" % v


def legenda_run(run_dir, chave, obs):
    """Legenda de um run conforme --legend-key."""
    rid = os.path.basename(os.path.abspath(run_dir))
    if chave == "run_id":
        return rid
    if chave == "razao":
        m1, m2 = obs.get("m1"), obs.get("m2")
        if isinstance(m1, float) and isinstance(m2, float) and m1 != 0:
            return r"m\s2\N/m\s1\N = " + _g(round(m2 / m1, 6))
        agr._aviso("%s: m1/m2 ausentes; usando run_id na legenda" % rid)
        return rid
    if chave == "zeta":
        if isinstance(obs.get("zeta"), float):
            return r"\f{Symbol}z\f{} = " + _g(obs["zeta"])
        agr._aviso("%s: zeta ausente; usando run_id na legenda" % rid)
        return rid
    v = obs.get(chave.lower())
    if v is None:
        agr._aviso("%s: chave '%s' ausente em observaveis.txt; usando run_id" % (rid, chave))
        return rid
    return "%s = %s" % (chave, _g(v) if isinstance(v, float) else v)


def _salva(out, series, xl, yl, png, **kw):
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    try:
        agr.write_agr(out, series, xl, yl, **kw)
    except (ValueError, OSError) as e:
        print("[plot_compare] ERRO: %s" % e, file=sys.stderr)
        return 1
    print(out)
    if png:
        p = os.path.splitext(out)[0] + ".png"
        ok, msg = pr.para_png(out, p)
        if ok:
            print(p)
        else:
            print("[plot_compare] ERRO: PNG não gerado: %s" % msg, file=sys.stderr)
            return 1
    return 0


def cmd_overlay(a):
    xl, yl = pr.QUANTIDADES[a.quantity][1:3]
    series, vlines = [], []
    for run in a.runs:
        run = run.rstrip("/")
        dados = pr.carrega_quantidade(run, a.quantity)
        if dados is None:
            continue
        obs = pr.observaveis(run)
        i = len(series)
        series.append((dados[0], dados[1], legenda_run(run, a.legend_key, obs), None))
        if a.vlines:
            if a.quantity != "u":
                agr._aviso("--vlines só vale para a grandeza u; ignorado")
            elif isinstance(obs.get("xi_medio"), float):
                est = agr.estilo_padrao(i)
                vlines.append((obs["xi_medio"],
                               dict(color=est["color"], linestyle=2, linewidth=2.0), ""))
            else:
                agr._aviso("%s: xi_medio ausente; sem reta vertical" % run)
    if not series:
        print("[plot_compare] ERRO: nenhum run com dados para '%s'" % a.quantity,
              file=sys.stderr)
        return 1
    return _salva(a.out, series, xl, yl, a.png, vlines=vlines)


def _num(s):
    try:
        return agr._para_float(s)
    except (ValueError, AttributeError):
        return None


def cmd_conv(a):
    try:
        with open(a.csv, newline="", errors="replace") as f:
            linhas = list(csv.DictReader(f))
    except OSError as e:
        print("[plot_compare] ERRO: %s" % e, file=sys.stderr)
        return 1
    if not linhas:
        print("[plot_compare] ERRO: CSV vazio: %s" % a.csv, file=sys.stderr)
        return 1
    cols = list(linhas[0].keys())
    for c in [a.param, a.y] + [f.split("=")[0] for f in a.filter if "=" in f]:
        if c not in cols:
            print("[plot_compare] ERRO: coluna '%s' ausente no CSV (colunas: %s)"
                  % (c, ", ".join(cols)), file=sys.stderr)
            return 1
    for f in a.filter:
        if "=" not in f:
            print("[plot_compare] ERRO: filtro inválido '%s' (use col=valor)" % f,
                  file=sys.stderr)
            return 1
        c, v = f.split("=", 1)
        vn = _num(v)
        def ok(l, c=c, v=v, vn=vn):
            ln = _num(l.get(c))
            return (ln is not None and vn is not None and abs(ln - vn) <= 1e-9 * max(1, abs(vn))) \
                or (l.get(c) or "").strip() == v
        linhas = [l for l in linhas if ok(l)]
    pts = []
    for l in linhas:
        x, y = _num(l.get(a.param)), _num(l.get(a.y))
        if x is None or y is None:
            agr._aviso("linha run_id=%s com %s/%s vazio ou inválido: ignorada"
                       % (l.get("run_id", "?"), a.param, a.y))
            continue
        pts.append((x, y))
    if not pts:
        print("[plot_compare] ERRO: nenhuma linha restante após filtros", file=sys.stderr)
        return 1
    pts.sort()
    xs = [p[0] for p in pts]
    if len(set(xs)) < len(xs):
        agr._aviso("valores repetidos de '%s' após filtros; refine com --filter" % a.param)
    ys = [p[1] for p in pts]
    xl = a.param
    yl = Y_ROTULOS[a.y]
    estilo = dict(color=agr.COR_TIJOLO, linestyle=1, linewidth=1.5,
                  symbol="circle", symbol_size=0.8)
    kw = {}
    if a.logerr:
        if a.ref is not None:
            ref, usados = a.ref, list(zip(xs, ys))
        else:
            if len(pts) < 2:
                print("[plot_compare] ERRO: --logerr precisa de >= 2 pontos", file=sys.stderr)
                return 1
            ref, usados = ys[-1], list(zip(xs[:-1], ys[:-1]))
        xs = [p[0] for p in usados]
        ys = [abs(p[1] - ref) for p in usados]
        yl = r"|\f{Symbol}D\f{}%s|" % Y_ROTULOS[a.y]  # |y_N - y_ref|
        kw["logy"] = True
        leg = "ref = " + _g(ref)
    else:
        leg = ""
    return _salva(a.out, [(xs, ys, leg, estilo)], xl, yl, a.png, **kw)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    o = sub.add_parser("overlay")
    o.add_argument("--runs", nargs="+", required=True)
    o.add_argument("--quantity", required=True, choices=list(pr.QUANTIDADES))
    o.add_argument("--legend-key", default="run_id")
    o.add_argument("--out", required=True)
    o.add_argument("--png", action="store_true")
    o.add_argument("--vlines", action="store_true")
    c = sub.add_parser("conv")
    c.add_argument("--csv", required=True)
    c.add_argument("--param", required=True)
    c.add_argument("--y", required=True, choices=list(Y_ROTULOS))
    c.add_argument("--out", required=True)
    c.add_argument("--filter", nargs="*", action="extend", default=[])
    c.add_argument("--logerr", action="store_true")
    c.add_argument("--ref", type=float, default=None)
    c.add_argument("--png", action="store_true")
    a = ap.parse_args(argv)
    return cmd_overlay(a) if a.cmd == "overlay" else cmd_conv(a)


if __name__ == "__main__":
    sys.exit(main())
