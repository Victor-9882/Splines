#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""plot_run.py -- gera os gráficos .agr (e opcionalmente PNG) de um run.

Uso:  python3 plot_run.py <run_dir> [--png]

Lê <run_dir>/derived/*.dat e <run_dir>/derived/observaveis.txt e escreve em
<run_dir>/plots/agr/: u_xi, g_gamma, g_z, psi_gamma, psi_z e dperp (se existir).
Também exporta as definições das grandezas, usadas por plot_compare.py.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import agr  # noqa: E402

# nome -> (arquivo em derived/, rótulo x, rótulo y, nome do .agr)
G = r"\f{Symbol}g"
XI = r"\f{Symbol}x"
ZZ = r"\f{Times-Italic}z"
QUANTIDADES = {
    "u":         ("dist_long.dat",     XI, r"\f{Times-Italic}u", "u_xi"),
    "g_gamma":   ("autovetoresG.dat",  G,  r"\f{Times-Italic}g", "g_gamma"),
    "g_z":       ("autovetoresZ.dat",  ZZ, r"\f{Times-Italic}g", "g_z"),
    "psi_gamma": ("psi.dat",           G,  r"\f{Symbol}y", "psi_gamma"),
    "psi_z":     ("psi_z.dat",         ZZ, r"\f{Symbol}y", "psi_z"),
    "dperp":     ("dist_trans.dat",    G,
                  r"\f{Times-Italic}D\s\f{Symbol}^\N\f{}(\f{Symbol}g\f{})", "dperp"),
}
# legenda de <xi>: "<" e ">" em Symbol (angleleft/angleright)
LEG_XI_MEDIO = r"\f{Symbol}\#{e1}x\#{f1}\f{} = %.4f"


def ler_chave_valor(path):
    """Lê arquivo 'chave = valor' (comentários # ou ! ignorados).
    Retorna dict chave(minúscula) -> float ou str.  {} se não existir."""
    d = {}
    if not os.path.isfile(path):
        return d
    with open(path, errors="replace") as f:
        for linha in f:
            s = linha.strip()
            if not s or s[0] in "#!" or "=" not in s:
                continue
            k, v = s.split("=", 1)
            v = v.strip().strip("'\"")
            try:
                d[k.strip().lower()] = agr._para_float(v)
            except ValueError:
                d[k.strip().lower()] = v
    return d


def ler_params_nml(run_dir):
    """Leitura leve de input/params.nml (apenas pares chave = valor)."""
    p = os.path.join(run_dir, "input", "params.nml")
    d = {}
    if not os.path.isfile(p):
        return d
    with open(p, errors="replace") as f:
        txt = f.read()
    for k, v in re.findall(r"(\w+)\s*=\s*([^,/\s]+)", txt):
        v = v.strip("'\"")
        try:
            d[k.lower()] = agr._para_float(v)
        except ValueError:
            d[k.lower()] = v
    return d


def observaveis(run_dir):
    """observaveis.txt com fallback para params.nml (chaves ausentes)."""
    d = ler_params_nml(run_dir)
    d.update(ler_chave_valor(os.path.join(run_dir, "derived", "observaveis.txt")))
    return d


def carrega_quantidade(run_dir, q):
    """Retorna (x, y) da grandeza q do run, ou None (com aviso) se não houver."""
    arq = os.path.join(run_dir, "derived", QUANTIDADES[q][0])
    if not os.path.isfile(arq):
        agr._aviso("%s não encontrado: %s pulado" % (arq, q))
        return None
    x, y = agr.read_columns(arq, (0, 1))
    if not x:
        agr._aviso("%s sem dados numéricos: %s pulado" % (arq, q))
        return None
    return x, y


def para_png(agr_path, png_path):
    """Converte .agr em PNG via gracebat.  Retorna (ok, mensagem)."""
    exe = shutil.which("gracebat")
    if not exe:
        return False, "gracebat não encontrado no PATH"
    png_path = os.path.abspath(png_path)
    pasta = os.path.dirname(png_path)
    os.makedirs(pasta, exist_ok=True)
    # o gracebat aborta ("buffer overflow") com caminhos de ~256 bytes ou mais:
    # roda dentro da pasta do PNG, com caminhos relativos curtos
    r = subprocess.run([exe, "-hardcopy", "-hdevice", "PNG", "-printfile",
                        os.path.basename(png_path), os.path.relpath(os.path.abspath(agr_path), pasta)],
                       capture_output=True, text=True, cwd=pasta)
    saida = (r.stdout + r.stderr).strip()
    ok = r.returncode == 0 and os.path.isfile(png_path)
    if saida:
        agr._aviso("gracebat (%s): %s" % (os.path.basename(agr_path), saida))
    return ok, saida


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("--png", action="store_true", help="gera também PNG via gracebat")
    a = ap.parse_args(argv)
    run = a.run_dir.rstrip("/")
    if not os.path.isdir(run):
        print("[plot_run] ERRO: run_dir inexistente: %s" % run, file=sys.stderr)
        return 1
    if not os.path.isdir(os.path.join(run, "derived")):
        print("[plot_run] ERRO: %s não tem a pasta derived/ (rodar o programa/pós-processamento antes)"
              % run, file=sys.stderr)
        return 1
    dir_agr = os.path.join(run, "plots", "agr")
    dir_png = os.path.join(run, "plots", "png")
    os.makedirs(dir_agr, exist_ok=True)
    obs = observaveis(run)
    gerados = []
    for q, (_, xl, yl, nome) in QUANTIDADES.items():
        dados = carrega_quantidade(run, q)
        if dados is None:
            continue
        x, y = dados
        est = dict(color=agr.COR_TIJOLO, linestyle=1, linewidth=2.0)
        vl = []
        if q == "u" and isinstance(obs.get("xi_medio"), float):
            vl.append((obs["xi_medio"],
                       dict(color=agr.COR_AZUL, linestyle=3, linewidth=1.5),
                       LEG_XI_MEDIO % obs["xi_medio"]))
        out = os.path.join(dir_agr, nome + ".agr")
        try:
            agr.write_agr(out, [(x, y, "", est)], xl, yl, vlines=vl)
        except (ValueError, OSError) as e:
            agr._aviso("%s: %s" % (nome, e))
            continue
        gerados.append(out)
        if a.png:
            png = os.path.join(dir_png, nome + ".png")
            ok, msg = para_png(out, png)
            if ok:
                gerados.append(png)
            else:
                agr._aviso("PNG de %s não gerado: %s" % (nome, msg))
    if not gerados:
        print("[plot_run] ERRO: nada foi gerado para %s" % run, file=sys.stderr)
        return 1
    print("\n".join(gerados))
    return 0


if __name__ == "__main__":
    sys.exit(main())
