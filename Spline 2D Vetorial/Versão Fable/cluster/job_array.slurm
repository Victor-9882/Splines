#!/bin/bash
#=======================================================================
# Array de jobs: uma linha de cluster/manifest.csv por tarefa.
#
#   NÃO SUBMETIDO. Antes de usar, ajuste as linhas marcadas com AJUSTAR
#   (partição, tempo, módulos) para o cluster escolhido e confira o
#   intervalo do --array com o número de linhas do manifest.
#
#   sbatch --array=1-<N>%<simultâneos> cluster/job_array.slurm
#
# Cada tarefa escreve só em Testes/<teste>/<run_id>/ (run_id único por
# linha do manifest) e roda com --no-csv: não há arquivo compartilhado
# entre jobs. Depois: python3 cluster/collect.py Testes/<teste> --plots
#=======================================================================
#SBATCH --job-name=splines_bs
#SBATCH --array=1-2                    # AJUSTAR: 1-<linhas do manifest>[%<máx. simultâneos>]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16             # AJUSTAR conforme os nós
#SBATCH --mem=4G                       # matrizes são pequenas (ver cluster/plano.md)
#SBATCH --time=12:00:00                # AJUSTAR: maior t_est_h do manifest com folga
#SBATCH --partition=AJUSTAR            # AJUSTAR: nome da partição/fila
#SBATCH --output=cluster/logs/%x_%A_%a.out
#SBATCH --error=cluster/logs/%x_%A_%a.err

set -euo pipefail

# Raiz do projeto (pasta que contém src/, tools/, cluster/). Submeta a partir dela.
RAIZ="${RAIZ:-${SLURM_SUBMIT_DIR:-$PWD}}"
MANIFEST="${MANIFEST:-$RAIZ/cluster/manifest.csv}"
BIN="${BIN:-$RAIZ/bin}"
TAREFA="${SLURM_ARRAY_TASK_ID:?defina SLURM_ARRAY_TASK_ID (índice da linha do manifest)}"
NTHREADS="${SLURM_CPUS_PER_TASK:-4}"

# ---------------- Ambiente (AJUSTAR para o cluster) -------------------
# Opção 1: módulos
#   module purge
#   module load gcc/13 intel-oneapi-mkl python/3
# Opção 2: oneAPI instalado em /opt (como na máquina local)
if [ -z "${MKLROOT:-}" ] && [ -f /opt/intel/oneapi/setvars.sh ]; then
    set +u; source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1 || true; set -u
fi

export OMP_NUM_THREADS="$NTHREADS"
export MKL_NUM_THREADS="$NTHREADS"
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=512M

cd "$RAIZ"
mkdir -p cluster/logs

# Linha do manifest: idx,teste,run_id,parametros,...   (parametros = "CHAVE=VALOR CHAVE=VALOR ...")
LINHA=$(awk -F, -v n="$TAREFA" 'NR > 1 && $1 == n' "$MANIFEST")
if [ -z "$LINHA" ]; then
    echo "ERRO: índice $TAREFA não existe em $MANIFEST" >&2
    exit 1
fi
TESTE=$(echo "$LINHA" | cut -d, -f2)
RUN_ID=$(echo "$LINHA" | cut -d, -f3)
PARAMS=$(echo "$LINHA" | cut -d, -f4)

echo "tarefa $TAREFA | teste $TESTE | run_id $RUN_ID | threads $NTHREADS | host $(hostname) | $(date '+%F %T')"

# Confere que o run_id do manifest é o que o run.py vai gerar (evita colisão/erro de digitação)
# shellcheck disable=SC2086
GERADO=$(python3 tools/run.py --teste "$TESTE" --print-run-id --set $PARAMS)
if [ "$GERADO" != "$RUN_ID" ]; then
    echo "ERRO: run_id do manifest ($RUN_ID) difere do gerado ($GERADO)" >&2
    exit 1
fi

# shellcheck disable=SC2086
python3 tools/run.py --teste "$TESTE" --bin "$BIN" --threads "$NTHREADS" \
        --no-csv --skip-existing --set $PARAMS

echo "fim $(date '+%F %T')"
