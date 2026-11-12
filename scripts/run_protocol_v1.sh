#!/usr/bin/env bash
# Protocol v1 (docs/protocol_matched_v1.md), one part per call. CPU only. Run each call in its
# own byobu session, and check the machine's load first (other users' GPU jobs come first).
#   scripts/run_protocol_v1.sh small <seed> [<seed> ...]    the nine small datasets, about 2 h per seed
#   scripts/run_protocol_v1.sh pbmc <seed> <stage> [b]      stage = full | mica | pca | score
#   scripts/run_protocol_v1.sh chooser <seed> [<seed> ...]  chooseR's R code on the nine small datasets
#   scripts/run_protocol_v1.sh chooser-pbmc <seed>          chooseR's R code on PBMC14K
# Confirmatory runs pass --require-tag: they refuse to start unless the code they execute is the
# code tagged matched-protocol-v1. Set PRETAG=1 to run the pre-tag checks without it.
# A small dataset whose provenance file shows a finished, tag-verified run at B=100 is skipped,
# so a stopped loop can be restarted with the same command. A failure prints a FAILED line.
set -u
cd "$(dirname "$0")/.."
source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
# The package is installed in editable mode from the main checkout. Put this checkout's src first,
# so a run from another checkout uses its own code (each provenance file records which ran).
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
guard=--require-tag
[ "${PRETAG:-0}" = 1 ] && guard=
part=$1
shift
SMALL="paper_yan paper_goolam paper_buettner paper_pollen paper_chung paper_usoskin paper_kolodziejczyk paper_klein paper_zeisel"
finished() {  # finished <provenance json>: a tag-verified run at B=100 wrote it
  [ -f "$1" ] && python -c "import json,sys; d=json.load(open(sys.argv[1])); sys.exit(0 if d['git'].get('code_matches_tag') and d.get('b')==100 else 1)" "$1"
}
run() {  # run <label> <command...>: stream the filtered log, report a failure
  "${@:2}" 2>&1 | grep --line-buffered -vE "warn|Warn|dispatch|return "
  local st=${PIPESTATUS[0]}
  [ "$st" = 0 ] || echo "FAILED $1 (exit $st) $(date)"
  return "$st"
}
case $part in
  small)
    for s in "$@"; do
      ok=1
      for d in $SMALL; do
        if [ -n "$guard" ] && finished "results/v1/provenance/${d}_seed${s}.json"; then
          echo "######## v1 small $d seed $s already finished, skipped"; continue
        fi
        echo "######## v1 small $d seed $s $(date)"
        run "small $d seed $s" python -u scripts/run_protocol_small.py "$d" "$s" 100 $guard || ok=0
      done
      [ $ok = 1 ] && echo "V1 SMALL SEED $s DONE $(date)" || echo "V1 SMALL SEED $s INCOMPLETE $(date)"
    done ;;
  pbmc)
    s=$1 stage=$2 b=${3:-100}
    echo "######## v1 pbmc14k seed $s stage $stage $(date)"
    run "pbmc seed $s $stage" python -u scripts/run_protocol_pbmc.py "$s" "$stage" "$b" $guard ${FRESH:+--fresh} \
      && echo "V1 PBMC SEED $s STAGE $stage DONE $(date)" ;;
  chooser)
    for s in "$@"; do
      for d in $SMALL; do
        cc=0
        case $d in paper_klein|paper_zeisel) cc=5 ;; esac
        echo "######## v1 chooser $d seed $s $(date)"
        run "chooser $d seed $s" python -u scripts/run_chooser_r.py "$d" "$s" --workers 8 --crosscheck $cc $guard
      done
    done ;;
  chooser-pbmc)
    s=$1
    echo "######## v1 chooser paper_pbmc14k seed $s $(date)"
    run "chooser pbmc seed $s" python -u scripts/run_chooser_r.py paper_pbmc14k "$s" --workers 4 --crosscheck 2 $guard ;;
  *) echo "unknown part $part"; exit 1 ;;
esac
