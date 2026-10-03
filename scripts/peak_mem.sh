#!/usr/bin/env bash
# Run a command and record the peak resident memory of it and all its children.
# Usage: scripts/peak_mem.sh <log_file> <command...>   (appends PEAK_RSS_GB and WALL_SECONDS)
log=$1; shift
start=$(date +%s); "$@" & pid=$!; peak=0
while kill -0 $pid 2>/dev/null; do
  rss=$(ps -o rss= --ppid $pid -p $pid 2>/dev/null | awk '{s+=$1} END {print s+0}')
  kids=$(pgrep -P $pid | xargs -r -I{} pgrep -P {} | xargs -r ps -o rss= -p 2>/dev/null | awk '{s+=$1} END {print s+0}')
  tot=$((rss + kids)); (( tot > peak )) && peak=$tot; sleep 10
done
wait $pid; code=$?
echo "PEAK_RSS_GB $(awk -v k=$peak 'BEGIN {printf "%.1f", k/1048576}')" >> "$log"
echo "WALL_SECONDS $(( $(date +%s) - start ))" >> "$log"
exit $code
