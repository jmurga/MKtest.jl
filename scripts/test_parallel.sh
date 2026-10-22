#!/usr/bin/env bash
# Run the test groups of test/runtests.jl as parallel single-thread processes, one
# `Pkg.test(test_args = [group])` each. Same tests as `Pkg.test()`, about a third of the
# wall time: most of the suite is compilation, which each process pays only for what it
# uses. The frozen fixtures need one thread, which every process gets.
#   scripts/test_parallel.sh [group ...]      (default: every group; logs in the printed folder)
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
groups=("$@")
[[ ${#groups[@]} -eq 0 ]] && groups=(rates abc mk engine extension)
logdir="$(mktemp -d "${TMPDIR:-/tmp}/mktest_tests_XXXX")"
echo "test_parallel: ${groups[*]} in $repo (logs in $logdir)"

pids=()
for g in "${groups[@]}"; do
    (cd "$repo" && JULIA_NUM_THREADS=1 julia --project=. -e "using Pkg; Pkg.test(test_args = [\"$g\"])" \
        > "$logdir/$g.log" 2>&1) &
    pids+=("$!")
done

fail=0
for i in "${!pids[@]}"; do
    g="${groups[$i]}"
    if wait "${pids[$i]}"; then
        echo "test_parallel: PASSED $g  $(grep -E '^[A-Za-z].*\|' "$logdir/$g.log" | grep -vE 'Test Summary|Pass +Total' | head -1 | sed 's/  */ /g')"
    else
        echo "test_parallel: FAILED $g (log $logdir/$g.log)"
        grep -E "Test Failed|Error During Test|ERROR" "$logdir/$g.log" | head -10
        fail=1
    fi
done
exit "$fail"
