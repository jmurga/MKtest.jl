#!/usr/bin/env bash
# Run what .github/workflows/CI.yml runs, on a clean export of a commit, before pushing:
#   1. tests:  remove Manifest.toml, resolve, build, Pkg.test() with JULIA_NUM_THREADS=1
#   2. docs:   docs project with the package developed, docs/make.jl, then doctest
# The commit is checked out in a temporary detached worktree (tracked files only, plus the
# git metadata Documenter needs for its source links), so an untracked docs/Manifest.toml
# or a local edit cannot make the check pass where CI fails.
#   scripts/check_ci.sh [rev]        (default HEAD; logs kept in the printed folder)
set -euo pipefail

rev="${1:-HEAD}"
repo="$(git rev-parse --show-toplevel)"
sha="$(git -C "$repo" rev-parse --short "$rev")"
work="$(mktemp -d "${TMPDIR:-/tmp}/mktest_ci_${sha}_XXXX")"
src="$work/MKtest.jl"
git -C "$repo" worktree add --quiet --detach "$src" "$rev"
trap 'git -C "$repo" worktree remove --force "$src"' EXIT
echo "check_ci: $sha checked out in $src (logs in $work)"

cd "$src"
rm -f Manifest.toml

echo "check_ci: tests ..."
if ! JULIA_NUM_THREADS=1 julia --project=. -e '
        using Pkg
        Pkg.Registry.update()
        Pkg.instantiate()
        Pkg.build()
        Pkg.test()' > "$work/tests.log" 2>&1; then
    tail -40 "$work/tests.log"
    echo "check_ci: FAILED (tests), $sha; full log $work/tests.log"
    exit 1
fi
grep -E "Test Summary|tests passed" "$work/tests.log" | tail -8

echo "check_ci: docs ..."
if ! GKSwstype=100 julia --project=docs -e '
        using Pkg
        Pkg.develop(PackageSpec(path = pwd()))
        Pkg.instantiate()
        include("docs/make.jl")
        using Documenter: DocMeta, doctest
        DocMeta.setdocmeta!(MKtest, :DocTestSetup, :(using MKtest); recursive = true)
        doctest(MKtest)' > "$work/docs.log" 2>&1; then
    grep -E "Error|ERROR" "$work/docs.log" | head -20
    echo "check_ci: FAILED (docs), $sha; full log $work/docs.log"
    exit 1
fi

echo "check_ci: PASSED, $sha (tests and docs as in CI.yml)"
