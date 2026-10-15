#!/usr/bin/env bash
# Verbatim gate for the grapes-engine restructure (Task 1, decisions D1-D4).
#
# main 614be3c7 had a flat src/. This branch moves that code into src/abcmk/, src/data/,
# src/utils.jl and src/alpha_b.jl by whole functions, without changing a line. The gate
# proves it without Julia: the multiset of non-blank lines of main's src/*.jl (except
# MKtest.jl) must equal the multiset of non-blank lines of every branch src/**/*.jl outside
# src/ml/ (and except MKtest.jl), up to the lines listed in test/verbatim_allowlist.txt:
#   "- <line>"  a line of main allowed to disappear
#   "+ <line>"  a line allowed to appear
# ext/ must be byte-identical to main, and main's test/runtests.jl must be a prefix of the
# branch's (new testsets are only appended).
#
# Usage: bash test/verbatim_gate.sh   (from the repository root; needs git)
set -euo pipefail

BASE="614be3c7"
ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"
ALLOW="test/verbatim_allowlist.txt"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
export LC_ALL=C

nonblank() { grep -v '^[[:space:]]*$' || true; }

# main side
: > "$TMP/main.txt"
while IFS= read -r f; do
    [[ "$f" == "src/MKtest.jl" ]] && continue
    git show "$BASE:$f" | nonblank >> "$TMP/main.txt"
done < <(git ls-tree -r --name-only "$BASE" -- src | grep '\.jl$')

# branch side
: > "$TMP/branch.txt"
while IFS= read -r f; do
    nonblank < "$f" >> "$TMP/branch.txt"
done < <(find src -name '*.jl' ! -path 'src/ml/*' ! -path 'src/MKtest.jl' | sort)

sort "$TMP/main.txt" > "$TMP/main.sorted"
sort "$TMP/branch.txt" > "$TMP/branch.sorted"
comm -23 "$TMP/main.sorted" "$TMP/branch.sorted" > "$TMP/removed"
comm -13 "$TMP/main.sorted" "$TMP/branch.sorted" > "$TMP/added"

: > "$TMP/allow_removed"
: > "$TMP/allow_added"
if [[ -f "$ALLOW" ]]; then
    { grep '^- ' "$ALLOW" || true; } | sed 's/^- //' | nonblank | sort > "$TMP/allow_removed"
    { grep '^+ ' "$ALLOW" || true; } | sed 's/^+ //' | nonblank | sort > "$TMP/allow_added"
fi

status=0
report() {
    local label="$1" got="$2" want="$3"
    echo "== $label: $(wc -l < "$got") line(s) (allowed: $(wc -l < "$want"))"
    sed 's/^/    /' "$got"
    if ! cmp -s "$got" "$want"; then
        echo "   FAIL: $label differs from the allow-list"
        diff "$want" "$got" | sed 's/^/    /' || true
        status=1
    fi
}
report "lines of main removed or changed" "$TMP/removed" "$TMP/allow_removed"
report "lines added" "$TMP/added" "$TMP/allow_added"

echo "== ext/ byte-identical to $BASE"
while IFS= read -r f; do
    if ! git show "$BASE:$f" | cmp -s - "$f"; then
        echo "   FAIL: $f differs"; status=1
    fi
done < <(git ls-tree -r --name-only "$BASE" -- ext)
extra_ext="$(comm -13 <(git ls-tree -r --name-only "$BASE" -- ext | sort) <(find ext -type f | sort))"
if [[ -n "$extra_ext" ]]; then echo "   FAIL: new files in ext/: $extra_ext"; status=1; fi

echo "== main's test/runtests.jl is a prefix of the branch's"
main_rt="$TMP/runtests.main"
git show "$BASE:test/runtests.jl" > "$main_rt"
if ! head -c "$(wc -c < "$main_rt")" test/runtests.jl | cmp -s - "$main_rt"; then
    echo "   FAIL: test/runtests.jl changes main's content"; status=1
fi

if [[ "$status" -eq 0 ]]; then echo "GATE PASS"; else echo "GATE FAIL"; fi
exit "$status"
