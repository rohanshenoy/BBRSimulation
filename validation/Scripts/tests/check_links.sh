#!/bin/bash
# check_links.sh REPO [FILE...] — every relative markdown link in FILE (default:
# every tracked *.md) names a tracked file (staged counts), and every in-page
# #anchor matches the GitHub slug of a heading in the same file. Prints one
# "ok|url|BAD" line per link, then LINKS OK or LINKS BAD; exit 0 or 1.
# Not checked: anchors into other files (file.md#x), unicode in slugs.
cd "$1" || exit 1; shift; bad=0
if [ $# -eq 0 ]; then
  files=(); while IFS= read -r f; do files+=("$f"); done < <(git ls-files '*.md')
  [ ${#files[@]} -gt 0 ] || { echo "BAD  no tracked *.md files in $PWD"; echo "LINKS BAD"; exit 1; }
  set -- "${files[@]}"
fi
# -H: with a single FILE grep would omit the file name, and the link would be
# resolved against REPO instead of the file's directory.
while IFS= read -r l; do
  f="${l%%:*}"; p="${l#*](}"
  case "$p" in http://*|https://*) echo "url  $p"; continue;; esac
  t="$(dirname "$f")/$p"
  if [ -e "$t" ] && [ -n "$(git ls-files -- "$t")" ]; then echo "ok   $f -> $p"
  else echo "BAD  $f -> $p ($t)"; bad=1; fi
done < <(grep -oHE '\]\(([^)#]+)' "$@")
# In-page anchors: GitHub slug of every heading.
for f in "$@"; do
  for a in $(grep -oE '\]\(#[^)]+\)' "$f" | sed -E 's/^\]\(#//; s/\)$//'); do
    if grep -E '^#+ ' "$f" | sed -E 's/^#+ //' | tr 'A-Z' 'a-z' | sed -E 's/[^a-z0-9 _-]//g; s/ /-/g' | grep -qx "$a"; then
      echo "ok   $f #$a"; else echo "BAD  $f #$a"; bad=1; fi
  done
done
[ $bad -eq 0 ] && echo "LINKS OK" || { echo "LINKS BAD"; exit 1; }
