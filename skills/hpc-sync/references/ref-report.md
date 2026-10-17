# The ref report

One screen that proves a sync, pasted into the run notes. Produce it after the push and
after the clone's fast-forward, never before.

In the clone you pushed from, with the cluster remote's name in place of
`<cluster-remote>` and the synced branches in place of `<branches>`:

```bash
git fetch -q origin && git fetch -q <cluster-remote> && {
  echo "== branches (after push)"
  for b in <branches>; do
    o="$(git rev-parse "origin/$b")"; c="$(git rev-parse "<cluster-remote>/$b" 2>/dev/null || echo none)"
    [[ "$o" == "$c" ]] && v=equal || v=DIFFERENT
    printf '%-6s %s %s %s\n' "$b" "${o:0:12}" "${c:0:12}" "$v"
  done
  echo "== cluster-only refs still present"
  git ls-remote --heads --tags <cluster-remote> | awk '{print substr($1,1,12), $2}'
}
```

On the cluster, append the clone's state:

```bash
source "$HVANTK_SITE_PROFILE" && echo "== working clone" && git -C "$HVANTK_REPO" log --oneline -1 && git -C "$HVANTK_REPO" status --porcelain | wc -l
```

A report where any branch reads `DIFFERENT`, or where the clone's commit is not the
branch's hash, is not a finished sync.

## Reading the fast-forward check

| result of `git merge-base --is-ancestor <cluster>/<b> origin/<b>` | meaning | action |
| --- | --- | --- |
| exit 0 | the cluster's commit is an ancestor of GitHub's: pushing is a fast-forward | push by refspec |
| exit 1 | the cluster has commits GitHub does not | stop; show `git log origin/<b>..<cluster>/<b>` to the user |
| error: not a valid object | the cluster has no such branch yet | the push creates it; still verify afterwards |
