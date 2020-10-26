
Build documentation from in here (`doc`) by running
```bash
pdoc --html --force --template-dir ./templates ../vorts
```

This generates the site in `./html/vorts`.

Push to `gh-pages` branch on GitHub using [`gph-import`](https://github.com/c-w/ghp-import)
```bash
ghp-import -n -p -f ./html/vorts
# -n  add .nojekyll
# -p  push (to origin/gh-pages by default)
# -f  force push
```
