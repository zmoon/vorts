
Build documentation from in here (`doc`) by running
```bash
pdoc --html --force --template-dir ./templates ../vorts
```

This generates the site in `./html/vorts`.

Alternatively, use
```bash
pdoc --http : --template-dir ./doc/templates ./vorts
```
to run a local HTTP server.
