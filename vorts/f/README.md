
`m_vorts` uses some [Fortran 2003 features](http://fortranwiki.org/fortran/show/Fortran+2003).

## Compiling

Requires:
- Meson
- Ninja
- a Fortran compiler

```
meson setup _build
meson compile -C _build
```

### Windows

You can download [gfortran binaries for Windows](https://gcc.gnu.org/wiki/GFortranBinariesWindows),
but in order to also have `make`, MinGW via [MSYS2](https://www.msys2.org/) or [WSL](https://docs.microsoft.com/en-us/windows/wsl/about) might be an easier option.
In order to run the MinGW-compiled executable from Python on Windows, MinGW tools have to be on your PATH.

### Linux

On Debian-like distros, `gfortran` can be obtained using `apt`.

On Linux or Mac, [Homebrew](https://brew.sh/) can be used to (potentially) obtain more recent versions of GNU Fortran
by installing the GNU compiler collection formula: `brew install gcc` (see [formulae](https://formulae.brew.sh/formula/gcc)).
