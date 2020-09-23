# ECHO RESONANCE Microbiome Paper 1

| | |
|------------:|:----------|
| **Authors** | [![Kevin Bonham][kevin-badge]][kevin-url] [![Vanja Klepac-Ceraj][vanja-badge]][vanja-url] |
| **DOIs**    | Repo: [![repo][repo-badge]][repo-url] Data: [![data][data-badge]][data-url] |

[kevin-badge]: https://img.shields.io/badge/Author-Kevin%20Bonham%2C%20PhD-blueviolet
[kevin-url]: http://nequals.me
[vanja-badge]: https://img.shields.io/badge/Author-Vanja%20Klepec--Ceraj%2C%20PhD-blueviolet
[vanja-url]: https://www.vkclab.com/
[vanja-badge]: https://img.shields.io/badge/Author-Vanja%20Klepec--Ceraj%2C%20PhD-blueviolet
[vanja-url]: https://www.vkclab.com/
[repo-badge]: https://zenodo.org/badge/222533623.svg
[repo-url]: https://zenodo.org/badge/latestdoi/222533623
[data-badge]: https://zenodo.org/badge/DOI/10.5281/zenodo.3633793.svg
[data-url]: https://doi.org/10.5281/zenodo.3633793

## Installation Instructions

All analysis code for this project is written in [`julia`][1],
(tested with version 1.3 on mac and linux).
In order to run it, follow the instructions below.

### Install julia

First, download the appropriate binaries for your system
from the [julia website][2]
and follow the installation instructions found there.
Code has been tested against julia v1.5.

Once completed,
you should be able to execute `julia` from the command prompt
to open the julia REPL.

```raw
$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.6.0-DEV.917 (2020-09-14)
 _/ |\__'_|_|_|\__'_|  |  Commit a0a68a54d6 (8 days old master)
|__/                   |


julia>
```

Type `exit()` and then `enter` from the julia REPL
to get back to your terminal's command prompt.

### Clone or download this repository

If you have `git` installed,
the easiest way to obtain this code is to clone the repository.

```sh
$ git clone https://github.com/Klepac-Ceraj-Lab/resonance_paper1.git
Cloning into 'resonance_paper1'...
remote: Enumerating objects: 289, done.
remote: Counting objects: 100% (289/289), done.
remote: Compressing objects: 100% (175/175), done.
remote: Total 289 (delta 184), reused 213 (delta 109), pack-reused 0
Receiving objects: 100% (289/289), 9.01 MiB | 9.26 MiB/s, done.
Resolving deltas: 100% (184/184), done.
$ cd resonance_paper1
$ git checkout v0.1.0
```

The final command checks out the version of the repository
from the initial submission.
You may also use the `master` branch,
but this is not guaranteed to work.
That said, if you find any problems with either the release
or with the master branch,
please [open an issue][3] and I will attempt to solve it as soon as possible.

Alternatively, the current release of these analysis notebooks
can be downloaded from the [releases page][4].
Simply download and unpack the archive file.

TODO: add correct link once release is made

```sh
$ curl -o repo.tar.gz
#...
$ cd resonance_paper1
```

### Instantiate the `analysis/` directory

The root directory has two files that enable easy replication
of the julia project environment and its dependencies,
`Project.toml` and `Manifest.toml`.
To use them, start a julia REPL,
activate the environment, and `instantiate`:

```raw
$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.6.0-DEV.917 (2020-09-14)
 _/ |\__'_|_|_|\__'_|  |  Commit a0a68a54d6 (8 days old master)
|__/                   |

[ Info: No project file, using default
[ Info: Loading OhMyREPL
[ Info: loading Revise
julia> using Pkg

julia> Pkg.activate(".")
Activating environment at `~/repos/lab/testing/resonance_paper1/Project.toml`

julia> Pkg.instantiate()
```

Then, take a look at the files in the `notebooks/` folder
to get started.
The `.jl` files are julia files, written in [Literate.jl][5] style.
are meant to be viewed in order.
All code is executed with the working directory set
to the root of this repository.

[1]: http://julialang.org
[2]: https://julialang.org/downloads/
[3]: https://github.com/Klepac-Ceraj-Lab/resonance_paper1/issues
[4]: https://github.com/Klepac-Ceraj-Lab/resonance_paper1/releases
[5]: https://fredrikekre.github.io/Literate.jl/stable/
