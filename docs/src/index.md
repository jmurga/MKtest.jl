# MKtest.jl

MKtest.jl is a Julia package including a fast Approximate Bayesian Computation version of the McDonald-Kreitman test (ABC-MK) presented in [Uricchio et al. (2019)](https://doi.org/10.1038/s41559-019-0890-6). The new ABC-MK implementation significantly improves the efficiency of the population genetics inferences. Following [Uricchio et al.(2019)](https://doi.org/10.1038/s41559-019-0890-6), the analytical estimations were used to explore the effect of background selection and selective interference on weakly beneficial alleles. Nonetheless, we developed a more straightforward and computationally efficient ABC-based inference procedure that accounts for the DFE of deleterious and beneficial alleles and partial recombination between selected genomic elements. Our approach estimates $\alpha$, $\alpha_W$, $\alpha_S$, and the Gamma distribution DFE parameters. 

In addition, the package automatizes other MK-like analyses parsing polymorphic and divergence data as well as including several extensions such as [Grapes](https://doi.org/10.1371/journal.pgen.1005774), [aMK](https://doi.org/10.1073/pnas.1220835110), [imputedMK](https://doi.org/10.1093/g3journal/jkac206) or [fwwMK](https://doi.org/10.1038/4151024a).



## Scratch installation
To install our module from scratch, we highly recommend using [LTS official Julia binaries](https://julialang.org/downloads/). Once you have installed Julia in your system, consider to install some important dependencies to automatize your pipelines. We have prepared a file to install them.

```bash
curl -o julia_dependencies.jl https://raw.githubusercontent.com/jmurga/MKtest.jl/master/scripts/julia_dependencies.jl
julia julia_dependencies.jl
```

You can easily install our Julia package manually executing

```bash
julia -e 'using Pkg;Pkg.add(PackageSpec(path="https://github.com/jmurga/MKtest.jl"))'
```

Or from Pkg REPL (by pressing `]` at Julia interpreter):

```julia
add https://github.com/jmurga/MKtest.jl
```

### ABCreg
We have linked [ABCreg](https://github.com/molpopgen/ABCreg) with Julia to perform ABC inference. Nonetheless others ABC softwares could be used (like [abc (R package)](https://doi.org/10.1111/j.2041-210X.2011.00179.x) or [ABCToolBox](https://doi.org/10.1186/1471-2105-11-116)). If you are going to use ABCreg to directly make inference from our software please [cite the publication](https://doi.org/10.1186/1471-2156-10-35) and compile it in your system. Anyway, once you get the summary statistic files you can use any other ABC software.

ABCreg needs *GSL* and *libz* to work. Please install both libraries before compile the software:

```bash
# Linux debian-based installation
sudo apt install libgsl-dev libz-dev build-essential git
# MacOS installation
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install GSL zlib git
```

```bash
git clone https://github.com/molpopgen/ABCreg.git ABCreg
cd ABCreg/src && make
```

## Docker installation
The Docker image is based on Debian and includes all the software needed to run the pipeline. You can access Julia or Jupyter pulling the image from [Docker Hub](https://hub.docker.com/r/jmurga/mktest). 

Please, remember to use ```julia -t``` to use multiple threads and parallelize the estimation. Although optional, note that you can run Julia using a pre-compiled image inside docker to avoid package loading latency. To save the results, you can link the folder ```/analysis``` with any folder in your ```${HOME}``` directory.

```bash
MYPATH="/home/jmurga/temporal_analysis/"
# Pull the image
docker pull jmurga/mktest

# Runninh temporal docker container linking some local volume to export data. Consider to create a container.
# Using -t 8 to parallelize using 8 threads
# Using -J mktest.so to avoid loading package latency
docker run -it -v --rm ${MYPATH}:/analysis/folder jmurga/mktest julia -t 8 -J mktest.so

# Run jupyter notebook from docker image. Change the port if 8888 is already used
docker run -it --rm -v ${MYPATH}:/analysis/folder -p 8888:8888 jmurga/mktest /bin/bash -c "jupyter-lab --ip='*' --port=8888 --no-browser --allow-root"
```

To use our command-line interface, just run

```bash
docker run -it -v ${MYPATH}:/analysis/ jmurga/mktest julia -t 8 -J mktest.so /analysis/abcmk_cli.jl
```