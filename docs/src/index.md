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

To parallelize your estimation please remember to use `-t/--threads` option when starting Julia. Set `-t/threads` to `auto`  to set the threads to the number of CPU threads.You also can export the environment variable `JULIA_NUM_THREADS` into your `.bashrc` to avoid thread indication when starting Julia. 

```bash
julia -t 8
# or
julia -t auto
# or
printf "\nexport JULIA_NUM_THREADS=8\n" >> ~/.bashrc && source ~/.bashrc
```

To check the number of avaibles threads you can run the command `Threads.nthreads()` into your Julia session

Julia still have some latency problem loading packages. Please consider to create a custom [sysimages](https://github.com/JuliaLang/PackageCompiler.jl) for reduced latency. 

```bash 
curl -o mktest/julia_dependencies.jl https://raw.githubusercontent.com/jmurga/MKtest.jl/master/scripts/precompile_mktest.jl 
julia -e 'using PackageCompiler;PackageCompiler.create_sysimage(:MKtest, sysimage_path="/mktest/mktest.so", precompile_execution_file="mktest/precompile_mktest.jl")'
```

You can run the sysimage using `-J` option when starting Julia

```bash
julia -t8 -J mktest/mktest.so
```


### ABCreg
We have linked [ABCreg](https://github.com/molpopgen/ABCreg) with Julia to perform ABC inference. Nonetheless others ABC softwares could be used (like [abc (R package)](https://doi.org/10.1111/j.2041-210X.2011.00179.x) or [ABCToolBox](https://doi.org/10.1186/1471-2105-11-116)). If you are going to use ABCreg to directly make inference from our software please [cite the publication](https://doi.org/10.1186/1471-2156-10-35)


## Docker installation
The Docker image is based on Debian and includes all the software needed to run the pipeline. You can access Julia or Jupyter pulling the image from [Docker Hub](https://hub.docker.com/r/jmurga/mktest) or build yourself using the [Dockerfile on github](https://github.com/jmurga/MKtest.jl/blob/main/scripts/docker/Dockerfile). 

By default the Docker image will use any available CPU thread in your machine. Please, remember to use ```julia -t``` to use multiple threads and parallelize the estimation. To save the results, you can link the folder ```/mktest``` with any folder in your ```${HOME}``` directory.

```bash
MYPATH="/home/jmurga/temporal_mktest/"
# Pull the image
docker pull jmurga/mktest

# Runninh temporal docker container linking some local volume to export data. Consider to create a container.
docker run -it -v --rm ${MYPATH}:/mktest/folder jmurga/mktest abcmk_cli.jl
# Using -t 8 to parallelize using 8 threads
# Using -J mktest.so to avoid loading package latency
docker run -it -v --rm ${MYPATH}:/mktest/folder jmurga/mktest julia -J mktest.so -t8 abcmk_cli.jl

# Run jupyter notebook from docker image. Change the port if 8888 is already used
docker run -it --rm -v ${MYPATH}:/mktest/folder -p 8888:8888 jmurga/mktest /bin/bash -c "jupyter-lab --ip='*' --port=8888 --no-browser --allow-root"
```

To use our command-line interface, just run

```bash
docker run --rm -it -v ${MYPATH}:/mktest/ jmurga/mktest
```

## Singularity installation
The Docker image is based on Debian and includes all the software needed to run the pipeline. You can access Julia pulling the image from [Sylabs Cloud Library](https://cloud.sylabs.io/library/jmurga/default/mktest). 

```bash
singularity pull mktest.sif library://jmurga/default/mktest

singularity run mktest.sif julia ulia -J mktest.so -t8 abcmk_cli.jl
```