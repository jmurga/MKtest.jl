FROM julia

ENV PATH /root/mambaforge/bin:$PATH

RUN apt-get update\
    && apt-get install -y libglib2.0-0 libxext6 wget bzip2 ca-certificates curl git vim make build-essential libgsl-dev libz-dev gzip parallel && \
    apt-get clean  && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /mktest

RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
RUN bash Mambaforge-$(uname)-$(uname -m).sh -b

ENV LD_LIBRARY_PATH /root/mambaforge/lib:$LD_LIBRARY_PATH

RUN mamba update -n base -c defaults conda
RUN mamba init bash

RUN mamba install -c conda-forge jupyter -y

ADD data/* mktest/
ADD scripts/abcmk_cli.jl /mktest/abcmk_cli.jl
ADD scripts/precompile_mktest.jl /mktest/precompile_mktest.jl

RUN julia -e 'using Pkg;Pkg.add(["CondaPkg","CSV","DataFrames","PackageCompiler","Fire"])'

RUN julia -e 'using Pkg;Pkg.add(PackageSpec(path="https://github.com/jmurga/MKtest.jl"))'
RUN julia -e 'using CondaPkg;CondaPkg.add("grapes-static", channel = "genomedk");CondaPkg.add("bedtools", channel = "bioconda")'

RUN julia -e 'using PackageCompiler;PackageCompiler.create_sysimage(:MKtest, sysimage_path="/mktest/mktest.so", precompile_execution_file="/mktest/precompile_mktest.jl")'

RUN git clone https://github.com/molpopgen/ABCreg.git /mktest/ABCreg
RUN cd /mktest/ABCreg/src && make

ENV PATH /mktest/ABCreg/src:$PATH

WORKDIR /mktest

ENV NTHREADS=1

CMD ["julia","-J","mktest.so","-t","auto"]