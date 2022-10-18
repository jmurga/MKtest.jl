FROM julia

ENV PATH /root/mambaforge/bin:$PATH

RUN apt-get update\
    && apt-get install -y libglib2.0-0 libxext6 wget bzip2 ca-certificates curl git vim make build-essential libgsl-dev libz-dev gzip parallel && \
    apt-get clean  && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /analysis

RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
RUN bash Mambaforge-$(uname)-$(uname -m).sh -b

ENV LD_LIBRARY_PATH /root/mambaforge/lib:$LD_LIBRARY_PATH

RUN mamba update -n base -c defaults conda
RUN mamba init bash

RUN /root/mambaforge/bin/mamba install -c conda-forge jupyter -y

ADD data/* analysis/
ADD scripts/abcmk_cli.jl /analysis/abcmk_cli.jl
ADD scripts/precompile_mktest.jl /analysis/precompile_mktest.jl

RUN julia -e 'using Pkg;Pkg.add(["IJulia","CSV","DataFrames","PackageCompiler","Fire"])'

RUN julia -e 'using Pkg;Pkg.add(PackageSpec(path="https://github.com/jmurga/MKtest.jl"))'

RUN julia -e 'using PackageCompiler;PackageCompiler.create_sysimage(:MKtest, sysimage_path="/analysis/mktest.so", precompile_execution_file="/analysis/precompile_mktest.jl")'

RUN git clone https://github.com/molpopgen/ABCreg.git /analysis/ABCreg
RUN cd /analysis/ABCreg/src && make

ENV PATH /analysis/ABCreg/src:$PATH

WORKDIR /analysis
