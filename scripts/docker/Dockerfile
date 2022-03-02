FROM julia

RUN apt-get update\
    && apt-get install -y libglib2.0-0 libxext6 wget bzip2 ca-certificates curl git vim make build-essential libgsl-dev libz-dev gzip parallel r-base && \
    apt-get clean  && \
    rm -rf /var/lib/apt/lists/*


ENV PATH /opt/conda/bin:$PATH

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

RUN mkdir -p /analysis
RUN conda install -c conda-forge mamba -y
RUN mamba install -c conda-forge jupyterlab r-essentials r-r.utils r-locfit numpy numba scipy -y

RUN conda init bash 

RUN mkdir -p /analysis

ADD data/tgp.txt /analysis/tgp.txt
ADD data/dgn_ral.txt /analysis/dgn_ral.txt
ADD data/dgn_zi.txt /analysis/dgn_zi.txt

RUN julia -e 'using Pkg;Pkg.add(["Fire","IJulia","PackageCompiler","CSV","JLD2","DataFrames","RCall","Plots","StatsBase","PyCall"])'

RUN julia -e 'using Pkg;Pkg.add(PackageSpec(path="https://github.com/jmurga/MKtest.jl"))'

ADD scripts/precompile_mktest.jl /analysis/precompile_mktest.jl
ADD scripts/abcmk_cli.jl /analysis/abcmk_cli.jl

RUN julia -e 'using PackageCompiler;PackageCompiler.create_sysimage(:MKtest, sysimage_path="/analysis/mktest.so", precompile_execution_file="/analysis/precompile_mktest.jl")'

RUN git clone https://github.com/molpopgen/ABCreg.git /analysis/ABCreg
RUN cd /analysis/ABCreg/src && make

ENV PATH /analysis/ABCreg/src:$PATH

WORKDIR /analysis

CMD ["julia","-J","mktest.so"]
