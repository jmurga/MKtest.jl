FROM julia

RUN apt-get update\
    && apt-get install -y libglib2.0-0 libxext6 wget bzip2 ca-certificates curl git vim make build-essential libgsl-dev libz-dev gzip parallel && \
    apt-get clean  && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /analysis

ADD data/tgp.txt /analysis/tgp.txt
ADD data/dgn_ral.txt /analysis/dgn_ral.txt
ADD data/dgn_zi.txt /analysis/dgn_zi.txt

RUN julia -e 'using Pkg;Pkg.add(["IJulia","CSV","DataFrames"])'

RUN julia -e 'using Pkg;Pkg.add(PackageSpec(path="https://github.com/jmurga/MKtest.jl"))'

ADD scripts/abcmk_cli.jl /analysis/abcmk_cli.jl

RUN julia -e 'using PackageCompiler;PackageCompiler.create_sysimage(:MKtest, sysimage_path="/analysis/mktest.so", precompile_execution_file="/analysis/precompile_mktest.jl")'

RUN git clone https://github.com/molpopgen/ABCreg.git /analysis/ABCreg
RUN cd /analysis/ABCreg/src && make

ENV PATH /analysis/ABCreg/src:$PATH

WORKDIR /analysis

CMD ["julia","-J","mktest.so"]
