FROM julia

RUN apt-get update\
    && apt-get install -y libglib2.0-0 libxext6 wget bzip2 ca-certificates curl git vim make build-essential libgsl-dev libz-dev gzip parallel && \
    apt-get clean  && \
    rm -rf /var/lib/apt/lists/*

RUN useradd -ms /bin/bash mktest && usermod -a -G sudo mktest

RUN chown -R mktest:mktest /usr/local/julia/

USER mktest

WORKDIR /home/mktest/

ADD data/* /home/mktest/
ADD scripts/abcmk_cli.jl /home/mktest/abcmk_cli.jl
ADD scripts/precompile_mktest.jl /home/mktest/precompile_mktest.jl

RUN julia -e 'using Pkg;Pkg.add(["CondaPkg","Fire","JLD2", "DataFrames", "CSV","StatsBase"]);Pkg.add(PackageSpec(path="https://github.com/jmurga/MKtest.jl"))'

RUN julia -e 'using CondaPkg;CondaPkg.add("grapes-static", channel = "genomedk");CondaPkg.add("bedtools", channel = "bioconda")'

RUN julia -e 'using CondaPkg,Pkg;CondaPkg.add.(["r-essentials","r-locfit","r-sparsem","r-quantreg","r-nnet"], channel = "conda-forge");r = CondaPkg.which("R");ENV["R_HOME"] = replace(r,"bin"=>"lib");Pkg.add("RCall");using RCall;R"""install.packages("abc",repos="http://cran.r-project.org",dependencies = TRUE)"""'

RUN git clone https://github.com/molpopgen/ABCreg.git /home/mktest/ABCreg
RUN cd /home/mktest/ABCreg/src && make

ENV PATH /home/mktest/ABCreg/src:$PATH

CMD ["julia"]
