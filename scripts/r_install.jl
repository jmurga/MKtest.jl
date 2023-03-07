using CondaPkg,Pkg
function r_install()
    if (isnothing(CondaPkg.which("R")))
        
        CondaPkg.add.(["r-essentials","r-locfit","r-sparsem","r-quantreg","r-nnet"], channel = "conda-forge");
        r = CondaPkg.which("R")

        ENV["R_HOME"] = replace(r,"bin"=>"lib")
        Pkg.add("RCall")
        using RCall
        R"install.packages('abc',repos='http://cran.r-project.org',dependencies = TRUE)"

    end
end