# =============================================================================
# DFEResult: what dfe_ml returns
#
# Fixed-schema long tables keyed by `fit` (one id per fitted model x dataset x settings),
# so results of different models, datasets or calls combine with `vcat` and join on `fit`
# (or on the readable `dataset` / `model` columns). Engine-specific tables go in `extra`.
# =============================================================================

"""
    DFEResult

The result of `dfe_ml` (and of the fitting form of `dfe_model_comparison`). Tables, all
`DataFrame`s keyed by the integer `fit` with readable `dataset` and `model` columns:

- `inference`: one row per fit. `model`, `config` (row of `config`), `n` (chromosomes),
  `converged`, `loglik`, `n_params`, `theta`, then the summary columns (`DFEMAlpha()`:
  `alpha`, `alpha_ci`, `omega_a`, `omega_a_ci`, `omega_na`, `omega_na_ci`, `omega_obs`,
  `alpha_div`, `alpha_dfe`). Intervals are `(lo, hi)` tuples; NaN when not computed.
- `params`: one row per fitted parameter: `parameter`, `kind` (`"dfe"` or `"nuisance"`),
  `value`, `lower`, `upper`, `free`, `at_bound` (a free parameter on its bound is not an
  estimate). The same columns for every model, so models of any shape stack.
- `sfs`: one row per SFS bin: `j`, `pn_obs`, `ps_obs`, `pn_exp`, `ps_exp` (fitted
  expectation), `r_j` (the per-bin nuisance, 1 under `Equilibrium()` and for `j = 1`).
- `divergence`: one row per fit: `Dn_obs`, `Ds_obs`, `Dn_total` (alpha's denominator: the
  observed Dn under `SaturatedDivergence`), `Dn_pred` (the DFE's prediction), `lambda`
  (neutral substitutions per site), and the expected Dn by S class, GRAPES' `prop_subst_*`:
  `Dn_sdel` (S <= -nearly_neutral), `Dn_wdel` (up to 0), `Dn_wadv` (up to nearly_neutral),
  `Dn_sadv` (above). `alpha = Dn_sadv / Dn_total`.
- `config`: one row per distinct setting: the axis settings, `nearly_neutral`,
  `alpha_ci_threshold`, `ci_method`, `nb_rand_start`, `fixed`, `bounds`, `fold`,
  `project_sfs`, `mktest_version`.
- `extra`: `OrderedDict{Symbol,DataFrame}` of engine-specific tables, same keys (e.g.
  `:exterior` when the expectation adds fitdadi's exterior masses).

Combine with `vcat(r1, r2, ...)`; save with `write_dfe_result` / `read_dfe_result` (CSV)
or `JLD2.save`. Accessors: `params_wide`, `param_values`, `expected_sfs`, `rj`.
"""
struct DFEResult
    inference::DataFrame
    params::DataFrame
    sfs::DataFrame
    divergence::DataFrame
    config::DataFrame
    extra::OrderedDict{Symbol,DataFrame}
end

# The settings that define the likelihood function: fits that share them (and the data)
# can be compared by likelihood. Optimizer, starts, fixed values and bounds do not change
# the function, only where on it a fit ends.
const _LIKELIHOOD_SETTINGS = (:expectation, :integrator, :quadrature, :s_grid, :dominance,
                              :demography, :ne_model, :div_model, :kernel, :fold, :project_sfs)

_key_cols(i, dataset, model) = (fit = i, dataset = dataset, model = model)

function _assemble_result(pieces::Vector, datasets::AbstractVector, config::OrderedDict{Symbol,Any})
    inf = DataFrame[]
    par = DataFrame[]
    sfs = DataFrame[]
    dvg = DataFrame[]
    ext = OrderedDict{Symbol,Vector{DataFrame}}()
    for (i, p) in enumerate(pieces)
        ds = datasets[i]
        model = p.inference[:model]
        row = OrderedDict{Symbol,Any}(:fit => i, :dataset => ds, :model => model, :config => 1)
        for (k, v) in p.inference
            k === :model && continue
            row[k] = v
        end
        push!(inf, DataFrame([NamedTuple{Tuple(keys(row))}(Tuple(values(row)))]))
        push!(par, DataFrame([merge(_key_cols(i, ds, model), r) for r in p.params]))
        nb = length(p.sfs.j)
        push!(sfs, DataFrame(merge((fit = fill(i, nb), dataset = fill(ds, nb), model = fill(model, nb)), p.sfs)))
        push!(dvg, DataFrame([merge(_key_cols(i, ds, model), p.divergence)]))
        for (name, t) in p.extra
            push!(get!(ext, name, DataFrame[]), DataFrame([merge(_key_cols(i, ds, model), t)]))
        end
    end
    cfg = DataFrame([NamedTuple{(:config, keys(config)...)}((1, values(config)...))])
    return DFEResult(vcat(inf...; cols = :union), vcat(par...; cols = :union), vcat(sfs...; cols = :union),
                     vcat(dvg...; cols = :union), cfg,
                     OrderedDict{Symbol,DataFrame}(k => vcat(v...; cols = :union) for (k, v) in ext))
end

"""
    vcat(r1::DFEResult, r2::DFEResult, ...) -> DFEResult

Stack results: `fit` ids are renumbered so they stay unique, and config rows with identical
settings are merged into one id.
"""
function Base.vcat(r1::DFEResult, rs::DFEResult...)
    all = (r1, rs...)
    tabs = Dict(k => DataFrame[] for k in (:inference, :params, :sfs, :divergence))
    ext = OrderedDict{Symbol,Vector{DataFrame}}()
    cfg = DataFrame()
    fo = 0
    for r in all
        # config ids: reuse an identical row already collected
        remap = Dict{Int,Int}()
        for c in eachrow(r.config)
            settings = NamedTuple(c[Not(:config)])
            hit = nrow(cfg) == 0 ? nothing :
                  findfirst(e -> isequal(NamedTuple(e[Not(:config)]), settings), eachrow(cfg))
            if hit === nothing
                nid = nrow(cfg) + 1
                cfg = vcat(cfg, DataFrame([merge((config = nid,), settings)]); cols = :union)
                remap[c.config] = nid
            else
                remap[c.config] = cfg.config[hit]
            end
        end
        for k in (:inference, :params, :sfs, :divergence)
            t = copy(getfield(r, k))
            t.fit .+= fo
            k === :inference && (t.config = [remap[c] for c in t.config])
            push!(tabs[k], t)
        end
        for (name, t) in r.extra
            t2 = copy(t)
            t2.fit .+= fo
            push!(get!(ext, name, DataFrame[]), t2)
        end
        fo += nrow(r.inference) == 0 ? 0 : maximum(r.inference.fit)
    end
    cat(k) = vcat(tabs[k]...; cols = :union)
    return DFEResult(cat(:inference), cat(:params), cat(:sfs), cat(:divergence), cfg,
                     OrderedDict{Symbol,DataFrame}(k => vcat(v...; cols = :union) for (k, v) in ext))
end

function Base.show(io::IO, ::MIME"text/plain", r::DFEResult)
    n = nrow(r.inference)
    models = unique(r.inference.model)
    println(io, "DFEResult: $n fit(s), $(length(unique(r.inference.dataset))) dataset(s), models ",
            join(models, ", "), ", $(nrow(r.config)) config(s)")
    for k in (:inference, :params, :sfs, :divergence, :config)
        t = getfield(r, k)
        println(io, "  .", rpad(string(k), 11), nrow(t), " x ", ncol(t))
    end
    isempty(r.extra) || print(io, "  .extra      ", join(keys(r.extra), ", "))
    if n > 0 && hasproperty(r.inference, :alpha)
        println(io)
        show(io, select(r.inference, intersect([:fit, :dataset, :model, :converged, :loglik, :alpha, :alpha_ci],
                                               propertynames(r.inference))); summary = false, eltypes = false)
    end
end
Base.show(io::IO, r::DFEResult) = print(io, "DFEResult($(nrow(r.inference)) fits)")

"""
    params_wide(r::DFEResult) -> DataFrame

One row per fit (`fit`, `dataset`, `model`), one column per parameter (`missing` where a
model does not have it).
"""
params_wide(r::DFEResult) = unstack(r.params, [:fit, :dataset, :model], :parameter, :value)

"""
    param_values(r::DFEResult, name) -> Vector

The value of parameter `name` for every fit, in `r.inference` order (`missing` where a
model does not have it).
"""
function param_values(r::DFEResult, name::Union{Symbol,AbstractString})
    s = String(name)
    d = Dict(row.fit => row.value for row in eachrow(r.params) if row.parameter == s)
    return [get(d, f, missing) for f in r.inference.fit]
end

"""
    expected_sfs(r::DFEResult, fit = 1) -> Matrix{Float64}

`[j pn_obs ps_obs pn_exp ps_exp]` of one fit.
"""
function expected_sfs(r::DFEResult, fit::Integer = first(r.inference.fit))
    t = r.sfs[r.sfs.fit .== fit, :]
    return Matrix{Float64}(t[:, [:j, :pn_obs, :ps_obs, :pn_exp, :ps_exp]])
end

"""
    rj(r::DFEResult, fit = 1) -> Vector{Float64}

The fitted per-bin nuisance `r_j` of one fit.
"""
rj(r::DFEResult, fit::Integer = first(r.inference.fit)) = r.sfs.r_j[r.sfs.fit .== fit]

# --- CSV I/O -------------------------------------------------------------------

_istuplecol(v) = eltype(v) <: Tuple || (eltype(v) <: Union{Missing,Tuple} && any(x -> x isa Tuple, v))

function _split_intervals(t::DataFrame)
    out = DataFrame()
    for c in propertynames(t)
        v = t[!, c]
        if _istuplecol(v)
            out[!, Symbol(c, :_low)] = [ismissing(x) ? missing : x[1] for x in v]
            out[!, Symbol(c, :_high)] = [ismissing(x) ? missing : x[2] for x in v]
        else
            out[!, c] = v
        end
    end
    return out
end

function _join_intervals(t::DataFrame)
    out = DataFrame()
    cols = propertynames(t)
    for c in cols
        s = string(c)
        if endswith(s, "_low") && Symbol(s[1:end-4], :_high) in cols
            base = Symbol(s[1:end-4])
            out[!, base] = [(a, b) for (a, b) in zip(t[!, c], t[!, Symbol(base, :_high)])]
        elseif endswith(s, "_high") && Symbol(s[1:end-5], :_low) in cols
            continue
        else
            out[!, c] = t[!, c]
        end
    end
    return out
end
