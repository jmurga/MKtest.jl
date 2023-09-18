_projection_cache = {}
"""
Log of N choose k.
"""
function _lncomb(N::Int64,k::Union{UnitRange{Int64},StepRange{Int64, Int64}})
    out = zeros(length(k))
    for (i,v) ∈ enumerate(k)
        out[i] = loggamma(N+1) - loggamma(v+1) - loggamma(N-v+1)
    end
    return out
end


function _lncomb(N::Int64,k::Int64)
    return loggamma(N+1) - loggamma(k+1) - loggamma(N-k+1)
end

"""
Coefficients for projection from a different fs size.

proj_to: Numper of samples to project down to.
proj_from: Numper of samples to project from.
hits: Number of derived alleles projecting from.
"""
def _cached_projection!(proj_to::Int, proj_from::Int, hits,_projection_cache::Dict(Tuple,Float64)):
    key = (proj_to, proj_from, hits)
    @assert proj_to < proj_from

    # We set numpy's error reporting so that it will ignore underflows, 
    # because those just imply that contrib is 0.
    proj_hits = 0:proj_to
    # For large sample sizes, we need to do the calculation in logs, and it
    # is accurate enough for small sizes as well.
    lncontrib = _lncomb(proj_to,proj_hits)
    @. lncontrib += _lncomb(proj_from-proj_to,hits-proj_hits)
    @. lncontrib -= _lncomb(proj_from, hits)
    contrib = sum(@. exp(lncontrib))
    _projection_cache[key] = contrib
    return contrib

def project(s, ns):
    """
    Project to smaller sample size.

    ns: Sample sizes for new spectrum.
    """
    if len(ns) != s.Npop:
        raise ValueError('Requested sample sizes not of same dimension '
                         'as spectrum. Perhaps you need to marginalize '
                         'over some populations first?')
    if np.any(np.asarray(ns) > np.asarray(s.sample_sizes)):
        raise ValueError('Cannot project to a sample size greater than '
                         'original. Original size is %s and requested size '
                         'is %s.' % (s.sample_sizes, ns))


    output = s.copy()

    if proj != s.sample_sizes[0]:
        output = output._project_one_axis(proj, axis)

    output.pop_ids = s.pop_ids
    output.extrap_x = s.extrap_x



def _project_one_axis(s, n, axis=0):
    """
    Project along a single axis.
    """
    # This gets a little tricky with fancy indexing to make it work
    # for fs with arbitrary number of dimensions.
    if n > s.sample_sizes[axis]:
        raise ValueError('Cannot project to a sample size greater than original. Called sizes were from %s to %s.'  % (s.sample_sizes[axis], n))

    newshape = list(s.shape)
    newshape[axis] = n+1
    # Create a new empty fs that we'll fill in below.
    pfs = Spectrum(np.zeros(newshape), mask_corners=False)

    # Set up for our fancy indexes. These slices are currently like
    # [:,:,...]
    from_slice = [slice(None) for ii in range(s.Npop)]
    to_slice = [slice(None) for ii in range(s.Npop)]
    proj_slice = [nuax for ii in range(s.Npop)]

    proj_from = s.sample_sizes[axis]
    # For each possible number of hits.
    for hits=0:proj_from
        # Adjust the slice in the array we're projecting from.
        from_slice[axis] = slice(hits, hits+1)
        # These are the least and most possible hits we could have in the
        #  projected fs.
        least, most = max(n - (proj_from - hits), 0), min(hits,n)
        to_slice[axis] = slice(least, most+1)
        # The projection weights.
        proj = _cached_projection!(n, proj_from, hits,_projection_cache)
        proj_slice[axis] = slice(least, most+1)
        # Do the multiplications
        pfs.data[tuple(to_slice)] += s.data[tuple(from_slice)] * proj[tuple(proj_slice)]
        pfs.mask[tuple(to_slice)] = np.logical_or(pfs.mask[tuple(to_slice)],s.mask[tuple(from_slice)])

    return pfs
