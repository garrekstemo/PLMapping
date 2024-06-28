module PLMapping

export
    lorentzian,
    trapezoid,
    find_maxima,
    flip_rows,
    find_area

using DelimitedFiles
using LsqFit
using Peaks
using Optim
using SignalFiltering

include("analysis.jl")
include("cleaning.jl")


end # module PLMapping
