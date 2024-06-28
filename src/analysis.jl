function lorentzian(ω, p)
    A, ω₀, γ = p
    return @. A * γ / ((ω - ω₀)^2 +  γ^2) / π
end

function trapezoid(x, y)
    area = 0
    i = 1
    while i < length(y)
        area += (y[i] + y[i+1]) * (x[i+1] - x[i]) / 2.0
        i += 1
    end
    return area
end

function find_maxima(dir)

    max_vals = []
    spectra = []
    files = readdir(dir)

    for file in files
        filepath = joinpath(dir, file)
        spec = readdlm(filepath, skipstart=1)
        x = spec[:, 1]
        y = spec[:, 2]
        max = findfirst(x -> x == maximum(y), y)
        push!(spectra, spec)
        push!(max_vals, max)
    end

    return spectra, max_vals
end

function find_area(spectrum, background, λ_0)
    area = 0.0
    fit = nothing
    x = spectrum[:, 1]
    y = spectrum[:, 2]
    y = y .- background
    y = median_filter(y, 5)
    if y == zeros(length(y))
        return area, fit, y
    end

    p0 = [5000.0, λ_0, 15.0]
    # fit = curve_fit(lorentzian, x, y, p0)
    fit = optimize(b -> squared_error(b, x, y), p0)
    params = Optim.minimizer(fit)
    if discard_poor_fit(params, λ_0)
        fit = nothing
    end
    area = trapezoid(x, lorentzian(x, params))
    # if discard_poor_fit(fit, λ_0)
    #     fit = nothing
    # end
    return area, fit, y
end

function discard_poor_fit(params, center)
    if params[2] < center - 50 || params[2] > center + 50 || params[3] < 0
        return true
    end
    return false
end

function squared_error(p, X, Y)
    error = 0.0
    for i in eachindex(X)
        pred_i = lorentzian(X[i], p)
        error += (Y[i] - pred_i)^2
    end
    return error
end