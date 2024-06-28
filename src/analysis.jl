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
    y = median_filter(y .- background, 5)
    if y == zeros(length(y))
        return area, fit, y
    end

    y_max = maximum(y)
    γ_guess = 15.0
    p0 = [2000.0, λ_0, γ_guess]
    fit = optimize(b -> squared_error(b, x, y), p0)
    params = Optim.minimizer(fit)
    if discard_poor_fit(params, λ_0)
        fit = nothing
    else
        area = trapezoid(x, lorentzian(x, params))
    end
    return area, fit, y
end

function discard_poor_fit(params, center)
    if params[1] < 0 || params[2] < center - 70 || params[2] > center + 70 || params[3] < 0 || params[3] > 50
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