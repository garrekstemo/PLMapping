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
    x = spectrum[:, 1]
    y = spectrum[:, 2]
    y = median_filter(y, 5)
    y = y .- background
    if y == zeros(length(y))
        return area
    end

    pks, vals = findmaxima(y, 40)
    center_guess = x[pks[argmax(vals)]]

    p0 = [maximum(y), center_guess, 1]
    if  center_guess > (λ_0 - 50) && center_guess < (λ_0 + 50)
        fit = curve_fit(lorentzian, x, y, p0)
        area = trapezoid(x, lorentzian(x, fit.param))
    end
    return area, fit
end
