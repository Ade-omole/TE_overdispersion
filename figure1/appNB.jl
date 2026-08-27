# Individual-based TE simulation with approximate and negative-binomial ODE closures.

using Random, Distributions, CSV, DataFrames, Statistics

# Configuration
const FITNESS_TYPE = "exp"  # "exp" => w(n) = exp(-s*n^2); "syn" => w(n) = 1 - s*n^2
const FREE_RECOMBINATION = true  # true: independent assortment; false: Poisson(R) crossovers

const n_init_param = 10
const population_size = 10^4
const generations = 30000
const R = 10.0^(-1)  # used only when FREE_RECOMBINATION = false

const u = 0.026
const v = 0.006
const s = 1 / population_size
const dt = 1.0

const UV_EVENT_MODEL = :binomial  # :poisson or :binomial

if UV_EVENT_MODEL == :poisson
    transposition_k(n::Int) = n > 0 ? rand(Poisson(u * n * dt)) : 0
    excision_k(n::Int) = n > 0 ? rand(Poisson(v * n * dt)) : 0
elseif UV_EVENT_MODEL == :binomial
    transposition_k(n::Int) = n > 0 ? rand(Binomial(n, clamp(u * dt, 0.0, 1.0))) : 0
    excision_k(n::Int) = n > 0 ? rand(Binomial(n, clamp(v * dt, 0.0, 1.0))) : 0
else
    error("UV_EVENT_MODEL must be :poisson or :binomial (got $UV_EVENT_MODEL)")
end

μ_app = 10.0
σ²_app = 10.0
μ_nb = 10.0
σ²_nb = 10.0
μ_c = 10.0

const PRINT_EVERY = 200
const WRITE_CSV = true

const te_root = dirname(@__DIR__)
const csv_dir = joinpath(te_root, "csv_files")
mkpath(csv_dir)

function create_individual()
    n_init = rand(Poisson(n_init_param))
    chrom1 = Float64[]
    chrom2 = Float64[]
    for _ in 1:n_init
        push!(rand(1:2) == 1 ? chrom1 : chrom2, rand())
    end
    return (chrom1, chrom2)
end

function apply_excision!(chrom1::Vector{Float64}, chrom2::Vector{Float64})
    n_total = length(chrom1) + length(chrom2)
    num_del_total = excision_k(n_total)
    if num_del_total <= 0 || n_total == 0
        return nothing
    end

    n1 = length(chrom1)
    n2 = length(chrom2)
    frac1 = n1 / n_total
    n1_del = min(rand(Binomial(num_del_total, frac1)), n1)
    n2_del = min(num_del_total - n1_del, n2)

    function delete_uniform_indices!(chrom::Vector{Float64}, d::Int)
        n = length(chrom)
        if d <= 0 || n == 0
            return nothing
        elseif d >= n
            empty!(chrom)
            return nothing
        end
        k = n - d
        sample_keep = k <= d
        m = sample_keep ? k : d
        selected = Set{Int}()
        sizehint!(selected, m)
        for j in (n - m + 1):n
            t = rand(1:j)
            if t in selected
                push!(selected, j)
            else
                push!(selected, t)
            end
        end
        write_idx = 0
        @inbounds for i in 1:n
            keep = sample_keep ? (i in selected) : !(i in selected)
            if keep
                write_idx += 1
                chrom[write_idx] = chrom[i]
            end
        end
        resize!(chrom, write_idx)
        return nothing
    end

    delete_uniform_indices!(chrom1, n1_del)
    delete_uniform_indices!(chrom2, n2_del)
    return nothing
end

function apply_transposition!(chrom1::Vector{Float64}, chrom2::Vector{Float64}, n_initial::Int)
    k = transposition_k(n_initial)
    for _ in 1:k
        push!(rand(1:2) == 1 ? chrom1 : chrom2, rand())
    end
    return nothing
end

function meiosis(chrom1::Vector{Float64}, chrom2::Vector{Float64})
    if FREE_RECOMBINATION
        g1 = Float64[]
        g2 = Float64[]
        sizehint!(g1, length(chrom1) + length(chrom2))
        sizehint!(g2, length(chrom1) + length(chrom2))
        @inbounds for p in chrom1
            (rand(Bool) ? push!(g1, p) : push!(g2, p))
        end
        @inbounds for p in chrom2
            (rand(Bool) ? push!(g1, p) : push!(g2, p))
        end
        return (g1, g2)
    end

    n_cross = rand(Poisson(R))
    if n_cross == 0
        return (chrom1, chrom2)
    end
    xovers = sort(rand(n_cross))
    g1 = Float64[]
    g2 = Float64[]
    sizehint!(g1, length(chrom1) + length(chrom2))
    sizehint!(g2, length(chrom1) + length(chrom2))
    @inbounds for p in chrom1
        seg = searchsortedfirst(xovers, p)
        isodd(seg) ? push!(g1, p) : push!(g2, p)
    end
    @inbounds for p in chrom2
        seg = searchsortedfirst(xovers, p)
        iseven(seg) ? push!(g1, p) : push!(g2, p)
    end
    return (g1, g2)
end

@inline function fitness_value(x_total::Int)
    if FITNESS_TYPE == "exp"
        return exp(-s * x_total^2)
    else
        return max(0.0, 1.0 - s * x_total^2)
    end
end

function form_offspring(parent1::Tuple{Vector{Float64}, Vector{Float64}},
                        parent2::Tuple{Vector{Float64}, Vector{Float64}})
    g1_p1, g2_p1 = meiosis(parent1[1], parent1[2])
    g1_p2, g2_p2 = meiosis(parent2[1], parent2[2])
    chrom1 = rand(1:2) == 1 ? copy(g1_p1) : copy(g2_p1)
    chrom2 = rand(1:2) == 1 ? copy(g1_p2) : copy(g2_p2)
    n0 = length(chrom1) + length(chrom2)
    apply_excision!(chrom1, chrom2)
    apply_transposition!(chrom1, chrom2, n0)
    return (chrom1, chrom2)
end

function generate_next_generation!(population::Vector{Tuple{Vector{Float64}, Vector{Float64}}},
                                  fitness::Vector{Float64},
                                  cum_probs::Vector{Float64})
    tot = 0.0
    @inbounds for i in 1:population_size
        p = population[i]
        n = length(p[1]) + length(p[2])
        fi = fitness_value(n)
        fitness[i] = fi
        tot += fi
    end
    tot = tot > 0 ? tot : 1.0

    acc = 0.0
    inv_tot = 1.0 / tot
    @inbounds for i in 1:population_size
        acc += fitness[i] * inv_tot
        cum_probs[i] = acc
    end
    cum_probs[end] = 1.0

    new_population = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, population_size)
    @inbounds for i in 1:population_size
        i1 = searchsortedfirst(cum_probs, rand())
        i2 = searchsortedfirst(cum_probs, rand())
        i1 = i1 > population_size ? population_size : i1
        i2 = i2 > population_size ? population_size : i2
        new_population[i] = form_offspring(population[i1], population[i2])
    end
    return new_population
end

function compute_sim_moments(population::Vector{Tuple{Vector{Float64}, Vector{Float64}}})
    sum_n = 0.0
    @inbounds for i in 1:population_size
        p = population[i]
        sum_n += (length(p[1]) + length(p[2]))
    end
    μ = sum_n / population_size

    sum_sq = 0.0
    @inbounds for i in 1:population_size
        p = population[i]
        d = (length(p[1]) + length(p[2])) - μ
        sum_sq += d * d
    end
    σ² = max(0.0, sum_sq / population_size)
    σ²_safe = max(eps(), σ²)

    sum_c3 = 0.0
    sum_c4 = 0.0
    @inbounds for i in 1:population_size
        p = population[i]
        d = (length(p[1]) + length(p[2])) - μ
        d2 = d * d
        sum_c3 += d2 * d
        sum_c4 += d2 * d2
    end
    m3 = sum_c3 / population_size
    m4 = sum_c4 / population_size
    skew = m3 / (σ²_safe^1.5)
    exkurt = (m4 / (σ²_safe^2)) - 3.0
    p_sim = μ / (σ² + eps())
    return μ, σ², skew, exkurt, p_sim, m3
end

function main()
    println("Simulation parameters:")
    println("  Initial TEs per individual: Poisson($n_init_param)")
    println("  Population size (N): $population_size")
    println("  Generations: $generations")
    println("  Transposition rate (u): $u, Excision rate (v): $v")
    println("  u/v event model: $UV_EVENT_MODEL")
    println("  Selection coefficient (s): $s")
    println("  Recombination: $(FREE_RECOMBINATION ? "free assortment" : "Poisson(R=$R) crossovers")")
    println("  Fitness: $(FITNESS_TYPE == "exp" ? "exp(-s * x_total^2)" : "1 - s * x_total^2")")
    println()

    μ_app_state = μ_app
    σ²_app_state = σ²_app
    μ_nb_state = μ_nb
    σ²_nb_state = σ²_nb
    μ_c_state = μ_c

    population = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, population_size)
    for i in 1:population_size
        population[i] = create_individual()
    end
    fitness = Vector{Float64}(undef, population_size)
    cum_probs = Vector{Float64}(undef, population_size)

    gen_col = collect(1:generations)
    mean_col = Vector{Float64}(undef, generations)
    var_col = Vector{Float64}(undef, generations)
    nb_mean_col = Vector{Float64}(undef, generations)
    nb_var_col = Vector{Float64}(undef, generations)
    charles_mean_col = Vector{Float64}(undef, generations)
    sim_p_col = Vector{Float64}(undef, generations)
    app_mean_col = Vector{Float64}(undef, generations)
    app_var_col = Vector{Float64}(undef, generations)
    app_p_col = Vector{Float64}(undef, generations)
    nb_p_col = Vector{Float64}(undef, generations)
    sim_skew_col = Vector{Float64}(undef, generations)
    nb_skew_col = Vector{Float64}(undef, generations)
    sim_exkurt_col = Vector{Float64}(undef, generations)
    nb_exkurt_col = Vector{Float64}(undef, generations)
    sim_m3_col = Vector{Float64}(undef, generations)

    for gen in 1:generations
        population = generate_next_generation!(population, fitness, cum_probs)

        μ_sim, σ²_sim, skew_sim, exkurt_sim, p_sim, m3 = compute_sim_moments(population)
        mean_col[gen] = μ_sim
        var_col[gen] = σ²_sim
        sim_p_col[gen] = p_sim
        sim_skew_col[gen] = skew_sim
        sim_exkurt_col[gen] = exkurt_sim
        sim_m3_col[gen] = m3

        dln_omega_dmu(μ, s_val) = FITNESS_TYPE == "exp" ? (-2 * s_val * μ) : ((-2 * s_val * μ) / max(eps(), (1 - s_val * μ^2)))
        d2ln_omega_dmu2(μ, s_val) = FITNESS_TYPE == "exp" ? (-2 * s_val) : ((-2 * s_val - (2 * s_val^2 * μ^2)) / max(eps(), ((1 - (s_val * μ^2))^2)))

        β₁_app = dln_omega_dmu(μ_app_state, s)
        β₂_app = d2ln_omega_dmu2(μ_app_state, s)
        dμ_dt_app = (u - v) * μ_app_state + σ²_app_state * β₁_app
        dσ²_dt_app = (2u + 0.5) * μ_app_state + (u - v - 0.5) * σ²_app_state +
                     0.5 * (β₁_app + β₂_app * μ_app_state) * σ²_app_state

        μ_app_state = max(0.0, μ_app_state + dμ_dt_app * dt)
        σ²_app_state = max(eps(), σ²_app_state + dσ²_dt_app * dt)
        p_app = μ_app_state / (σ²_app_state + eps())

        app_mean_col[gen] = μ_app_state
        app_var_col[gen] = σ²_app_state
        app_p_col[gen] = p_app

        σ²_app_safe = max(eps(), σ²_app_state)
        p_app_mom = clamp(μ_app_state / σ²_app_safe, eps(), 1.0 - eps())
        ρ_app = (2 - p_app_mom) / (p_app_mom + eps())
        α_app = (6 * (1 - p_app_mom) + p_app_mom^2) / (p_app_mom^2 + eps())
        skew_app = ρ_app / sqrt(σ²_app_safe)
        exkurt_app = α_app / σ²_app_safe

        β₁_nb = dln_omega_dmu(μ_nb_state, s)
        β₂_nb = d2ln_omega_dmu2(μ_nb_state, s)

        σ²_nb_safe = max(eps(), σ²_nb_state)
        p_nb = clamp(μ_nb_state / σ²_nb_safe, eps(), 1.0 - eps())
        ρ_nb = (2 - p_nb) / (p_nb + eps())
        α_nb = (6 * (1 - p_nb) + p_nb^2) / (p_nb^2 + eps())
        γ_nb = ρ_nb / sqrt(σ²_nb_safe)
        Ex_κ_nb = α_nb / σ²_nb_safe

        dμ_dt_nb = (u - v) * μ_nb_state + σ²_nb_state * β₁_nb + (1/2) * ρ_nb * σ²_nb_state * β₂_nb
        dσ²_dt_nb = (2u + 0.5) * μ_nb_state + (u - v - 0.5) * σ²_nb_state +
                    (1/2) * (β₁_nb + β₂_nb * μ_nb_state) * σ²_nb_state +
                    (1/4) * (2 * β₁_nb + β₂_nb) * ρ_nb * σ²_nb_state +
                    (1/4) * α_nb * β₂_nb * σ²_nb_state

        dμ_dt_c = μ_c_state * (u - v) + μ_c_state * dln_omega_dmu(μ_c_state, s)

        μ_nb_state = max(0.0, μ_nb_state + dμ_dt_nb * dt)
        σ²_nb_state = max(eps(), σ²_nb_state + dσ²_dt_nb * dt)
        μ_c_state += dμ_dt_c * dt

        nb_mean_col[gen] = μ_nb_state
        nb_var_col[gen] = σ²_nb_state
        nb_p_col[gen] = p_nb
        nb_skew_col[gen] = γ_nb
        nb_exkurt_col[gen] = Ex_κ_nb
        charles_mean_col[gen] = μ_c_state

        if gen == 1 || (gen % PRINT_EVERY == 0)
            first_chr = length(population[1][1])
            println(
                "Gen $gen: first_chr=$first_chr | " *
                "SIM μ=$(round(μ_sim, digits=3)) σ²=$(round(σ²_sim, digits=3)) skew=$(round(skew_sim, digits=3)) exK=$(round(exkurt_sim, digits=3)) | " *
                "NB  μ=$(round(μ_nb_state, digits=3)) σ²=$(round(σ²_nb_state, digits=3)) skew=$(round(γ_nb, digits=3)) exK=$(round(Ex_κ_nb, digits=3)) | " *
                "APP μ=$(round(μ_app_state, digits=3)) σ²=$(round(σ²_app_state, digits=3)) skew=$(round(skew_app, digits=3)) exK=$(round(exkurt_app, digits=3))"
            )
        end
    end

    if WRITE_CSV
        fitness_suffix = FITNESS_TYPE == "exp" ? "exp" : "syn"
        R_label = FREE_RECOMBINATION ? "free" : R
        csv_file_path = joinpath(
            csv_dir,
            "appNB_s:$(s)_N:$(population_size)_beta:$(u)_delta:$(v)_R:$(R_label)_$(fitness_suffix).csv",
        )

        df = DataFrame(
            generation = gen_col,
            Mean_X_Count = mean_col,
            Approx_Mean = app_mean_col,
            NB_Mean = nb_mean_col,
            Sim_Variance = var_col,
            Approx_Variance = app_var_col,
            NB_Variance = nb_var_col,
            Charles_Mean = charles_mean_col,
            Sim_P = sim_p_col,
            Approx_P = app_p_col,
            NB_P = nb_p_col,
            Sim_Skewness = sim_skew_col,
            NB_Skewness = nb_skew_col,
            Sim_Third_Moment = sim_m3_col,
            Sim_Excess_kurtosis = sim_exkurt_col,
            NB_Excess_Kurtosis = nb_exkurt_col,
        )
        CSV.write(csv_file_path, df)
        println("Wrote: $csv_file_path")
    end

    return nothing
end

@time main()
