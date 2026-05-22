# Sweep transposition rate u with fixed v.
# Writes CSV under csv_files/.

using Random, Distributions, CSV, DataFrames, Statistics, Printf

# Configuration
const FITNESS = "exp"
const FREE_RECOMBINATION = true  # true: free recombination; false: map-length R

const population_size = 10000
const num_diploid_pairs = 1
const generations = 50000
const n_init_param = 10
const R = 10.0
const s_val = 1.0 / population_size
const dt = 1.0

const UV_EVENT_MODEL = :binomial  # :poisson or :binomial

if UV_EVENT_MODEL == :poisson
    transposition_k(u::Float64, n::Int) = n > 0 ? rand(Poisson(u * n * dt)) : 0
    excision_k(v::Float64, n::Int) = n > 0 ? rand(Poisson(v * n * dt)) : 0
elseif UV_EVENT_MODEL == :binomial
    # Each copy experiences an event independently with probability p ≈ rate*dt.
    transposition_k(u::Float64, n::Int) = n > 0 ? rand(Binomial(n, clamp(u * dt, 0.0, 1.0))) : 0
    excision_k(v::Float64, n::Int) = n > 0 ? rand(Binomial(n, clamp(v * dt, 0.0, 1.0))) : 0
else
    error("UV_EVENT_MODEL must be :poisson or :binomial (got $UV_EVENT_MODEL)")
end


const v_const = 0.0001
const u_values = [0.0005, 0.001, 0.0015, 0.0022, 0.003, 0.004, 0.005, 0.006, 0.008, 0.01, 0.012, 0.016, 0.02, 0.026, 0.032, 0.04, 0.05, 0.07, 0.08]
const v_values = fill(v_const, length(u_values))

const te_root = dirname(@__DIR__)
const csv_dir = joinpath(te_root, "csv_files")
const CSV_BASENAME = "sweep_u_roze_more_values_vConst:0.0001_5_$UV_EVENT_MODEL"
mkpath(csv_dir)

const t0_pw = 9000
const k_pw = 100

# Simulation

function create_individual_roze()
    n_init = rand(Poisson(n_init_param))
    chrom1, chrom2 = Float64[], Float64[]
    for _ in 1:n_init
        push!(rand(1:2) == 1 ? chrom1 : chrom2, rand())
    end
    return (chrom1, chrom2)
end

# Faster excision: no shuffle(1:n) and no repeated deleteat!
function apply_excision_roze!(chrom1::Vector{Float64}, chrom2::Vector{Float64}, v_val::Float64, dt_step::Float64)
    n_total = length(chrom1) + length(chrom2)
    num_del_total = excision_k(v_val, n_total)
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

function apply_transposition_roze!(chrom1::Vector{Float64}, chrom2::Vector{Float64}, u_val::Float64, dt_step::Float64, n_initial::Int)
    k = transposition_k(u_val, n_initial)
    for _ in 1:k
        push!(rand(1:2) == 1 ? chrom1 : chrom2, rand())
    end
    return nothing
end

function meiosis_roze(chrom1::Vector{Float64}, chrom2::Vector{Float64})
    if FREE_RECOMBINATION
        # Independent assortment per TE, without concat+shuffle (faster, exact).
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
    g1, g2 = Float64[], Float64[]
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

function fitness_val(x_total, fit_type::String)
    fit_type == "exp" ? exp(-s_val * x_total^2) : max(0.0, 1.0 - s_val * x_total^2)
end

function form_offspring_roze(parent1, parent2, u_val::Float64, v_val::Float64)
    g1_p1, g2_p1 = meiosis_roze(parent1[1], parent1[2])
    g1_p2, g2_p2 = meiosis_roze(parent2[1], parent2[2])
    chrom1 = rand(1:2) == 1 ? copy(g1_p1) : copy(g2_p1)
    chrom2 = rand(1:2) == 1 ? copy(g1_p2) : copy(g2_p2)
    n0 = length(chrom1) + length(chrom2)
    apply_excision_roze!(chrom1, chrom2, v_val, dt)
    apply_transposition_roze!(chrom1, chrom2, u_val, dt, n0)
    return (chrom1, chrom2)
end

function dln_omega_dmu(μ, s, ft::String)
    ft == "exp" ? -2.0 * s * μ : ((-2.0 * s * μ) / max(eps(), 1.0 - s * μ^2))
end

function d2ln_omega_dmu2(μ, s, ft::String)
    ft == "exp" ? -2.0 * s : ((-2.0 * s - 2.0 * s^2 * μ^2) / max(eps(), (1.0 - s * μ^2)^2))
end

function run_simulation_one_u(u_val::Float64, v_val::Float64)
    denom_app = 2.0 * s_val * (2.0 * u_val + 2.0 * v_val + 1.0)
    μ_app = denom_app > eps() ? ((u_val - v_val) * (1.0 - 2.0 * (u_val - v_val))) / denom_app : 0.0
    μ_app = max(0.0, μ_app)
    σ²_app = max(eps(), (u_val - v_val) / (2.0 * s_val))
    μ_nb = σ²_nb = 10.0

    # Typed population + reusable buffers
    population = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, population_size)
    for i in 1:population_size
        population[i] = create_individual_roze()
    end
    fitness = Vector{Float64}(undef, population_size)
    cum_probs = Vector{Float64}(undef, population_size)

    pw_mean = Float64[]
    pw_var = Float64[]
    pw_skew = Float64[]
    pw_exkurt = Float64[]

    for gen in 1:generations
        # Parent sampling
        tot = 0.0
        @inbounds for i in 1:population_size
            p = population[i]
            n = length(p[1]) + length(p[2])
            fi = fitness_val(n, FITNESS)
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

        new_pop = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, population_size)
        @inbounds for i in 1:population_size
            i1 = searchsortedfirst(cum_probs, rand())
            i2 = searchsortedfirst(cum_probs, rand())
            i1 = i1 > population_size ? population_size : i1
            i2 = i2 > population_size ? population_size : i2
            new_pop[i] = form_offspring_roze(population[i1], population[i2], u_val, v_val)
        end
        population = new_pop

        # Moments without allocating counts vector
        sum_n = 0.0
        @inbounds for i in 1:population_size
            p = population[i]
            sum_n += (length(p[1]) + length(p[2]))
        end
        μ_sim = sum_n / population_size

        sum_sq = 0.0
        @inbounds for i in 1:population_size
            p = population[i]
            d = (length(p[1]) + length(p[2])) - μ_sim
            sum_sq += d * d
        end
        σ²_sim = max(0.0, sum_sq / population_size)
        σ²_safe = max(eps(), σ²_sim)

        sum_c3 = 0.0
        sum_c4 = 0.0
        @inbounds for i in 1:population_size
            p = population[i]
            d = (length(p[1]) + length(p[2])) - μ_sim
            d2 = d * d
            sum_c3 += d2 * d
            sum_c4 += d2 * d2
        end
        m3 = sum_c3 / population_size
        m4 = sum_c4 / population_size
        skew = m3 / (σ²_safe^1.5)
        exkurt = (m4 / (σ²_safe^2)) - 3.0

        if gen >= t0_pw && (gen - t0_pw) % k_pw == 0
            push!(pw_mean, μ_sim)
            push!(pw_var, σ²_sim)
            push!(pw_skew, skew)
            push!(pw_exkurt, exkurt)
        end

        # NB ODE update
        β₁ = dln_omega_dmu(μ_nb, s_val, FITNESS)
        β₂ = d2ln_omega_dmu2(μ_nb, s_val, FITNESS)
        p_nb = clamp(μ_nb / max(eps(), σ²_nb), eps(), 1.0 - eps())
        ρ = (2.0 - p_nb) / p_nb
        α = (6.0 * (1.0 - p_nb) + p_nb^2) / max(eps(), p_nb^2)
        dμ = (u_val - v_val) * μ_nb + σ²_nb * β₁ + 0.5 * ρ * σ²_nb * β₂
        dσ² = (2u_val + 0.5) * μ_nb + (u_val - v_val - 0.5) * σ²_nb +
              0.5 * (β₁ + β₂ * μ_nb) * σ²_nb +
              0.25 * (2β₁ + β₂) * ρ * σ²_nb + 0.25 * α * β₂ * σ²_nb
        μ_nb = max(0.0, μ_nb + dμ * dt)
        σ²_nb = max(eps(), σ²_nb + dσ² * dt)

        print("\r  generation ", gen, " / ", generations)
    end
    println()

    sim_mean_pw = mean(pw_mean)
    sim_var_pw = mean(pw_var)
    sim_skew_pw = mean(pw_skew)
    sim_exkurt_pw = mean(pw_exkurt)
    sim_p_pw = sim_mean_pw / max(eps(), sim_var_pw)

    σ²_nb_last = max(eps(), σ²_nb)
    p_nb_last = clamp(μ_nb / σ²_nb_last, eps(), 1.0 - eps())
    ρ_last = (2.0 - p_nb_last) / p_nb_last
    α_last = (6.0 * (1.0 - p_nb_last) + p_nb_last^2) / max(eps(), p_nb_last^2)
    nb_skew_last = ρ_last / sqrt(σ²_nb_last)
    nb_exkurt_last = α_last / σ²_nb_last

    p_app_last = μ_app / max(eps(), σ²_app)

    summary_row = (
        u = u_val,
        v = v_val,
        u_minus_v = u_val - v_val,
        fit = FITNESS,
        sim_μ = sim_mean_pw,
        app_μ = μ_app,
        nb_μ = μ_nb,
        sim_σ² = sim_var_pw,
        app_σ² = σ²_app,
        nb_σ² = σ²_nb,
        sim_p = sim_p_pw,
        app_p = p_app_last,
        nb_p = p_nb_last,
        sim_skew = sim_skew_pw,
        nb_skew = nb_skew_last,
        sim_exkurt = sim_exkurt_pw,
        nb_exkurt = nb_exkurt_last,
    )

    println("  ┌─────────────┬──────────────┬──────────────┬──────────────┬──────────────┬──────────────┐")
    println("  │             │     Mean     │   Variance   │      p       │   Skewness   │  Ex. Kurtosis│")
    println("  ├─────────────┼──────────────┼──────────────┼──────────────┼──────────────┼──────────────┤")
    @printf("  │ Simulation  │ %12.2f │ %12.2f │ %12.4f │ %12.2f │ %12.2f │\n",
            sim_mean_pw, sim_var_pw, sim_p_pw, sim_skew_pw, sim_exkurt_pw)
    @printf("  │ Approx      │ %12.2f │ %12.2f │ %12.4f │      —       │      —       │\n",
            μ_app, σ²_app, p_app_last)
    @printf("  │ Theory (NB) │ %12.2f │ %12.2f │ %12.4f │ %12.2f │ %12.2f │\n",
            μ_nb, σ²_nb, p_nb_last, nb_skew_last, nb_exkurt_last)
    println("  └─────────────┴──────────────┴──────────────┴──────────────┴──────────────┴──────────────┘")

    return summary_row
end

summary_path = joinpath(csv_dir, "$(CSV_BASENAME)_$(FITNESS).csv")
n_total = length(u_values)
for (i, u_val) in enumerate(u_values)
    v_val = v_values[i]
    println("\n[$i/$n_total]  Running u = $u_val  (v = $v_val, u-v = $(u_val - v_val), fitness = $FITNESS, u/v model = $UV_EVENT_MODEL) ...")
    summary_row = run_simulation_one_u(u_val, v_val)
    df_row = DataFrame([summary_row])
    if i == 1
        CSV.write(summary_path, df_row)
    else
        CSV.write(summary_path, df_row; append=true)
    end
end
println("\nSummary CSV written: $summary_path")

