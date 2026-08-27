# Sweep map length R under the Roze (2023) ectopic-pair fitness w = exp(-beta * n_ecto).

using Random, Distributions, CSV, DataFrames, Statistics, Printf, QuadGK

# Configuration
const FITNESS = "exp" # ectopic roze
const THEORY_FITNESS = "exp"   # for μ_app, σ²_app, and NB ODE: ln w = -s n²
const u_val = 10^-2
const v_val = u_val / 100
const α_for_rho = 0.0   # α in ρ formula: ρ ≈ 1 + (ε₁/(1-uε₁)) * ((u+v+α)/2). Set 0 or change as needed.

const population_size = 10^4
const num_diploid_pairs = 1
const generations = 30000
const n_init_param = 10
# Simulation: β in w = exp(-β · n_ecto) (α = 0).
const β_ect = u_val / 10
# Theory (approx + NB) uses s in w = exp(-s n^2), independent of beta_ect
const s_val = u_val / 20
const dt = 1.0

# u/v event model (switch with one constant)
const UV_EVENT_MODEL = :binomial  # :poisson or :binomial

if UV_EVENT_MODEL == :poisson
    transposition_k(n::Int) = n > 0 ? rand(Poisson(u_val * n * dt)) : 0
    excision_k(n::Int) = n > 0 ? rand(Poisson(v_val * n * dt)) : 0
elseif UV_EVENT_MODEL == :binomial
    transposition_k(n::Int) = n > 0 ? rand(Binomial(n, clamp(u_val * dt, 0.0, 1.0))) : 0
    excision_k(n::Int) = n > 0 ? rand(Binomial(n, clamp(v_val * dt, 0.0, 1.0))) : 0
else
    error("UV_EVENT_MODEL must be :poisson or :binomial (got $UV_EVENT_MODEL)")
end

const R_values = [10.0^k for k in -5:3]

const te_root = dirname(@__DIR__)
const csv_dir = joinpath(te_root, "csv_files")
const CSV_BASENAME = "sweep_R_ectFit_burnin_binomial_pre"
mkpath(csv_dir)s

const t0_pw = 9000
const k_pw = 100

# no burn-in

# Roze-style simulation
function create_individual_roze()
    n_init = rand(Poisson(n_init_param))
    chrom1, chrom2 = Float64[], Float64[]
    for _ in 1:n_init
        push!(rand(1:2) == 1 ? chrom1 : chrom2, rand())
    end
    return (chrom1, chrom2)
end

function apply_excision_roze!(chrom1, chrom2)
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

function apply_transposition_roze!(chrom1, chrom2, n_initial::Int)
    k = transposition_k(n_initial)
    for _ in 1:k
        push!(rand(1:2) == 1 ? chrom1 : chrom2, rand())
    end
    return nothing
end

function meiosis_roze(chrom1, chrom2, R_val::Float64)
    n_cross = rand(Poisson(R_val))
    if n_cross == 0
        return (chrom1, chrom2)
    end
    xovers = n_cross > 0 ? sort(rand(n_cross)) : Float64[]
    g1, g2 = Float64[], Float64[]
    for p in chrom1
        seg = searchsortedfirst(xovers, p)
        isodd(seg) ? push!(g1, p) : push!(g2, p)
    end
    for p in chrom2
        seg = searchsortedfirst(xovers, p)
        iseven(seg) ? push!(g1, p) : push!(g2, p)
    end
    return (g1, g2)
end

"""Count insertion sites shared by both homologs (exact float equality), after sorting each chromosome."""
function count_homozygous_te_sites(chrom1::Vector{Float64}, chrom2::Vector{Float64})
    (isempty(chrom1) || isempty(chrom2)) && return 0
    s1 = sort(chrom1)
    s2 = sort(chrom2)
    n_hom = 0
    a, b = 1, 1
    n1, n2 = length(s1), length(s2)
    @inbounds while a <= n1 && b <= n2
        if s1[a] < s2[b]
            a += 1
        elseif s1[a] > s2[b]
            b += 1
        else
            n_hom += 1
            a += 1
            b += 1
        end
    end
    return n_hom
end

"""Roze fitness with α = 0: w = exp(-β · n_ecto), n_ecto = C(n,2) - n_hom."""
function fitness_ectopic_roze(chrom1::Vector{Float64}, chrom2::Vector{Float64}, β::Float64)
    n = length(chrom1) + length(chrom2)
    n < 2 && return 1.0
    n_hom = count_homozygous_te_sites(chrom1, chrom2)
    n_ecto = div(n * (n - 1), 2) - n_hom
    return exp(-β * n_ecto)
end

function form_offspring_roze(parent1, parent2, R_val::Float64)
    g1_p1, g2_p1 = meiosis_roze(parent1[1], parent1[2], R_val)
    g1_p2, g2_p2 = meiosis_roze(parent2[1], parent2[2], R_val)
    chrom1 = rand(1:2) == 1 ? copy(g1_p1) : copy(g2_p1)
    chrom2 = rand(1:2) == 1 ? copy(g1_p2) : copy(g2_p2)
    n0 = length(chrom1) + length(chrom2)   # pre-transposition/excision count (= m')
    apply_excision_roze!(chrom1, chrom2)
    apply_transposition_roze!(chrom1, chrom2, n0)
    return (chrom1, chrom2, n0)
end

function dln_omega_dmu(μ, s, ft::String)
    ft == "exp" ? -2.0 * s * μ : ((-2.0 * s * μ) / max(eps(), 1.0 - s * μ^2))
end

function d2ln_omega_dmu2(μ, s, ft::String)
    ft == "exp" ? -2.0 * s : ((-2.0 * s - 2.0 * s^2 * μ^2) / max(eps(), (1.0 - s * μ^2)^2))
end

## E₁ (ε₁) for linear genetic map of length R Morgans (Roze appendix)
function compute_E1(R_val::Float64, u_par::Float64)
    ρ_R = u_par / max(R_val, 1e-12)
    if R_val <= 1.0
        half_plus_rho = 0.5 + ρ_R
        if ρ_R <= 0.0 || half_plus_rho <= 0.0
            return 0.0
        end
        return (2.0 / R_val) * ((1.0 + 2.0 * ρ_R) * (log(half_plus_rho) - log(max(ρ_R, 1e-12))) - 1.0)
    else
        integrand(x) = (R_val - x) / (max(eps(), (1.0 - exp(-2.0 * x)) / 2.0 + 2.0 * u_par))
        res, _ = quadgk(integrand, 0.0, R_val; order=7)
        return (2.0 / (R_val^2)) * res
    end
end

function compute_rho_numerical(E1::Float64, u_par::Float64, v_par::Float64, α_par::Float64)
    denom = 1.0 - u_par * E1
    denom = abs(denom) < 1e-10 ? sign(denom) * 1e-10 : denom
    return 1.0 + (E1 / denom) * ((u_par + v_par + α_par) / 2.0)
end

# Run one full simulation for a given R
function run_simulation_one_R(R_val::Float64)
    # Approx mean/variance for theory fitness w = exp(-s n²) (not ectopic).
    denom_app = 2.0 * s_val * (2.0 * u_val + 2.0 * v_val + 1.0)
    μ_app = denom_app > eps() ? ((u_val - v_val) * (1.0 - 2.0 * (u_val - v_val))) / denom_app : 0.0
    μ_app = max(0.0, μ_app)
    σ²_app = max(eps(), (u_val - v_val) / (2.0 * s_val))
    μ_nb = σ²_nb = 10.0
    population = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, population_size)
    for i in 1:population_size
        population[i] = create_individual_roze()
    end

    fitness = Vector{Float64}(undef, population_size)
    cum_probs = Vector{Float64}(undef, population_size)
    m0_vals = Vector{Int}(undef, population_size)
    pw_mean = Float64[]
    pw_var = Float64[]
    pw_skew = Float64[]
    pw_exkurt = Float64[]

    # Pre-transposition/excision (m') moments, recorded at the same generations
    pw_mean_m0 = Float64[]
    pw_var_m0 = Float64[]
    pw_skew_m0 = Float64[]
    pw_exkurt_m0 = Float64[]

    for gen in 1:generations
        tot = 0.0
        @inbounds for i in 1:population_size
            p = population[i]
            fi = fitness_ectopic_roze(p[1], p[2], β_ect)
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
            r1 = rand()
            r2 = rand()
            i1 = searchsortedfirst(cum_probs, r1)
            i2 = searchsortedfirst(cum_probs, r2)
            i1 = i1 > population_size ? population_size : i1
            i2 = i2 > population_size ? population_size : i2
            chrom1, chrom2, n0 = form_offspring_roze(population[i1], population[i2], R_val)
            new_pop[i] = (chrom1, chrom2)
            m0_vals[i] = n0
        end
        population = new_pop

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

        # Same moment computation, applied to the pre-transposition/excision counts (m')
        sum_n_m0 = 0.0
        @inbounds for i in 1:population_size
            sum_n_m0 += m0_vals[i]
        end
        μ_sim_m0 = sum_n_m0 / population_size

        sum_sq_m0 = 0.0
        @inbounds for i in 1:population_size
            d = m0_vals[i] - μ_sim_m0
            sum_sq_m0 += d * d
        end
        σ²_sim_m0 = max(0.0, sum_sq_m0 / population_size)
        σ²_safe_m0 = max(eps(), σ²_sim_m0)

        sum_c3_m0 = 0.0
        sum_c4_m0 = 0.0
        @inbounds for i in 1:population_size
            d = m0_vals[i] - μ_sim_m0
            d2 = d * d
            sum_c3_m0 += d2 * d
            sum_c4_m0 += d2 * d2
        end
        m3_m0 = sum_c3_m0 / population_size
        m4_m0 = sum_c4_m0 / population_size
        skew_m0 = m3_m0 / (σ²_safe_m0^1.5)
        exkurt_m0 = (m4_m0 / (σ²_safe_m0^2)) - 3.0

        if gen >= t0_pw && (gen - t0_pw) % k_pw == 0
            push!(pw_mean, μ_sim)
            push!(pw_var, σ²_sim)
            push!(pw_skew, skew)
            push!(pw_exkurt, exkurt)

            push!(pw_mean_m0, μ_sim_m0)
            push!(pw_var_m0, σ²_sim_m0)
            push!(pw_skew_m0, skew_m0)
            push!(pw_exkurt_m0, exkurt_m0)
        end

        β₁ = dln_omega_dmu(μ_nb, s_val, THEORY_FITNESS)
        β₂ = d2ln_omega_dmu2(μ_nb, s_val, THEORY_FITNESS)
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

    sim_mean_pw_m0 = mean(pw_mean_m0)
    sim_var_pw_m0 = mean(pw_var_m0)
    sim_skew_pw_m0 = mean(pw_skew_m0)
    sim_exkurt_pw_m0 = mean(pw_exkurt_m0)
    sim_p_pw_m0 = sim_mean_pw_m0 / max(eps(), sim_var_pw_m0)

    σ²_nb_last = max(eps(), σ²_nb)
    p_nb_last = clamp(μ_nb / σ²_nb_last, eps(), 1.0 - eps())
    ρ_last = (2.0 - p_nb_last) / p_nb_last
    α_last = (6.0 * (1.0 - p_nb_last) + p_nb_last^2) / max(eps(), p_nb_last^2)
    nb_skew_last = ρ_last / sqrt(σ²_nb_last)
    nb_exkurt_last = α_last / σ²_nb_last

    p_app_last = μ_app / max(eps(), σ²_app)

    E1_val = compute_E1(R_val, u_val)
    rho_num = compute_rho_numerical(E1_val, u_val, v_val, α_for_rho)

    inv_p_sim = 1.0 / max(eps(), sim_p_pw)
    inv_p_app = 1.0 / max(eps(), p_app_last)
    inv_p_nb = 1.0 / max(eps(), p_nb_last)
    inv_p_sim_m0 = 1.0 / max(eps(), sim_p_pw_m0)   # rho at the pre-transposition/excision stage

    summary_row = (
        R = R_val,
        u = u_val,
        v = v_val,
        beta_ect = β_ect,
        s_theory = s_val,
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
        E1 = E1_val,
        rho_numerical = rho_num,
        inv_p_sim = inv_p_sim,
        inv_p_app = inv_p_app,
        inv_p_nb = inv_p_nb,
        # Pre-transposition/excision (m') simulation moments
        sim_μ_m0 = sim_mean_pw_m0,
        sim_σ²_m0 = sim_var_pw_m0,
        sim_p_m0 = sim_p_pw_m0,
        sim_skew_m0 = sim_skew_pw_m0,
        sim_exkurt_m0 = sim_exkurt_pw_m0,
        inv_p_sim_m0 = inv_p_sim_m0,
    )

    println("  ┌─────────────┬──────────────┬──────────────┬──────────────┬──────────────┬──────────────┐")
    println("  │             │     Mean     │   Variance   │      p       │   Skewness   │  Ex. Kurtosis│")
    println("  ├─────────────┼──────────────┼──────────────┼──────────────┼──────────────┼──────────────┤")
    @printf("  │ Simulation  │ %12.2f │ %12.2f │ %12.4f │ %12.2f │ %12.2f │\n",
            sim_mean_pw, sim_var_pw, sim_p_pw, sim_skew_pw, sim_exkurt_pw)
    @printf("  │ Pre-transp. │ %12.2f │ %12.2f │ %12.4f │ %12.2f │ %12.2f │\n",
            sim_mean_pw_m0, sim_var_pw_m0, sim_p_pw_m0, sim_skew_pw_m0, sim_exkurt_pw_m0)
    @printf("  │ Approx      │ %12.2f │ %12.2f │ %12.4f │      —       │      —       │\n",
            μ_app, σ²_app, p_app_last)
    @printf("  │ Theory (NB) │ %12.2f │ %12.2f │ %12.4f │ %12.2f │ %12.2f │\n",
            μ_nb, σ²_nb, p_nb_last, nb_skew_last, nb_exkurt_last)
    println("  └─────────────┴──────────────┴──────────────┴──────────────┴──────────────┴──────────────┘")
    println("  E₁ = $E1_val,  ρ (numerical) = $rho_num")
    println("  1/p: sim = $inv_p_sim,  app = $inv_p_app,  nb = $inv_p_nb  (approx/NB: exp(-s·n²), s = $s_val)")
    println("  1/p (pre-transposition/excision): sim = $inv_p_sim_m0")

    return summary_row
end

# Main sweep
summary_path = joinpath(csv_dir, "$(CSV_BASENAME)_$(FITNESS).csv")
n_total = length(R_values)
for (i, R_val) in enumerate(R_values)
    println("\n[$i/$n_total]  Running R = $R_val  (u = $u_val, v = $v_val, sim β_ect = $β_ect, theory s = $s_val, rate model = $UV_EVENT_MODEL)")
    summary_row = run_simulation_one_R(R_val)
    df_row = DataFrame([summary_row])
    if i == 1
        CSV.write(summary_path, df_row)
    else
        CSV.write(summary_path, df_row; append=true)
    end
end
println("\nSummary CSV written: $summary_path")
