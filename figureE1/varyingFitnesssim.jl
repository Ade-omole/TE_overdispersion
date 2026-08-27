# Sweep fitness curvature gamma at fixed u, v and s under free recombination.

using Random, Distributions, CSV, DataFrames, Statistics, Printf

# Configuration
const FITNESS = "exp"
const GAMMAS = [4.0, 2.5]  # run order
const FREE_RECOMBINATION = true

const population_size = 10000
const generations = 50000
const n_init_param = 10 # 10
const s_val = 1.0 / population_size
const dt = 1.0
const v_const = 1.0e-4

const u_values = [0.0005, 0.001, 0.0022, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.026, 0.032]
const v_values = fill(v_const, length(u_values))

const te_root = dirname(@__DIR__)
const csv_dir = joinpath(te_root, "csv_files")
mkpath(csv_dir)
summary_path(gamma) = joinpath(csv_dir, "varyingFitness_sim_$(FITNESS)_gamma$(gamma).csv")

const t0_pw = 9000
const k_pw = 100

transposition_k(u::Float64, n::Int) = n > 0 ? rand(Binomial(n, clamp(u * dt, 0.0, 1.0))) : 0
excision_k(v::Float64, n::Int) = n > 0 ? rand(Binomial(n, clamp(v * dt, 0.0, 1.0))) : 0

fitness_val(n, gamma::Float64) = exp(-s_val * n^gamma)

# Simulation core (as varyingRatesim.jl)
function create_individual()
    n_init = rand(Poisson(n_init_param))
    chrom1, chrom2 = Float64[], Float64[]
    for _ in 1:n_init
        push!(rand(1:2) == 1 ? chrom1 : chrom2, rand())
    end
    return (chrom1, chrom2)
end

function apply_excision!(chrom1::Vector{Float64}, chrom2::Vector{Float64}, v_val::Float64)
    n_total = length(chrom1) + length(chrom2)
    num_del_total = excision_k(v_val, n_total)
    if num_del_total <= 0 || n_total == 0
        return nothing
    end
    n1, n2 = length(chrom1), length(chrom2)
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
            push!(selected, t in selected ? j : t)
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

function apply_transposition!(chrom1::Vector{Float64}, chrom2::Vector{Float64}, u_val::Float64, n_initial::Int)
    k = transposition_k(u_val, n_initial)
    for _ in 1:k
        push!(rand(1:2) == 1 ? chrom1 : chrom2, rand())
    end
    return nothing
end

function meiosis(chrom1::Vector{Float64}, chrom2::Vector{Float64})
    g1, g2 = Float64[], Float64[]
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

function form_offspring(parent1, parent2, u_val::Float64, v_val::Float64)
    g1_p1, g2_p1 = meiosis(parent1[1], parent1[2])
    g1_p2, g2_p2 = meiosis(parent2[1], parent2[2])
    chrom1 = rand(1:2) == 1 ? copy(g1_p1) : copy(g2_p1)
    chrom2 = rand(1:2) == 1 ? copy(g1_p2) : copy(g2_p2)
    n0 = length(chrom1) + length(chrom2)   # pre-transposition/excision count (= m')
    apply_excision!(chrom1, chrom2, v_val)
    apply_transposition!(chrom1, chrom2, u_val, n0)
    return (chrom1, chrom2, n0)
end

function moments(values)
    m = mean(values)
    d = values .- m
    var_ = max(0.0, mean(d .^ 2))
    var_safe = max(eps(), var_)
    skew = mean(d .^ 3) / (var_safe^1.5)
    exkurt = mean(d .^ 4) / (var_safe^2) - 3.0
    return m, var_, skew, exkurt
end

# Run one full simulation for a given (u, v, gamma)
function run_simulation_one_u(gamma::Float64, u_val::Float64, v_val::Float64)
    population = Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef, population_size)
    for i in 1:population_size
        population[i] = create_individual()
    end
    fitness = Vector{Float64}(undef, population_size)
    cum_probs = Vector{Float64}(undef, population_size)
    m0_vals = Vector{Int}(undef, population_size)

    pw_mean, pw_var, pw_skew, pw_exkurt = Float64[], Float64[], Float64[], Float64[]
    pw_mean_m0, pw_var_m0, pw_skew_m0, pw_exkurt_m0 = Float64[], Float64[], Float64[], Float64[]

    for gen in 1:generations
        tot = 0.0
        @inbounds for i in 1:population_size
            p = population[i]
            n = length(p[1]) + length(p[2])
            fi = fitness_val(n, gamma)
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

        new_pop = Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef, population_size)
        @inbounds for i in 1:population_size
            i1 = searchsortedfirst(cum_probs, rand())
            i2 = searchsortedfirst(cum_probs, rand())
            i1 = i1 > population_size ? population_size : i1
            i2 = i2 > population_size ? population_size : i2
            chrom1, chrom2, n0 = form_offspring(population[i1], population[i2], u_val, v_val)
            new_pop[i] = (chrom1, chrom2)
            m0_vals[i] = n0
        end
        population = new_pop

        if gen >= t0_pw && (gen - t0_pw) % k_pw == 0
            counts = [length(p[1]) + length(p[2]) for p in population]
            m, v_, sk, ek = moments(Float64.(counts))
            push!(pw_mean, m); push!(pw_var, v_); push!(pw_skew, sk); push!(pw_exkurt, ek)

            m0, v0, sk0, ek0 = moments(Float64.(m0_vals))
            push!(pw_mean_m0, m0); push!(pw_var_m0, v0); push!(pw_skew_m0, sk0); push!(pw_exkurt_m0, ek0)
        end

        print("\r  gamma=", gamma, " u=", u_val, "  generation ", gen, " / ", generations)
    end
    println()

    sim_mean_post, sim_var_post = mean(pw_mean), mean(pw_var)
    sim_skew_post, sim_exkurt_post = mean(pw_skew), mean(pw_exkurt)
    sim_rho_post = sim_var_post / max(eps(), sim_mean_post)

    sim_mean_pre, sim_var_pre = mean(pw_mean_m0), mean(pw_var_m0)
    sim_skew_pre, sim_exkurt_pre = mean(pw_skew_m0), mean(pw_exkurt_m0)
    sim_rho_pre = sim_var_pre / max(eps(), sim_mean_pre)

    summary_row = (
        u = u_val, v = v_val, gamma = gamma, fit = FITNESS,
        sim_mean_post = sim_mean_post, sim_var_post = sim_var_post, sim_rho_post = sim_rho_post,
        sim_skew_post = sim_skew_post, sim_exkurt_post = sim_exkurt_post,
        sim_mean_pre = sim_mean_pre, sim_var_pre = sim_var_pre, sim_rho_pre = sim_rho_pre,
        sim_skew_pre = sim_skew_pre, sim_exkurt_pre = sim_exkurt_pre,
    )

    @printf("  gamma=%.2f u=%.4f  rho_post=%.4f  rho_pre=%.4f\n", gamma, u_val, sim_rho_post, sim_rho_pre)
    return summary_row
end

# One gamma's full u-sweep -> one CSV, rewritten fully after every u
function run_gamma_sweep(gamma::Float64)
    n_total = length(u_values)
    rows = []
    path = summary_path(gamma)
    for (i, u_val) in enumerate(u_values)
        v_val = v_values[i]
        println("\n[$i/$n_total]  fit=$FITNESS gamma=$gamma  u=$u_val (v=$v_val) ...")
        push!(rows, run_simulation_one_u(gamma, u_val, v_val))
        CSV.write(path, DataFrame(rows))
        println("  ($i/$n_total rows) written to $path")
    end
    println("\nSummary CSV complete: $path")
    return path
end

# Main: run every gamma in GAMMAS, back-to-back, no input needed
println("Running gamma sweep for GAMMAS = $GAMMAS (fitness = $FITNESS), in order, one after another.")
for (gi, gamma) in enumerate(GAMMAS)
    println("\n=========== [$gi/$(length(GAMMAS))] gamma = $gamma ===========")
    try
        run_gamma_sweep(gamma)
    catch e
        @error "gamma=$gamma sweep failed; continuing with the next gamma" exception = (e, catch_backtrace())
    end
end
println("\nAll done.")
