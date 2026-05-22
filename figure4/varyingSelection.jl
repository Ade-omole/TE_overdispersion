# Sweep selection coefficient s in fitness w(n) = exp(-s n^2).
# Roze (2023): beta = 2s (CSV column beta). Keep s_list in sync with varyingSelection.py.

using Random, Distributions, CSV, DataFrames, Statistics, Printf

const FITNESS = "exp"
const EVENT_MODEL = :binomial  # :binomial or :poisson
const R_MAP = 10.0
const u_val = 1e-2
const v_val = 1e-4
const population_size = 10^4
const generations = 50000
const n_init_param = 10
const dt = 1.0
const t0_pw = 9000
const k_pw = 100
const PROGRESS_EVERY = 100  # 0 = no progress prints

const s_list = sort([3e-5, 5e-5, 1e-4, 2e-4, 3e-4, 5e-4, 1e-3, 2e-3])

const te_root = dirname(@__DIR__)
const csv_dir = joinpath(te_root, "csv_files")
const csv_path = joinpath(
    csv_dir,
    "varyingSelection_$(string(EVENT_MODEL))_$(FITNESS)_u$(u_val)_v$(v_val)_4.csv",
)

function transposition_k(n::Int)
    n <= 0 && return 0
    if EVENT_MODEL === :binomial
        return rand(Binomial(n, clamp(u_val * dt, 0.0, 1.0)))
    elseif EVENT_MODEL === :poisson
        return rand(Poisson(u_val * n * dt))
    else
        error("EVENT_MODEL must be :binomial or :poisson")
    end
end

function excision_k(n::Int)
    n <= 0 && return 0
    if EVENT_MODEL === :binomial
        return rand(Binomial(n, clamp(v_val * dt, 0.0, 1.0)))
    elseif EVENT_MODEL === :poisson
        return rand(Poisson(v_val * n * dt))
    else
        error("EVENT_MODEL must be :binomial or :poisson")
    end
end

fitness_val(n::Int, s::Float64, ft::String) = ft == "exp" ? exp(-s * Float64(n)^2) : max(0.0, 1.0 - s * Float64(n)^2)

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
    xovers = sort(rand(n_cross))
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

function form_offspring_roze(parent1, parent2, R_val::Float64)
    g1_p1, g2_p1 = meiosis_roze(parent1[1], parent1[2], R_val)
    g1_p2, g2_p2 = meiosis_roze(parent2[1], parent2[2], R_val)
    chrom1 = rand(1:2) == 1 ? copy(g1_p1) : copy(g2_p1)
    chrom2 = rand(1:2) == 1 ? copy(g1_p2) : copy(g2_p2)
    n0 = length(chrom1) + length(chrom2)
    apply_excision_roze!(chrom1, chrom2)
    apply_transposition_roze!(chrom1, chrom2, n0)
    return (chrom1, chrom2)
end

function run_one(s_val::Float64)
    population = [create_individual_roze() for _ in 1:population_size]
    pop_next = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, population_size)
    fitness = Vector{Float64}(undef, population_size)
    cum_probs = Vector{Float64}(undef, population_size)
    pw_mean = Float64[]
    pw_var = Float64[]
    pw_skew = Float64[]
    pw_exkurt = Float64[]

    for gen in 1:generations
        tot = 0.0
        @inbounds for i in 1:population_size
            p = population[i]
            n = length(p[1]) + length(p[2])
            fi = fitness_val(n, s_val, FITNESS)
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

        @inbounds for i in 1:population_size
            r1, r2 = rand(), rand()
            i1 = min(searchsortedfirst(cum_probs, r1), population_size)
            i2 = min(searchsortedfirst(cum_probs, r2), population_size)
            pop_next[i] = form_offspring_roze(population[i1], population[i2], R_MAP)
        end
        population, pop_next = pop_next, population

        sum_n = 0.0
        sum_n2 = 0.0
        @inbounds for i in 1:population_size
            p = population[i]
            n = length(p[1]) + length(p[2])
            sum_n += n
            sum_n2 += n * n
        end
        μ_sim = sum_n / population_size
        σ²_sim = max(0.0, sum_n2 / population_size - μ_sim * μ_sim)
        σ²_safe = max(eps(), σ²_sim)

        if gen >= t0_pw && (gen - t0_pw) % k_pw == 0
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
            push!(pw_mean, μ_sim)
            push!(pw_var, σ²_sim)
            push!(pw_skew, skew)
            push!(pw_exkurt, exkurt)
        end
        if PROGRESS_EVERY > 0 && gen % PROGRESS_EVERY == 0
            @printf("    gen %7d / %-7d  s=%.4g  μ_sim=%.2f  σ²=%.2f\n", gen, generations, s_val, μ_sim, σ²_sim)
        end
    end
    if PROGRESS_EVERY > 0 && generations % PROGRESS_EVERY != 0
        @printf("    gen %7d / %-7d  s=%.4g  μ_sim=%.2f  σ²=%.2f  (end)\n", generations, generations, s_val, μ_sim, σ²_sim)
    end

    m = mean(pw_mean)
    v = mean(pw_var)
    rho = v / max(eps(), m)
    sk = mean(pw_skew)
    ek = mean(pw_exkurt)
    return (sim_μ = m, sim_σ² = v, rho_sim = rho, sim_skew = sk, sim_exkurt = ek)
end

function main()
    mkpath(csv_dir)
    rows = []
    n = length(s_list)
    for i in 1:n
        s = s_list[i]
        β = 2.0 * s
        println("\n[$i/$n]  s = $s  (Roze β = 2s = $β)  R=$R_MAP  u=$u_val  v=$v_val")
        r = run_one(s)
        push!(rows, (
            beta = β,
            s = s,
            R = R_MAP,
            u = u_val,
            v = v_val,
            sim_μ = r.sim_μ,
            sim_σ² = r.sim_σ²,
            rho_sim = r.rho_sim,
            sim_skew = r.sim_skew,
            sim_exkurt = r.sim_exkurt,
        ))
    end
    df = DataFrame(rows)
    CSV.write(csv_path, df)
    println("\nWrote ", csv_path)
end

main()
