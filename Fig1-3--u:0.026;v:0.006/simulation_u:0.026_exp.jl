using Pkg
Pkg.add("Random")
Pkg.add("Distributions")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("Statistics")

using Random, Distributions, CSV, DataFrames, Statistics

######################################################### Parameters #######################################################

L = 50000  ## Length of each haploid chromosome (loci)
num_TEs = 5  ## Initial number of autonomous TEs ("X") per chromosome
population_size = 10000  ## Total population size (number of individuals)
num_diploid_pairs = 1  ## Each individual has 1 pair of diploid chromosomes (2 haploid chromosomes total)
generations = 20000  ## Number of generations to simulate

s = 1 / population_size ## Selection coefficient against each TE element
r = 1  ## Recombination rate (Poisson mean); r = 1(Poisson mean) <==> r = 0.5cM/Mb (recombination rate) i.e., free recombination 

u = 0.026
v = 0.006


##### Initial conditions
dt = 1.0                 ## Time step
μ_c = μ_app = μ_nb = σ²_app = σ²_nb = 9.8   ## Initial mean and variance of TE copy number per individual

## chromosome_length = 900
# u = -(1/4) * d2ln_omega_dmu2(μ_f, s) *  μ_f   
# v = 0.0004 

################################################ Simulations ###################################################

######## Function to create a single haploid chromosome with autonomous TEs
function create_haploid_chromosome()
    chromosome = fill('.', L)  ## Initialize chromosome with empty sites represented by '.'
    autonomous_sites = rand(1:L, num_TEs)  ## Randomly select sites for TEs ('X')
    for site in autonomous_sites
        chromosome[site] = 'X'  ## Place 'X' at the selected sites
    end
    return chromosome
end

######## Function to create an individual with 1 pair of diploid chromosomes
function create_individual()
    chromosomes = [create_haploid_chromosome() for _ in 1:(num_diploid_pairs * 2)]  ## Create two haploid chromosomes for the individual
    return chromosomes
end

######## Function to create the initial population
function create_population(population_size)
    population = [create_individual() for _ in 1:population_size]  ## Generate a population of individuals
    return population
end

######## Function to perform crossover between homologous haploid chromosomes
function crossover(chromosome1, chromosome2)
    num_crossovers = rand(Poisson(r))  ## Randomly determine the number of crossovers
    if num_crossovers == 0
        return chromosome1, chromosome2  ## Return chromosomes unchanged if no crossovers occur
    end

    crossover_positions = sort(rand(1:L, num_crossovers))  ## Get crossover positions in sorted order
    for pos in crossover_positions
        if pos < L
            ## Swap segments after each crossover point
            chromosome1[pos:end], chromosome2[pos:end] = chromosome2[pos:end], chromosome1[pos:end]
        end
    end
    return chromosome1, chromosome2
end

######## Function to update both chromosomes of an individual by duplicating and deleting TEs based on propensity
function update_chromosomes(chromosome1, chromosome2)
    ## Count the initial number of 'X' elements in both chromosomes
    x_old1 = count(c -> c == 'X', chromosome1)
    x_old2 = count(c -> c == 'X', chromosome2)
    x_old_total = x_old1 + x_old2

    ## Create copies of the chromosomes
    new_chromosome1 = copy(chromosome1)
    new_chromosome2 = copy(chromosome2)

    ## Calculate the number of 'X' to be deleted from each chromosome
    num_x_to_delete_total = rand(Poisson(v * x_old_total * dt))
    
    ## Distribute deletions proportionally between chromosomes
    if x_old_total > 0
        num_x_to_delete1 = rand(Binomial(num_x_to_delete_total, x_old1/x_old_total))
        num_x_to_delete2 = num_x_to_delete_total - num_x_to_delete1
    else
        num_x_to_delete1 = 0
        num_x_to_delete2 = 0
    end

    ## Function to delete elements from a chromosome
    function delete_elements!(chromosome, num_to_delete)
        positions = findall(c -> c == 'X', chromosome)
        if !isempty(positions) && num_to_delete > 0
            delete_positions = rand(positions, min(length(positions), num_to_delete))
            for pos in delete_positions
                chromosome[pos] = '.'
            end
        end
    end
    
    ######## Perform random deletion of 'X' on both chromosomes
    delete_elements!(new_chromosome1, num_x_to_delete1)
    delete_elements!(new_chromosome2, num_x_to_delete2)

    ######## Replication Process
    ######## Number of elements to add in total
    num_x_to_add_total = rand(Poisson(u * x_old_total * dt))
    
    ######## Find all available positions in both chromosomes
    available_positions1 = findall(c -> c == '.', new_chromosome1)
    available_positions2 = findall(c -> c == '.', new_chromosome2)
    all_available_positions = vcat(
        [(1, pos) for pos in available_positions1],
        [(2, pos) for pos in available_positions2]
    )
    
    ######## Place new elements randomly in available positions
    for _ in 1:num_x_to_add_total
        if !isempty(all_available_positions)
            idx = rand(1:length(all_available_positions))
            (chrom_num, pos) = all_available_positions[idx]
            if chrom_num == 1
                new_chromosome1[pos] = 'X'
            else
                new_chromosome2[pos] = 'X'
            end
            deleteat!(all_available_positions, idx)  # Remove this position from available
        end
    end

    return new_chromosome1, new_chromosome2
end

######## Function to calculate fitness and parent probabilities
function calculate_parent_probabilities(population)

    ######## Calculate fitness for each individual
    fitness = []
    for individual in population
        x_total = sum(count(c -> c == 'X', chromosome) for chromosome in individual) 
        # println("x_total: $x_total")    
        # fitness_value =  1 - (s * x_total^2) # Synergistic fitness function
        fitness_value = exp(-s * x_total^2)  # Exponential fitness function
        # fitness_value = 1 - (s * x_total) # Linear fitness function
        # fitness_value = 1 - (s * x_total^1.5)  # Charlesworth 1983 fitness function
        
        if fitness_value < 0
            println("fitness less than 0")
        end 

        # println("fitness_value: $fitness_value")

        push!(fitness, max(fitness_value, 0))  # Ensure fitness is non-negative
    end
    # print(fitness)
    ######## Calculate probabilities for each individual
    total_fitness = sum(fitness)
    probabilities = [f / total_fitness for f in fitness]

    ######## Calculate cumulative probabilities
    cumulative_probabilities = cumsum(probabilities)
    
    return cumulative_probabilities
end

######## Function to select parents based on cumulative probabilities
function select_parents(cumulative_probabilities)
    r1, r2 = rand(), rand()  # Generate two random numbers uniformly between 0 and 1
    parent1_index = findfirst(p -> p >= r1, cumulative_probabilities)
    parent2_index = findfirst(p -> p >= r2, cumulative_probabilities)
    # Handle the case where findfirst returns nothing
    if parent1_index === nothing
        parent1_index = length(cumulative_probabilities)
    end
    if parent2_index === nothing
        parent2_index = length(cumulative_probabilities)
    end
    return parent1_index, parent2_index
end

######## Function to perform recombination within a diploid pair of chromosomes
function recombine_diploid_pair(diploid_pair)
    chrom1, chrom2 = diploid_pair
    recombined_chrom1, recombined_chrom2 = crossover(copy(chrom1), copy(chrom2))  ## Perform crossover and return recombined chromosomes
    return recombined_chrom1, recombined_chrom2
end

######## Function to count the number of X elements in each chromosome of the offspring
function count_TEs_per_chromosome(offspring)
    X_counts = []
    for chromosome in offspring
        push!(X_counts, count(c -> c == 'X', chromosome))  ## Count 'X' elements
    end
    return X_counts
end

######## Function to generate the next generation by performing recombination and update
function generate_next_generation(gen)
    global population, μ_c, μ_app, μ_nb, m, σ²_app, σ²_nb # Use to track the mean/variance copies 

    new_population = []

    generation_X_count = 0  ## Track total X count for the generation

    ######## Create a DataFrame to store individual chromosome counts for the final generation
    individual_chromosome_counts = DataFrame(Individual_ID = Int[], X_1 = Int[], X_2 = Int[])

    ######## Calculate cumulative probabilities for parent selection
    cumulative_probabilities = calculate_parent_probabilities(population)

    for i in 1:population_size
        ## Select two parents based on fitness-weighted probabilities
        parent_indices = select_parents(cumulative_probabilities)
        parent1, parent2 = population[parent_indices[1]], population[parent_indices[2]]

        parent1_recombined = []
        parent2_recombined = []

        crossover_occurred = false  ## Track if crossover occurred

        ######## Recombine chromosomes from both parents
        recombined_chrom1, recombined_chrom2 = recombine_diploid_pair((parent1[1], parent1[2]))
        push!(parent1_recombined, recombined_chrom1)
        push!(parent1_recombined, recombined_chrom2)

        recombined_chrom3, recombined_chrom4 = recombine_diploid_pair((parent2[1], parent2[2]))
        push!(parent2_recombined, recombined_chrom3)
        push!(parent2_recombined, recombined_chrom4)

        ######## Check if crossover occurred
        if recombined_chrom1 != parent1[1] || recombined_chrom2 != parent1[2]
            crossover_occurred = true
        end
        if recombined_chrom3 != parent2[1] || recombined_chrom4 != parent2[2]
            crossover_occurred = true
        end
    
        offspring = []
        ######## Form offspring by selecting one chromosome from each recombined or original parent
        if crossover_occurred
            chrom_from_parent1 = rand(1:2)
            chrom_from_parent2 = rand(1:2)
            push!(offspring, parent1_recombined[chrom_from_parent1])
            push!(offspring, parent2_recombined[chrom_from_parent2])
        else
            chrom_from_parent1 = rand(1:2)
            chrom_from_parent2 = rand(1:2)
            push!(offspring, parent1[chrom_from_parent1])
            push!(offspring, parent2[chrom_from_parent2])
        end

        ######## Update the chromosomes of the offspring
        updated_offspring = collect(update_chromosomes(offspring[1], offspring[2]))
        X_counts = count_TEs_per_chromosome(updated_offspring)

        ######## Accumulate generation-level counts
        generation_X_count += sum(X_counts)

        ######## Add individual counts to the DataFrame
        push!(individual_chromosome_counts, (i, X_counts...))

        push!(new_population, updated_offspring)
    end

    ######## Update the global population for the next generation
    population = new_population

    ######## Track the count of the first chromosome of the first individual
    first_individual_first_chromosome_count = count(c -> c == 'X', population[1][1])

    ######## Compute moments from the population
    # Collect total X count per individual (sum of both chromosomes)
    total_X_counts = [sum(count_TEs_per_chromosome(population[i])) for i in 1:population_size]


    ################################ Compute central of tendencies from the simulation ####################################
    # Mean:
    # mean_X_count = generation_X_count / population_size
    mean_X_count = mean(total_X_counts)

    # Variance:
    # variance_X =  mean((total_X_counts .- mean_X_count).^2)
    variance_X = max(0.0, mean((total_X_counts .- mean_X_count).^2))
    # variance_X = var(total_X_counts, corrected=false)

    # Third central moment (θ)
    mean_centered_cubed = mean((total_X_counts .- mean_X_count).^3)

    # Fourth central moment:
    mean_centered_fourth = mean((total_X_counts .- mean_X_count).^4)

    # skewness = third central moment / variance^3
    skewness_X = mean_centered_cubed / (variance_X^(3/2) + eps())  # eps() to avoid div by zero

    # Skewness Julia package
    # skewness_confirm = skewness(total_X_counts)

    # Kurtosis (κ) = fourth central moment / variance^2
    # kurtosis_X = mean_centered_fourth / (variance_X^2 + eps())  # eps() to avoid div by zero
    excess_kurtosis_X = (mean_centered_fourth / (variance_X^2 + eps())) - 3

    ## Simulated p
    P_sim = mean_X_count / variance_X


    ############################################## Theoretical computations ###############################################
    ######## Function to calculate the derivative of ln(ω(μ(x)))
    function dln_omega_dmu(μ, s)
        # return (-2 * s * μ) / (1 - s * μ^2)  # Synergistic fitness function 
        return -2 * s * μ  # Exponential fitness function
        # return (-1.5 * s * μ^0.5) / (1 - s * μ^1.5)  # Charlesworth fitness function
        # return -s 
    end

    ########## Functions to compute the second derivatives of ln(ω(μ(x)))
    function d2ln_omega_dmu2(μ, s) 
        # return (-2 * s - (2 * s^2 * μ^2)) / ((1 - (s * μ^2))^2) # Synergistic fitness function
        return -2 * s    # Exponential fitness function
        # return -(0.75 * s * (1 + 2 * s * μ^1.5)) / (μ^0.5 * (1 - s * μ^1.5)^2) # Charlesworth 
    end
    
    # Define β₁ and β₂ for μ_d
    β₁_app = dln_omega_dmu(μ_app, s)
    β₂_app = d2ln_omega_dmu2(μ_app, s)

    β₁_nb = dln_omega_dmu(μ_nb, s)
    β₂_nb = d2ln_omega_dmu2(μ_nb, s)

    # Compute p for Approximate dynamics
    p_app = μ_app / σ²_app

    # Compute ρ_nb and α_nb needed for Negative Binomial dynamics
    p_nb = μ_nb / σ²_nb
    ρ_nb = (2 - p_nb) / p_nb 
    α_nb = (6 * (1 - p_nb) + p_nb^2) / (p_nb^2)

    # Skewness (γ_nb) from equation (28): γ_nb = ρ/σ_nb
    γ_nb = ρ_nb / sqrt(σ²_nb)

    # Excess kurtosis (κ_nb - 3) from equation (29): Ex_κ_nb = α/σ²_nb
    Ex_κ_nb = α_nb / σ²_nb


    ######### Approximate μ dynamics:  dμ/dt = (u - v)μ_d + σ²_d * β₁_d
    dμ_dt_app = (u - v) * μ_app + σ²_app * β₁_app

    ######### Approximate σ² dynamics:
    # dσ²/dt = (2u + 0.5)μ_d + (u - v - 0.5)σ²_d + 0.5 * (β₁_d + β₂_d * μ_d) * σ²_d
    dσ²_dt_app = (2u + 0.5) * μ_app + (u - v - 0.5) * σ²_app + 0.5 * (β₁_app + β₂_app * μ_app) * σ²_app

    ######## Negative Binomial dynamics 
    # dμ_nb/dt = (u - v)μ_nb + σ²_nb β₁(μ_nb) + (1/2)ρ σ²_nb β₂(μ_nb)
    dμ_dt_nb = (u - v) * μ_nb + σ²_nb * β₁_nb + (1/2) * ρ_nb * σ²_nb * β₂_nb

    # dσ²_nb/dt = (2u + 1/2)μ_nb + (u - v - 1/2)σ²_nb + (1/2)(β₁(μ_nb) + β₂(μ_nb)μ_nb)σ²_nb
    #                           + (1/4)(2β₁(μ_nb) + β₂(μ_nb))ρσ²_nb + (1/4)α β₂(μ_nb)σ²_nb
    dσ²_dt_nb = (2u + 0.5) * μ_nb + (u - v - 0.5) * σ²_nb + 
                (1/2) * (β₁_nb + β₂_nb * μ_nb) * σ²_nb +
                (1/4) * (2 * β₁_nb + β₂_nb) * ρ_nb * σ²_nb +
                (1/4) * α_nb * β₂_nb * σ²_nb

    ####### Classical Mean (Ch. & Ch.) dynamics
    ## dμ/dt  = μ * (u - v) + μ * β₁ 
    dμ_dt_c = μ_c * (u - v) + μ_c * dln_omega_dmu(μ_c, s)

    # Update μ_app and σ²_app
    μ_app += dμ_dt_app * dt
    σ²_app += dσ²_dt_app * dt

    # Update μ_nb and σ²_nb
    μ_nb += dμ_dt_nb * dt
    σ²_nb += dσ²_dt_nb * dt

    # Update Ch. & Ch. mean
    μ_c += dμ_dt_c * dt


    
    

    ############################################# Print and Data storage ##############################################
    println("Generation $gen: First Individual First Chromosome X Count = $first_individual_first_chromosome_count, 
             Mean X Count = $mean_X_count, Approx mean = $μ_app, NB mean = $μ_nb, 
             Sim Variance = $variance_X, Approx variance = $σ²_app, NB variance = $σ²_nb, Charles Mean = $μ_c, 
             Sim P = $P_sim, p_app = $p_app, p_nb = $p_nb,
             Sim Skewness = $skewness_X, γ_nb = $γ_nb, 
             Sim Third Moment = $mean_centered_cubed, 
             Sim Excess kurtosis = $excess_kurtosis_X, Ex_κ_nb = $Ex_κ_nb, ")
            

    ######## Create a DataFrame to store generation-level counts
    generation_data = DataFrame(generation = [gen], 
                                First_Individual_First_Chromosome_X_Count = [first_individual_first_chromosome_count],
                                Mean_X_Count = [mean_X_count],
                                Approx_Mean = [μ_app], 
                                NB_Mean = [μ_nb], 
                                Sim_Variance = [variance_X], 
                                Approx_Variance = [σ²_app],
                                NB_Variance = [σ²_nb], 
                                Charles_Mean = [μ_c],
                                Sim_P = [P_sim], 
                                Approx_P = [p_app], 
                                NB_P = [p_nb], 
                                Sim_Skewness = [skewness_X], 
                                NB_Skewness = [γ_nb], 
                                Sim_Third_Moment = [mean_centered_cubed],
                                Sim_Excess_kurtosis = [excess_kurtosis_X], 
                                NB_Excess_Kurtosis = [Ex_κ_nb])


    ######## Write generation-level counts to a CSV file for TEs dynamics
    csv_file_path = "path/s:$(s)_N:$(population_size)_r:$(r)_beta:$(u)_delta:$(v)_exp___.csv"
    # csv_file_path = "path/moran_s:$(s)_N:$(population_size)_r:$(r)_beta:$(u)_delta:$(v)_synergistic.csv"

    if gen == 1  ## Create the file and add the header for the first generation
        CSV.write(csv_file_path, generation_data)
    else  ## Append subsequent generations without headers
        CSV.write(csv_file_path, generation_data, append=true)
    end
end


####### Time the simulation
@time begin
    ####### Display simulation parameters #######
    println("Simulation parameters:")
    println("  Chromosome length (L): $L")
    println("  Initial TEs per chromosome: $num_TEs")
    println("  Population size (N): $population_size")
    println("  Number of generations: $generations")
    println("  Transposition rate (u): $u")
    println("  Excision rate (v): $v")
    println("  Selection coefficient (s): $s")
    println("  Recombination rate (r): $r")
    println("  Fitness function: exp(-s * x_total^2)  # exponential fitness")

    println()           
    println("--------------------------------") 
    println()       

    ## Create initial population and run the simulation
    population = create_population(population_size)
    for gen in 1:generations
        generate_next_generation(gen)
    end
end

































