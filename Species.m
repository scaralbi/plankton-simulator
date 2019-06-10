classdef Species
    % Class defining a single deme
    
    properties
        name = [];
        traits = [];
        population = [];
        mutation_rate = [];
    end
    
    methods
        % Constructor
        function obj = Species(input_name, input_traits, input_pop, mutation_rate)
            obj.name = input_name;
            obj.traits = input_traits;
            obj.population = input_pop;
            obj.mutation_rate = mutation_rate;
        end
            
        %Updater
        function obj = species_update(obj, global_pool)
            obj.name = obj.name;
            obj.traits = obj.traits;
            obj.mutation_rate = obj.mutation_rate;
            obj.population = global_pool;%global pool has to be indexed
        end
        
            
        
         %Evolver = randomly change traits
        function new_traits = evolve_traits(obj)
            
            old_traits = obj.traits;
            new_traits = [];
            for trait = 1:length(old_traits)
                LOWER_BOUND_EVOLUTION_RATE = -0.2; % TO BE REFINED
                UPPER_BOUND_EVOLUTION_RATE = 0.2; % TO BE REFINED
                update = (UPPER_BOUND_EVOLUTION_RATE-LOWER_BOUND_EVOLUTION_RATE).*rand + LOWER_BOUND_EVOLUTION_RATE;
                new_traits(trait) = max(min(old_traits(trait) + update, 1), 0);
            end
        end
        
        
        %implement speciation
        function [species_out1, species_out2] = evolve(species_in,N)
            % Define population split
            population1 = species_in.population * (1-species_in.mutation_rate);
            population2 = species_in.population * species_in.mutation_rate;
            % Generate mutated traits
            new_traits = evolve_traits(species_in);
            new_budget = sum(new_traits);
            
            % Generate new species
            species_out1 = species_update(species_in, population1);
            species_out2 = Species(N+1, new_traits, population2, species_in.mutation_rate);
        end
       
    end % methods
end %Class


