classdef Species
    % Class defining a single deme
    
    properties
        name = [];
        traits = [];
        population = [];
        budget = [];
    end
    
    methods
        % Constructor
        function obj = Species(input_name, input_traits, input_pop, input_budget)
            obj.name = input_name;
            obj.traits = input_traits;
            obj.population = input_pop;
            obj.budget = input_budget;
        end
            
        %Updater
        function obj = species_update(obj, global_pool)
            obj.name = obj.name;
            obj.traits = obj.traits;
            obj.budget = obj.budget;
            obj.population = global_pool;%global pool has to be indexed
        end
        
         %Evolver = randomly change traits
        function obj = evolve(obj)
            a = -0.1;
            b = 0.3;
                r = (b-a).*rand + a;
            obj.traits = obj.traits + r;
     
        end
       
    end % methods
end %Class


