classdef Deme
    % Class defining a single deme
    
    properties
        concentrations = [];
        fluxes = [];
        composition = [];
        local_pool = [];
        local_traits = [];
        capacity = [];
    end
    
    methods
        % Constructor
        function obj = Deme(input_conc, input_fluxes, input_species, input_capacity)
            obj.concentrations = input_conc;
            obj.fluxes = input_fluxes;
            obj.composition = input_species;
            obj.capacity = input_capacity;
            
            
            for i = 1: length(input_species)
                obj.local_pool(i) = input_species(i).population;
            end
            
            for i = 1: length(input_species)
                obj.local_traits(i,:) = input_species(i).traits;
            end    
        end
        
        %Copty onstructor with different supply
        function obj = Deme_copy(obj, new_supply)
            obj.fluxes = new_supply;
        end
            
        %num_species = compute number of species in a deme
        function n = num_species(obj)
            n = length(obj.composition);
        end
        
        %num_resources = compute number of resources within a deme
        function n = num_resources(obj)
            n = length(obj.concentrations);
        end
        
        
           
        %natural_selection = remove species below cutoff value (unfit strategies)
        function [obj, survivors] = natural_selection(obj, cutoff)
            
            for species = 1: num_species(obj)
                if obj.local_pool(species) < cutoff
                    obj.local_pool(species) = 0; %set population values below cutoff to 0
                    obj.local_traits(species, :) = 0; %as above but for traits
                    obj.composition(species).population = 0;
                end
            end
            
            survivors = [];
            
            for i= 1 : obj.num_species()
                if obj.local_pool(i) > 0
                survivors(i) = obj.composition(i).name;
                end
            end
            survivors = char(survivors);
        end
      
        function obj = simulate(obj, options)
            %% ODE Params
            ode_params = struct;
     
            ode_params.strategies = obj.local_traits;
            ode_params.v = options.v; 
            ode_params.K = options.K;  
            ode_params.mu = options.mu; 
            ode_params.d = options.d; 
            ode_params.supply = obj.fluxes;
            ode_params.num_species = obj.num_species();
            ode_params.num_resources = obj.num_resources();
            sim_time = options.t;
       
            %% ODE SOLVER
            simulation_time = 1 : options.time_step : sim_time ;
           
            initial_state = [obj.local_pool, obj.concentrations]; %array of initial conditions
            [ t, Y ] = ode15s(@(t, y)dynamics(t, y, ode_params), simulation_time, initial_state); %ODE solver
            
            %% Extract output data
            populations_log = Y(:,1:obj.num_species()); %retrieve n(t) from Y
            concentrations_log = Y(:, (obj.num_species() + 1) : (obj.num_species()+obj.num_resources()) ); %retrieve c(t) from Y
            
            %% Update deme
            
            obj.local_pool = populations_log(end, : );
            
            obj.concentrations = concentrations_log(end, :);
            
        end % simulate
         
    end % methods
end %Class


