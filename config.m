function configurator()

clear all
    name = 21:40;
%     replicates = 20;

    for file = 1:length(name)
    %% Plankton Parameters
        %Params    
        params = struct; %organise inputs and parameters in a structure array
        params.K = 1;
        
        params.capacity = 3	; %number of species                
        params.num_demes= 9; %number of species
        params.num_resources = 3;%number of resources (determined by the number of supply vectors)
        params.total_supply = 1; %sum of all supply vectors
        params.num_planktons = 1; %number of  plankton
        params.pop_cutoff = 0.00001; %cutoff value for natural seletion
        params.shuffle_amount = 0.3; %percentage of shuffled population, needed by shuffle funtion
        params.size = [3 3]; %rows and columns of demes matrix
        params.locality = 2; %exponent of distance for shuffle transition probabilities

        %Kinetic 
        params.v = ones(params.capacity);  %value for conversion of product into biomass = v (now set to 1 for everything but can change it)
        params.mu = 0.01; %degradation rate of resources
        params.d = 0.1; %death rate of cells


        %supply matrix (initialised with spatial heterogeneity) one for every deme
        params.supply = zeros(1,params.num_resources, params.num_demes);
%         for deme = 1 : params.num_demes()
%             params.supply(:,:,deme) = rand(1, params.num_resources);
%             params.supply(:,:,deme) = params.supply(:,:,deme)/(sum(params.supply(:,:,deme))/params.total_supply);
%         end
        
        params.supply(1,:,1) = params.total_supply * [0.8 0.1 0.1];
        params.supply(1,:,2) = params.total_supply * [0.7 0.2 0.1];
        params.supply(1,:,3) = params.total_supply * [0.45 0.45 0.1];
        params.supply(1,:,4) = params.total_supply * [0.4 0.2 0.4];
        params.supply(1,:,5) = params.total_supply * [0.3 0.4 0.3];
        params.supply(1,:,6) = params.total_supply * [0.15 0.7 0.15];
        params.supply(1,:,7) = params.total_supply * [0.1 0.1 0.8];
        params.supply(1,:,8) = params.total_supply * [0.1 0.2 0.7];
        params.supply(1,:,9) = params.total_supply * [0.1 0.45 0.45];


        params.homo_supply = params.total_supply * [1/3 1/3 1/3];
       
        
        
        
        
        %Initial concentration of resources (is the same as supply)
        params.initial_concentrations = params.supply;
        
        
%         params.initial_concentrations = zeros(1, params.num_resources, params.num_demes);
%         for deme = 1:params.num_demes
%             params.initial_concentrations(:,:,deme) = params.supply(:,:,deme);
%         end



        %Total Budget
        params.total_budget = [];   %fixed enzyme budget (E)
        for s = 1: params.capacity*params.num_demes
            r = normrnd(0,0.01); %randomly distributed variable (mean of zero ad standard deviation as in posfai et al, 2017, supplementary information)
            if r >= 0   
                params.total_budget(s) = 1 + r;
            elseif r<0
                params.total_budget(s) = 1 - r;
            end
        end
        
        params.initial_populations = 1; %initial conditions for cell population to solve ODE

        %Traits
        params.strategies = zeros(1, params.num_resources, params.capacity*params.num_demes);
        for s = 1:params.capacity*params.num_demes
            params.strategies(:,:,s) = rand(1, params.num_resources);
            params.strategies(:,:,s) = params.strategies(:,:,s)/(sum(params.strategies(:,:,s))/params.total_budget(s));
        end
        
        %Mutation Rate 
        params.mutation_rate = [];
        for s = params.capacity*params.num_demes
            LOWER_BOUND_MUTATION_RATE = 0; % TO BE REFINED
            UPPER_BOUND_MUTATION_RATE = 0.1; % TO BE REFINED
            r = (UPPER_BOUND_MUTATION_RATE-LOWER_BOUND_MUTATION_RATE).*rand + LOWER_BOUND_MUTATION_RATE;
            params.mutation_rate(s) = r;
        end
            

        %% Simulation Options

        %simulation options (kinetics)
        params.v = 1;
        params.mu = 0.001;
        params.d = 0.1;
        params.time_step = 1;

        %Simulation Parameters
        params.t = 2000; %steady state value
        params.T = 200; %number of iterations (big simulation time)


         %% SAVE CONFIG FILES

        file_name = strcat('config', num2str(name(file)), '.mat')


        cd /end/home/student2/Desktop/as6616/matlab/plankton-simulator/inputs %change directory to input directory

        save(file_name,'-struct','params'); 

        cd /end/home/student2/Desktop/as6616/matlab/plankton-simulator

    end

    
end
    
    
    