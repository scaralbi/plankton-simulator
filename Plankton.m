                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
       
classdef Plankton
%     Class defining the global pool 
    
    properties
        demes = [];
        cutoff = [];
        locality = [];
    end
    
    
    methods
        %Plankton Constructor
        function obj = Plankton(deme_list, size, params_cutoff, params_locality)
            obj.demes = reshape(deme_list, size); %make plankton matrix with specified size from vector of demes
            obj.cutoff = params_cutoff;
            obj.locality = params_locality;
        end
        
  
        %Plankton Updater
        function obj = Plankton_update(obj)
            obj.demes = obj.demes;
        end        
            
        %Compute number of demes
        function n = num_demes(obj)
            n = numel(obj.demes);
        end
        
        
        %num_species = compute number of species in the plankton
        function N = num_species(obj)
            dims = obj.dim_demes();
            rows = dims(1);
            cols = dims(2);
            species_list = [];
                for row = 1 :rows
                    for col = 1 : cols
                        deme_compos = obj.demes(row, col).composition;
                        for i = 1:length(deme_compos)
                            species_list(deme_compos(i).name) = 1;
                        end
                    end
                end
                N = sum(species_list);
        end
        
        
        %dim_demes function = Return dimensions of Demes matrix
        function n = dim_demes(obj)
            n = size(obj.demes);
        end
        
        
        %update global pool
        function global_pool = global_update(obj)
            dims = obj.dim_demes();
            rows = dims(1);
            cols = dims(2);
            N = obj.num_species();
            global_pool = zeros(1, N);

            for row = 1 :rows
                for col = 1 : cols
                    length_pool = length(obj.demes(row,col).local_pool);    
                    for s = 1:N
                        for i = 1:length_pool
                            if obj.demes(row,col).composition(i).name == s
                                global_pool(1, s) = global_pool(s) + obj.demes(row,col).local_pool(i);
                            end
                        end
                    end
                end
            end
        end
        
        %update local pool of populations
        function local_pool = local_update(obj)
        dims = obj.dim_demes();
        rows = dims(1);
        cols = dims(2);
        N = obj.num_species();
        local_pool = zeros(1, N, rows*cols);
        deme = 0;
            for row = 1 : rows
                for col = 1 : cols
                        deme = deme + 1;
                        for s = 1:N
                            for i = 1:length(obj.demes(row,col).local_pool)
                                if obj.demes(row,col).composition(i).name == s
                                    local_pool(1,s,deme) = local_pool(1,s,deme) + obj.demes(row,col).local_pool(i);
                                end
                                
                            end
                        end 
                end
            end
        
        end   
        
        
        %update local pool of concentrations
        function concentration_pool = concentration_update(obj)
        dims = obj.dim_demes();
        rows = dims(1);
        cols = dims(2);
        p = obj.demes(1,1).num_resources();
        concentration_pool = zeros(1, p, rows*cols);
        deme = 0;
            for row = 1 : rows
                for col = 1 : cols
                        deme = deme + 1;
                        for resource = 1:p
                                
                            concentration_pool(1,resource,deme) = obj.demes(row,col).concentrations(resource);
                                
                            
                        end 
                end
            end
        
        end   
        
        
            
        %simulate function = Simulate demes
        function obj = simulate(obj, options)
            dims = obj.dim_demes();
            rows = dims(1);
            cols = dims(2);
            for i = 1:rows
                for j = 1:cols
                    obj.demes(i,j) = obj.demes(i,j).simulate(options);
                end  
            end
        end
        
        %natural_selection = remove species below cutoff value (unfit strategies)
        function obj = natural_selection(obj)
            dims = obj.dim_demes();
            rows = dims(1);
            cols = dims(2);
            for row = 1:rows
                for col = 1:cols
                    for species = 1: obj.demes(row,col).num_species()
                        if  obj.demes(row,col).local_pool(species) < obj.cutoff
                            obj.demes(row,col).local_pool(species) = 0; %set population values below cutoff to 0
                            obj.demes(row,col).local_traits(species, :) = 0; %as above but for traits
                            obj.demes(row,col).composition(species).population = 0;
                        end
                    end
                end
            end 
        end
        
    %local shuffle 
    function [obj, shuffle_log] = local_shuffle(obj,amount)
        dims = obj.dim_demes(); %dimension of plankton
        rows = dims(1);
        cols = dims(2);
        demes = obj.num_demes()
        N = obj.num_species(); %number of species, needs to be generaised
        invaders_name = zeros(1,N,demes); %names of shuffled species (for indexing)
        invaders_pool = zeros(1,N,demes); %population values of shuffled species
        
        %%Invader Selector
        deme = 0; %init deme count
        for row= 1:rows 
            for col= 1:cols
                weights = zeros(1, length(obj.demes(row,col).local_pool)); %init weights for each species
                deme = deme + 1; %count which deme you are working on
                pool_sum = sum(obj.demes(row,col).local_pool); %sum of population values to normalise (so that sum of each row adds up to 1)
                for i = 1:length(obj.demes(row,col).local_pool)
                    weights(i) = obj.demes(row,col).local_pool(i) / pool_sum; %calculate probabilities (normalised populations)
                end
                r = randsample(1:length(obj.demes(row,col).local_pool), 1, true, weights); %samples randomply which species to be an invader weighted by population values
                invaders_name(1, obj.demes(row,col).composition(r).name, deme) = r; %return the index of the invader
                invaders_pool(1, obj.demes(row,col).composition(r).name, deme) = obj.demes(row,col).local_pool(r); %add invader population value into invader pool (useful to take it back)
            end
        end
        
        %%Deme Selector
        %linear to subsrcipt indexes
        IND = 1 : obj.num_demes();
        dims = obj.dim_demes();
        [col_index,row_index] = ind2sub(dims,IND);
        %row_index = [1 1 1 2 2 2 3 3 3]; %row indexes of each deme in the plankton (needs to be generalised)         
        %col_index = [1 2 3 1 2 3 1 2 3]; %column indexes
        Tprob = zeros(rows*cols, rows*cols); %init transition probabilities matrix
        Distance = zeros(rows*cols, rows*cols); %init Distance matrix

        %Calculate Euclidean Distance between Demes
        for pre_deme = 1 : rows*cols 
            x = row_index(pre_deme); %compute x coordinate of donor deme
            y = col_index(pre_deme); % compute y coordinate of donor deme 
            for post_deme = 1:rows*cols
                 X = row_index(post_deme); %compute x coordinate of recipient deme
                 Y = col_index(post_deme);%compute y coordinate of recipient deme   
                 Distance(pre_deme, post_deme) = sqrt((X - x)^2 + (Y- y)^2); %compute euclidean Distance (Distance = sqrt(X_diff^2 + Y_diff^2)     
            end
        end
        
         %Calculate Transition Probabilities from Euclidean Distances
         for row = 1 : rows*cols %iterate for rows
            for col = 1 : rows*cols %iterate for cols
                Tprob(row,col) = Distance(row,col)^obj.locality; %distance to the power of locality
                Tprob(row,col) = 1/Tprob(row,col); %probability of being shuffled is inversely poportional to distance 
                if Tprob(row,col) == inf %if infinity (same demes transition)
                    Tprob(row,col) = 0; %set to zero (no probs for same deme shuffles)
                end
            end
        end
        
        %normalise prob so that sum of all probs is equal to 1 
        for row = 1 : rows*cols
            sum_row = sum(Tprob(row,:));
            for col = 1 : rows*cols
                Tprob(row,col) = Tprob(row,col) / sum_row;
            end
        end
        
        %Draw Recipient Deme
        shuffle_log = zeros(demes,demes);
        deme = 0;    
        for row= 1:rows
            for col = 1:cols
                deme = deme + 1;
                donor = deme;
                rand_deme = randsample(1:rows*cols, 1, true, Tprob(deme,:)); %draw recipient deme
                recipient = rand_deme;
                rand_row = row_index(rand_deme);%recipient row
                rand_col = col_index(rand_deme);%recipient col
        %Update Recipient Properties            
            for i = 1 : N
                if invaders_name(1,i,deme) > 0
                   obj.demes(rand_row,rand_col).local_traits = [obj.demes(rand_row,rand_col).local_traits; obj.demes(row, col).local_traits(invaders_name(1,i,deme),:)]; %append row of traits for randomly drawn species
                   obj.demes(rand_row,rand_col).composition = [obj.demes(rand_row,rand_col).composition obj.demes(row, col).composition(invaders_name(1,i,deme))]; 
                   obj.demes(rand_row,rand_col).local_pool = [obj.demes(rand_row,rand_col).local_pool amount*obj.demes(row,col).local_pool(invaders_name(1,i,deme))]; %append population  for randomly drawn species
        
           %Update Donor Properties            
                   %substract invader population from original deme (conservation of mass)
                   obj.demes(row,col).local_pool(invaders_name(1,i,deme)) = obj.demes(row,col).local_pool(invaders_name(1,i,deme)) - amount*obj.demes(row,col).local_pool(invaders_name(1,i,deme)); 
                
                   shuffle_log(donor, recipient) = i;
                end
            end

            end
      
        end
  
    end
        
    %Well Mixer Function = make nutrient concentration homogenous in all demes
    function obj = mix(obj, homo_supply)
        dims = obj.dim_demes();
        rows = dims(1);
        cols = dims(2);
        for i = 1:rows
            for j = 1:cols
                obj.demes(i,j).fluxes = homo_supply;
            end  
        end
    end
        
    
 %UnMixer Function = make nutrient concentration heterogeneous in all demes
    function obj = unmix(obj, supply)
        dims = obj.dim_demes();
        rows = dims(1);
        cols = dims(2);
        deme = 0;
        for i = 1:rows
            for j = 1:cols
                deme = deme + 1
                obj.demes(i,j).fluxes = supply(1,:,deme);
            end  
        end
    end
    
    
     % EVOLVE
    % = Evolve demes
        function obj = evolve(obj)
     
            % 1. Unite all species into global pool
            species_list = obj.composition();

            % 2. Evolve all species
            mutatedSpecies_list = obj.composition();
            for i = 1:length(species_list)
                [old_species, new_species] = species_list(i).evolve(length(species_list)+i-1);
                species_list(i) = old_species;
                mutatedSpecies_list(i) = new_species;
            end

            % 3. Reassign old and new species to demes
            dims = obj.dim_demes();
            rows = dims(1);
            cols = dims(2);
            for i = 1:rows
                for j = 1:cols
                    deme_comp = obj.demes(i,j).composition;

                    for s = randi(length(deme_comp))
                        species_name = deme_comp(s).name;

                        contribution = 0;
                        if species_list(species_name).population + mutatedSpecies_list(species_name).population > 0
                            contribution = deme_comp(s).population / (species_list(species_name).population + mutatedSpecies_list(species_name).population);
                        end

                        old_species_left = species_list(species_name).population * contribution;
                        new_species_add = mutatedSpecies_list(species_name).population * contribution;
                        obj.demes(i,j).local_pool(species_name) = old_species_left;
                        obj.demes(i,j).local_pool = [obj.demes(i,j).local_pool new_species_add];
                        obj.demes(i,j).composition = [obj.demes(i,j).composition species_update(mutatedSpecies_list(species_name), new_species_add)];
                        obj.demes(i,j).local_traits = [obj.demes(i,j).local_traits; mutatedSpecies_list(species_name).traits];

                    end
                end  
            end
        end 
       
    end % methods   
end %Class