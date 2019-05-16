function return_vals = dynamics(t, current_vals, params)
    %Initialise lists for data  
    current_population =((current_vals(1 : params.num_species)));
    current_concentration = (current_vals(params.num_species+1 : length(current_vals)));

    %Initialise and fill Uptake Rate Function
    uprate = zeros(1, params.num_resources);
    for resource = 1 : params.num_resources
        concentration = current_concentration(resource);
        uprate(resource) = concentration / (params.K + concentration); % define uptake rates (Monod function) for each metabolite
    end

    %Initialise and fill Growth Rate Function
    growth_rate = zeros(params.num_species, params.num_resources);
    for species = 1 : params.num_species %define growth rate function for each cell type 
        for resource = 1 : params.num_resources %on each metabolite
%             growth_rate(species, resource) = params.v(species) * params.strategies(species, resource) * uprate(resource);
              growth_rate(species, resource) = params.v * params.strategies(species, resource) * uprate(resource);

        end
    end


    % Population differential

    dndt = zeros(1, params.num_species); %initialise list for dndt

    for species = 1 : params.num_species
        growth_rate_sum = 0;
        for resource = 1 : params.num_resources
            growth_rate_sum = growth_rate_sum + growth_rate(species, resource); %sum all growth rate contributions for each cell on each resource
        end
        population = current_population(species);
        deaths = params.d * population;
        dndt(species) = growth_rate_sum * population - deaths; %growth dynamics, as in (Ponsfai et al., 2018)
    end


    % Concentration differential
    dcdt = zeros(1, params.num_resources); 

    for resource = 1 : params.num_resources
        total_consumption = 0;
        for species = 1 : params.num_species
            consumption = current_population(species) * params.strategies(species, resource); %compute amount of resource consumed by each species
            total_consumption = total_consumption + consumption;
        end

        concentration = current_concentration(resource);
        degradation = params.mu * concentration;
        dcdt(resource) = params.supply(resource) - total_consumption * uprate(resource) - degradation; 
    end

    return_vals = transpose([dndt,  dcdt]);

end