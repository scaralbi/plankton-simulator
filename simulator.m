function plankton_sim()

 clear all
    
 %% Read Config Files
 
    cd /end/home/student2/Desktop/as6616/matlab/plankton-simulator/inputs %change diretory to config_files
  
    listing = dir('*config*.mat') %find config files in directory
    
    num_files = size(listing, 1) %count number of files (just number of rows)
    
    myVars = {'supply', 'homo_supply', 'num_resources', 'size', 'shuffle_amount', 'num_planktons','num_demes', 'total_supply', 'pop_cutoff', 'capacity', 'time_step', 'total_budget', 'initial_populations', 'locality', 'initial_concentrations', 't', 'T', 'v', 'K', 'mu', 'd', 'strategies'};
    t = datetime('now','Format','yy.MM.dd''_''HH.mm');

    S = char(t);
    folder_name = ['sym_',S]

    %iterate for every plankton
    for file = 1 :num_files
        
        cd /end/home/student2/Desktop/as6616/matlab/plankton-simulator/inputs %change diretory to config_files

        file_name = strcat('config', num2str(file), '.mat')
        params = load(file_name, myVars{:}); %import parameter structures
   

     %Inputs 
        cd /end/home/student2/Desktop/as6616/matlab/plankton-simulator %change diretory to config_files

        %Init species
        species = [];
        for i = 1:params.capacity*params.num_demes
            species = [species Species(i, params.strategies(:,:,i), params.initial_populations, params.total_budget)];
        end

        %Init demes
        demes = [];
        species_list = [];
        for i = 1:params.num_demes
            base_offset = params.capacity*(i-1);
            for s = 1 :params.capacity
               species_list = [species_list species(base_offset +s )];
            end
            demes = [demes Deme(params.initial_concentrations(:,:,i), params.supply(:,:,i), species_list, params.capacity)];
        end

        %Init Plankton
        plankton = Plankton(demes, params.size, params.pop_cutoff, params.locality);

        %Init Local and Global Pool - N = total number of species
        N = params.capacity*params.num_demes;
        p = params.num_resources;
        global_pool = zeros(params.T, N);
        temp_global_pool = zeros(params.T, N);
        global_pool(1,:) = params.initial_populations*ones(1,N);
        deme_pool = zeros(params.T, N, params.num_demes);
        resource_pool = zeros(params.T, p, params.num_demes);

        shuffle_history = zeros(params.num_demes,params.num_demes);

        %init initial populations in each deme
        
        for deme = 1:params.num_demes
            deme_pool(1, deme*params.capacity-params.capacity+1:deme*params.capacity, deme) = params.initial_populations;
        end        
        
        
        % Waitbar for execution time estimate
%         w = waitbar(0,'Wait...Plankton is growing');

        %% Iterate

        for iteration = 2:params.T
            
            

        
            plankton.demes(1,1).fluxes
            plankton.demes(1,1).concentrations
            
            %simulate plankton   
            plankton = simulate(plankton, params);

      

            %update global pool
            global_pool(iteration, :) = global_update(plankton);

            %update local pools for each deme         
            deme_pool(iteration, :, :) = local_update(plankton);
            
            %update resource pools for each deme         
            resource_pool(iteration, :, :) = concentration_update(plankton);

            %update species population 
            for s = 1:N
                species(s) = species_update(species(s), global_pool(iteration, s));
            end

            %set to 0 traits, localpool and composition of species below cutoff
            plankton = natural_selection(plankton);


            %update temp global pool
            temp_global_pool(iteration, :) = global_update(plankton);
            
            %update species population
            for s = 1:N
                species(s) = species_update(species(s) ,temp_global_pool(iteration, s));
            end

            %Update Plankton Plankton        
            plankton = Plankton_update(plankton);

            %shuffle plankton
            [plankton, shuffle_log] = local_shuffle(plankton, params.shuffle_amount);
            
            shuffle_history = cat(3,shuffle_history,shuffle_log);
            
            if 500 <= iteration && iteration <= 1000 
                plankton = mix(plankton, params.homo_supply);
                
            elseif iteration > 1000
                plankton = unmix(plankton, params.supply);
                
            end
            
      

            %show waiting bar
%             waitbar(iteration/params.T)
        end

     %% Outputs 
 shuffle_history
       
        
%         close(w)
        
        %% DATA
        
    %Compute surviving species
    surviving = 0;
    for s = 1:N
        if global_pool(end,s) > 3
        surviving = surviving + 1;
        end
    end
    
    %Compute the fittest species
    [M, I] = max(global_pool(end,:));
    best = I;
    best_trait = params.strategies(1,:,best);
 
    %Compute Shannon's Index
    H = zeros(1,N);
    Ntot = 0;
    for s = 1:N
        Ntot = Ntot + global_pool(params.T,s);
    end
    
    for s= 1:N 
       p = global_pool(params.T,s)/Ntot;
       H(1,s) = - p * log(p);
    end
   
    H(isnan(H))=0;
    
    Htot = sum(H(1,:))   ;
    
    %Output data
    Data = {'Demes', 'n', 'N', 'l', 'p', 'E', 'S', 'cutoff', 'amount', 'Surviving', 'Best Traits', 'H', 'shuffle_history', 'strategies', 'global pool', 'local pool';
            params.size, params.capacity, N, params.locality, params.num_resources, params.total_budget, params.total_supply, params.pop_cutoff, params.shuffle_amount, surviving, best_trait, Htot, shuffle_history, params.strategies, global_pool, deme_pool};    
    
        %% FIGURES
             
%    Global pool plot
    f1 = figure;
%     set(f1, 'Visible', 'off');
    h1 = plot(0:params.T-1, global_pool(:,:), 'LineWidth', 2); hold on
    set(h1, {'color'}, num2cell(hsv(N),2));
    xlabel('Time (iterations)');
    ylabel('Cell Abundance');
    ax = gca;
    ax.XLim = [0 params.T-1];
    ylim = max(global_pool(:));
    ax.YLim = [0 ylim];
    ax.FontSize = 14;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridLineStyle = ':';
    L = 1:N;
    legendCell = cellstr(num2str(L', '%-d'));
    hleg = legend(legendCell);
    htitle = get(hleg,'Title');
    hleg.Location = 'bestoutside';
    set(htitle,'String','Species');
    box off
    title('Simulated Plankton Dynamics'); hold off
    
    
%      plot demes dynamics
    f2 = figure('units','normalized','outerposition',[0 0 1 1]);
%     set(f2, 'Visible', 'off');
    for deme = 1:params.num_demes; hold on
        subplot(3, 3, deme);
        h2 = plot(0:params.T-1, deme_pool(:,:,deme), 'LineWidth', 2);
        set(h2, {'color'}, num2cell(hsv(N),2));
        xlabel('Time (iterations)');
        ylabel('Cell Abundance');
        ax = gca;
        ax.XLim = [0 params.T-1];
        ylim = max(deme_pool(:));
        ax.YLim = [0 ylim];
        ax.FontSize = 14;
        ax.YTick = 0 : 5: ylim;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridLineStyle = ':';
        box off
        title(['Deme' num2str(deme)]); 
    end
    hold off
    
    
    %      plot demes dynamics
    f3 = figure('units','normalized','outerposition',[0 0 1 1]);
%     set(f3, 'Visible', 'off');
    for deme = 1:params.num_demes; hold on
        subplot(3, 3, deme);
        h3 = plot(0:params.T-1, resource_pool(:,:,deme), 'LineWidth', 2);
        xlabel('Time (iterations)');
        ylabel('Nutrient Conentration');
        legend('i = 1', 'i = 2', 'i = 3')
        ax = gca;
        ax.XLim = [0 params.T-1];
        ylim = max(resource_pool(:));
        ax.YLim = [0 ylim];
        ax.FontSize = 14;
        ax.YTick = 0 : 5: ylim;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridLineStyle = ':';
        box off
        title(['Deme' num2str(deme)]); 
    end
    hold off
%     
% %     plot simplex
%     f4 = figure;
%     set(f4, 'Visible', 'off');
%     deme_color = summer(params.num_demes); %demes colromap
%     species_color = hsv(N); %species colormap
%     for deme= 1:params.num_demes
%         hAxis = subplot(3,3,deme);
%         su = scatter3(params.supply(1,1,deme)/params.total_supply, params.supply(1,2,deme)/params.total_supply, params.supply(1,3,deme)/params.total_supply, 100, 'd', 'filled', 'k'); hold on
%         f = fill3([1 0 0],[0 1 0],[0 0 1],deme_color(deme,:)); hold on 
%         for row = 1:params.size(1)
%             for col = 1:params.size(2)
%                 for i = 1 : N
%                     if deme_pool(end,i,deme) > params.pop_cutoff
% %                     h = scatter3(params.strategies(1,1,i),params.strategies(1,2,i), params.strategies(1,3,i), 70, 'o', 'filled'); hold on
%                     h = scatter3(species(i).traits(1,1),species(i).traits(1,2), species(i).traits(1,3), 70, 'o', 'filled'); hold on
%                     h.MarkerFaceColor = species_color(i,:);
%                     h.MarkerEdgeColor = species_color(i,:);                    
%                     end
%                 end
%             end
%         end
%         f.EdgeColor = 'none';
%         f.FaceAlpha = 0.3;
%         hAxis.XAxis.FirstCrossoverValue  = hAxis.YLim(1); 
%         hAxis.XAxis.SecondCrossoverValue = hAxis.ZLim(1);
%         hAxis.YAxis.FirstCrossoverValue  = hAxis.XLim(1);
%         hAxis.YAxis.SecondCrossoverValue = hAxis.ZLim(1);
%         hAxis.ZAxis.FirstCrossoverValue  = hAxis.XLim(1); 
%         hAxis.ZAxis.SecondCrossoverValue = hAxis.YLim(1);      
%         hAxis.YLim =[0 1];
%         hAxis.ZLim =[0 1]; 
%         hAxis.XLim =[0 1]; 
%         xlabel('\alpha_1');
%         ylabel('\alpha_2');
%         zlabel('\alpha_3');
%         hAxis.XDir = 'reverse';
%         hAxis.YDir = 'reverse';
%         grid on
%         box off
%         %axis off
%     end
%     hold off

 
%     Save files to outputs directory
    
    file_name = ['config' num2str(file) '_data'];
%     global_name = strcat('global', num2str(file),'.png');
%     local_name = strcat('local', num2str(file), '.png');
%     simplex_name = strcat('simplex', num2str(file), '.png');

    cd /end/home/student2/Desktop/as6616/matlab/plankton-simulator/outputs %open output directory
    
    if file == 1 
        mkdir (folder_name)
        cd (folder_name) %save data to  subdirectory

    else
        cd (folder_name) %save data to  subdirectory
    end
    

       
    save(file_name); %make new subdirectory in output directory



    cd /end/home/student2/Desktop/as6616/matlab/plankton-simulator %return to current directory

        
    end
    

end

  