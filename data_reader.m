
function datagather()
clear all

    cd /end/home/student2/Desktop/as6616/matlab/plankton-simulator/outputs/sym_19.05.15_12.45 %change diretory to config_files

    replicates = 20;



    listing = dir('*config*.mat'); %find config files in directory
    
    num_files = size(listing, 1) %count number of files (just number of rows)
%     myVars = {'params'};

    cases = num_files/replicates;
    H = [];
    num_resources = [];
    surviving = [];
    N = 27;
    rank_plot = zeros(N, 2, num_files);


    

    %iterate for every config file
    for file = 1:num_files
        file_name = strcat('config', num2str(file), '_data.mat')
%         data = load(file_name); %import parameter structures
        
        load(file_name, 'params', 'Data', 'global_pool', 'deme_pool', 'shuffle_history');
        H(file)  = Data{2,12};
        
        surviving(file) = 0;
        for s = 1:N
            if global_pool(end,s) > 3
                surviving(file) = surviving(file) + 1;
            end
        end
        
%         num_resources(file) = params.num_resources
        
%            %Rank Abundance Plot
%        rank_plot(:,1,file) = 1 : N;
%        tot_pop = sum(global_pool(end,:,end));
%        steady_pop = global_pool(end,:,end);
%        sort_pop = sort(steady_pop, 'descend');
%        rel_pop = sort_pop./tot_pop;
%        for i = 1: length(rel_pop)
%            rank_plot(i,2,file) = rel_pop(i);
%        end
% %        rank_plot(:,2,file) = rel_pop';   
%        Rank = 1:Data{2,3};
%        Rank = mat2cell(Rank,1,Data{2,3});
   
   
%     %compute invasion attempts for eah species
%     invasions = zeros(params.T,Data{2,3});
%     for T = 1:params.T
%         for species = 1:Data{2,3}
%             if find(shuffle_history(:,:,T)==species)
%                 [row, col] = find(shuffle_history(:,:,T)==species)
%                 invasions(T,species) = invasions(T,species)+1;
%             end
%         end
%     end
%     invasions
        
        
% % %         Global pool plot
%     f1 = figure;
% %     set(f1, 'Visible', 'off');
%     h1 = plot(0:params.T-1, global_pool(:,:), 'LineWidth', 2); hold on
%     set(h1, {'color'}, num2cell(hsv(Data{2,3}),2));
%     xlabel('Time (iterations)');
%     ylabel('Cell Abundance');
%     ax = gca;
%     ax.XLim = [0 params.T-1];
%     ylim = max(global_pool(:));
%     ax.YLim = [0 ylim];
%     ax.FontSize = 14;
%     ax.XGrid = 'on';
%     ax.YGrid = 'on';
%     ax.GridLineStyle = ':';
%     L = 1:Data{2,3};
%     legendCell = cellstr(num2str(L', '%-d'));
%     hleg = legend(legendCell);
%     htitle = get(hleg,'Title');
%     hleg.Location = 'bestoutside';
%     set(htitle,'String','Species');
%     box off
%     title('Simulated Plankton Dynamics'); hold off
%     
%     
% %      plot demes dynamics
%     f2 = figure('units','normalized','outerposition',[0 0 1 1]);
% %     set(f2, 'Visible', 'off');
%     for deme = 1:params.num_demes; hold on
%         subplot(5, 5, deme);
%         h2 = plot(0:params.T-1, deme_pool(:,:,deme), 'LineWidth', 2);
%         set(h2, {'color'}, num2cell(hsv(Data{2,3}),2));
%         xlabel('Time (iterations)');
%         ylabel('Cell Abundance');
%         ax = gca;
%         ax.XLim = [0 params.T-1];
%         ylim = max(deme_pool(:));
%         ax.YLim = [0 ylim];
%         ax.FontSize = 14;
%         ax.YTick = 0 : 5: ylim;
%         ax.XGrid = 'on';
%         ax.YGrid = 'on';
%         ax.GridLineStyle = ':';
%         box off
%         title(['Deme' num2str(deme)]); 
%     end
%     hold off
%     
%     plot simplex
    f3 = figure;
%     set(f3, 'Visible', 'off');
    iteration = 750;
    deme_color = summer(params.num_demes); %demes colromap
    species_color = hsv(Data{2,3}); %species colormap
    for deme= 1:params.num_demes
        hAxis = subplot(3,3,deme);
%         su = scatter3(params.supply(1,1,deme)/params.total_supply, params.supply(1,2,deme)/params.total_supply, params.supply(1,3,deme)/params.total_supply, 100, 'd', 'filled', 'k'); hold on
        su = scatter3(params.homo_supply(1,1)/params.total_supply, params.homo_supply(1,2)/params.total_supply, params.homo_supply(1,3)/params.total_supply, 100, 'd', 'filled', 'k'); hold on
        f = fill3([1 0 0],[0 1 0],[0 0 1],deme_color(deme,:)); hold on 
        for row = 1:params.size(1)
            for col = 1:params.size(2)
                for i = 1 : Data{2,3}
                    if deme_pool(iteration,i,deme) > 1
%                     h = scatter3(params.strategies(1,1,i),params.strategies(1,2,i), params.strategies(1,3,i), 70, 'o', 'filled'); hold on
                    h = scatter3(params.strategies(1,1,i), params.strategies(1,2,i), params.strategies(1,3,i), 70, 'o', 'filled'); hold on
                    h.MarkerFaceColor = species_color(i,:);
                    h.MarkerEdgeColor = species_color(i,:);                    
                    end
                end
            end
        end
        f.EdgeColor = 'none';
        f.FaceAlpha = 0.3;
        hAxis.XAxis.FirstCrossoverValue  = hAxis.YLim(1); 
        hAxis.XAxis.SecondCrossoverValue = hAxis.ZLim(1);
        hAxis.YAxis.FirstCrossoverValue  = hAxis.XLim(1);
        hAxis.YAxis.SecondCrossoverValue = hAxis.ZLim(1);
        hAxis.ZAxis.FirstCrossoverValue  = hAxis.XLim(1); 
        hAxis.ZAxis.SecondCrossoverValue = hAxis.YLim(1);      
        hAxis.YLim =[0 1];
        hAxis.ZLim =[0 1]; 
        hAxis.XLim =[0 1]; 
        xlabel('\alpha_1');
        ylabel('\alpha_2');
        zlabel('\alpha_3');
        hAxis.XDir = 'reverse';
        hAxis.YDir = 'reverse';
        grid on
        box off
        %axis off
    end
    hold off
       

    end
    
    cd /end/home/student2/Desktop/as6616/matlab/plankton-simulator
    
   
%     %Rank Abundance plot   
%    std_1 = [];
%    mean_1 = [];
%    
%    std_2 = [];
%    mean_2 = [];
%    
%    std_3 = [];
%    mean_3 = [];
%    
%    std_4 = [];
%    mean_4 = [];
%    
%     
%    std_5 = [];
%    mean_5 = [];
%    
%    for rank = 1 : 36
%            std_1(rank) = std(rank_plot(rank,2,1:20));
%            mean_1(rank) = mean(rank_plot(rank,2,1:20));
%            
%            std_2(rank) = std(rank_plot(rank,2,21:40));
%            mean_2(rank) = mean(rank_plot(rank,2,21:40));
%            
%            std_3(rank) = std(rank_plot(rank,2,41:60));
%            mean_3(rank) = mean(rank_plot(rank,2,41:60));
%            
%            std_4(rank) = std(rank_plot(rank,2,61:80));
%            mean_4(rank) = mean(rank_plot(rank,2,61:80));
%            
%            
%            std_5(rank) = std(rank_plot(rank,2,81:100));
%            mean_5(rank) = mean(rank_plot(rank,2,81:100));
%    end
   
%    mean_global = log(mean_global)
%    mean_local = log(mean_local)
%    std_global = log(std_global)
%    std_local = log(std_local)

%    
%    figure()
%    e1 = errorbar(mean_1,std_1, '-o','MarkerSize',7, 'MarkerEdgeColor','blue','MarkerFaceColor','blue', 'LineWidth', 1.5);hold on
%    e2 = errorbar(mean_2,std_2, '-s','MarkerSize',7, 'MarkerEdgeColor','red','MarkerFaceColor','red', 'LineWidth', 1.5);
%    e3 = errorbar(mean_3,std_3, '-d','MarkerSize',7, 'MarkerEdgeColor','g','MarkerFaceColor','g', 'LineWidth', 1.5);
%    e3.Color = 'g';
%    e4 = errorbar(mean_4,std_4, '-p','MarkerSize',7, 'MarkerEdgeColor','m','MarkerFaceColor','m', 'LineWidth', 1.5);
%    e4.Color = 'm'; 
%    e5 = errorbar(mean_5,std_5, '-h','MarkerSize',7, 'MarkerEdgeColor',[255/255,165/255,0],'MarkerFaceColor',[255/255,165/255,0], 'LineWidth', 1.5);
%    e5.Color = [255/255,165/255,0]; 
%    legend('i = 1', 'i = 2', 'i = 3', 'i = 4', 'i = 5');
%    title('Simulated Rank Abundance with Different Number of Resources')
%    xlabel('Species Rank')
%    ylabel('Relative Abundance')
%    ax1 = gca
%    ax1.YScale = 'log';
%    ax1.FontSize = 18;
%    ax1.XLim = [1 10];
%    ax1.XTick = 1:36
%    box off
%    grid on
%    hold off
% %    
%    
%  
%    for i = 21 : 40
%        semilogy(rank_plot(:,1,i), rank_plot(:,2,i));hold on
%    end
%     xlabel('Species Rank')
%     ylabel('Relative Abundance')
%     grid on
%     ax = gca
%     ax.XLim = [1 36]
%     ax.XTick = 1:36
%     box off
%     title('Rank Abundance Plot with Global Shuffle')
%     hold off
%   
    
    
%     
%     Rank_shade = zeros(27,20,2);
%     for obs = 1:27
%         Rank_shade(obs,:,1) = rank_plot(obs,2,1:20); %global shuffle
%         Rank_shade(obs,:,2) = rank_plot(obs,2,41:60) ;%local shuffle
% 
%     end
    
    
%     Rank_shade = log(Rank_shade);
%    figure()
% % Prepare data
%     y_global= Rank_shade(:,:,1)';
%     y_local= Rank_shade(:,:,2)';
% 
%     x = 1:27;
% 
%     % Make the plot
%     clf
%     s1 = shadedErrorBar(x,y_global,{@mean,@std},'lineprops', '-b', 'patchSaturation',0.2); hold on
%     s2 = shadedErrorBar(x,y_local,{@mean,@std}, 'lineprops', '-g', 'patchSaturation',0.2); 
%     set(gca, 'YScale', 'log')
%     ax = gca
%     ax.XLim = [1 27]
%     ax.YLim = [1.0000e-06 1]
%     s1.mainLine.LineWidth = 2;
%     s1.patch.FaceColor = 'b';
%     s2.mainLine.LineWidth = 2;
%     s2.patch.FaceColor = 'g';
% 
% 
%     grid on
% 
% 
%    
%    figure()
%    for i = 41 : 60
%        semilogy(rank_plot(:,1,i), rank_plot(:,2,i));hold on
%    end
%     xlabel('Species Rank')
%     ylabel('Relative Abundance')
%     grid on
%     ax = gca
%     ax.XLim = [1 27]
%     ax.XTick = 1:27
%     box off
%     title('Rank Abundance Plot with Local Shuffle')
%     hold off
%    



%     Entropy = zeros(replicates,cases);
%     for col = 1:cases
%         Entropy(:,col) = H(col*replicates - replicates +1 : col*replicates);
%     end
%     
%     Surviving = zeros(replicates,cases);
%     for col = 1:cases
%         Surviving(:,col) = surviving(col*replicates - replicates +1 : col*replicates);
%     end
%     
%           
% 
% 
% % 
% % 
%     f2 = figure()
%     subplot(1,2,1)
%     b1 = boxplot(Entropy, 'labels', {'1', '2', '3', '4', '5'});
%     set(b1,{'linew'},{2})
%     xlabel('Number of Resources')
%     ylabel('Entropy')
%     title('Ecosystem Entropy at Different Number of Resources') 
%     ax = gca;
%     ax.FontSize = 16;
%     box off
%     grid on
%     subplot(1,2,2)    
%     b2 = boxplot(Surviving, 'labels', {'1', '2', '3', '4', '5'});
%     set(b2,{'linew'},{2})
%     xlabel('Number of Resources')
%     ylabel('Number of Surviving Species')
%     title('Ecosystem Biodiversity at Different Number of Resources') 
%     ax = gca;
%     ax.FontSize = 16;
%     box off
%     grid on
%     
    
% cd /end/home/student2/Desktop/as6616/matlab/plankton-simulator %change diretory back to current folder

end
