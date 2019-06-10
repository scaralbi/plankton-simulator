
function datagather()
clear all

    cd /end/home/student2/Desktop/as6616/matlab/plankton-simulator/outputs/sym_19.06.04_19.46%change diretory to config_files

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
%     rank_homo_plot = zeros(N, 2, num_files);



    

    %iterate for every config file
    for file = 26%:num_files
        file_name = strcat('config', num2str(file), '_data.mat')
%         data = load(file_name); %import parameter structures
        
        load(file_name, 'params', 'Data', 'global_pool', 'deme_pool', 'resource_pool', 'shuffle_history');
%         Data
        surviving(file) = 0;
        for s = 1:Data{2,3}
            if global_pool(end,s) > 0.0001
                surviving(file) = surviving(file) + 1;
            end
        end
        H(file) = Data{2,11};
        num_resources(file) = params.num_resources;
        

      dims = Data{2,1};
% %            %Rank Abundance Plot
%        rank_plot(:,1,file) = 1 : N;
% % %        rank_homo_plot(:,1,file) = 1 : N;
% % % 
%        tot_pop = sum(global_pool(end,:,end));
% % %        tot_homo_pop = sum(global_pool(400,:,end));
% % % 
%        steady_pop = global_pool(end,:,end);
% % %        homo_pop = global_pool(400,:,end);
% % % 
%        sort_pop = sort(steady_pop, 'descend');
% % %        sort_homo_pop = sort(homo_pop, 'descend');
% % % 
%        rel_pop = sort_pop./tot_pop;
% % %        rel_homo_pop = sort_homo_pop./tot_homo_pop;
% % % 
%        for i = 1: length(rel_pop)
%            rank_plot(i,2,file) = rel_pop(i);
% % %            rank_homo_plot(i,2,file) = rel_homo_pop(i);
%        end
%        rank_plot(:,2,file) = rel_pop';   
%        Rank = 1:Data{2,3};
%        Rank = mat2cell(Rank,1,Data{2,3});
%               
% %       


%       for deme = 1:params.num_demes


%     %%Distance Plot 
%        dist_plot = zeros(params.T, N);
%        invasion_plot = zeros(params.T,2,dims(1)*dims(2));
% 
%         %linear to subsrcipt indexes
%         IND = 1 : dims(1)*dims(2);
%         dims = dims;
%         rows = dims(1);
%         cols = dims(2);
%         [col_index,row_index] = ind2sub(dims,IND);
%         Distance = zeros(rows*cols, rows*cols); %init Distance matrix
% 
%         %Calculate Euclidean Distance between Demes
%         for pre_deme = 1 : rows*cols 
%             x = row_index(pre_deme); %compute x coordinate of donor deme
%             y = col_index(pre_deme); % compute y coordinate of donor deme 
%             for post_deme = 1:rows*cols
%                  X = row_index(post_deme); %compute x coordinate of recipient deme
%                  Y = col_index(post_deme);%compute y coordinate of recipient deme   
%                  Distance(pre_deme, post_deme) = sqrt((X - x)^2 + (Y- y)^2); %compute euclidean Distance (Distance = sqrt(X_diff^2 + Y_diff^2)     
%             end
%         end
%         
%         
%     %calculate distance travelled
%     for T = 1:params.T
%         for row = 1:dims(1)*dims(2) %rows of shuffle history
%             for col = 1:dims(1)*dims(2) %cols of shhuffle history  
%                 if shuffle_history(row,col,T) > 0
%                     species = shuffle_history(row,col,T) ;
%                     invasion_plot(T, 1, row) = species; %add species ibnto donot log in speciefied donor deme = row
%                     invasion_plot(T, 2, col) = species; %same as above buit for recipient
%                     distance = Distance(row,col);
%                     dist_plot(T,species) = dist_plot(T,species) + distance;
%                 end
%             end
%         end
%     end
%     
%     
%             
%      %  plot inavsions
%     f2 = figure('units','normalized','outerposition',[0 0 1 1]);
%     for deme = 1:params.num_demes; hold on
%         subplot(3, 3, deme);    
%         for t = 1:params.T
%             if invasion_plot(t,1,deme) == 0
%                 size_l = 1;
%             else 
%                 leaver_color = invasion_plot(t,1,deme);
%                 size_l = 10;
%             end
%             
%             if invasion_plot(t,2,deme) == 0
%                 size_c = 1;
%             else
%                 size_c = 10;
%             end
%             
%             s1= scatter(t, invasion_plot(t,1,deme), size_l, 'r', 'filled'); hold on
%             s2 = scatter(t, invasion_plot(t,2,deme), size_c, 'b', 'filled');
%             s1.Marker = 'o';
%             s2.Marker = 'o';
%         end
% 
%         xlabel('Time (iterations)');
%         ylabel('Species');
%         ax = gca;
%         ax.XLim = [1 params.T];
%         ax.YLim = [0 27];
%         ax.FontSize = 10;
%         ax.YTick = 1 : 27;
%         xlabel('Time (iterations')
%         ylabel('Species')
%         % Color YTickLabels with colormap
%         cm = hsv(N);
%         for i = 1:N
%             ax.YTickLabel{i} = sprintf('\\color[rgb]{%f,%f,%f}%s', cm(i,:), ax.YTickLabel{i});
%         end
%         ax.XGrid = 'on';
%         ax.YGrid = 'on';
%         ax.GridLineStyle = ':';
%         box off
%         title(['Deme' num2str(deme)]); 
%     end
%     hold off
%             
%     tot_dist = [];
%     for species = 1 : N
%         tot_dist(species) = sum(dist_plot(:,species));
%     end
%     
% 
%        
       %Bar plot with distances
%     figure()
%     b = bar(tot_dist, 'FaceColor','flat')
%     color = hsv(N);
%     for k = 1:size(tot_dist,2)
%         b.CData(k,:) = color(k,:);
%     end
%     ax = gca;
%     ax.XTick = 1:27;
%     

%     figure()
%       for i = 1:27
%         plot(1:200, dist_plot(:,i)); hold on
%       end

%         Global pool plot
%     f1 = figure('units','normalized','outerposition',[0 0 1 1]);
% %     set(f1, 'Visible', 'off');
%     h1 = plot(0:params.T-1, global_pool(:,:), 'LineWidth', 2); hold on
%     set(h1, {'color'}, num2cell(hsv(Data{2,3}),2));
%     xlabel('Time (iterations)');
%     ylabel('Species Abundance');
%     ax = gca;
%     ax.XLim = [0 params.T-1];
%     ylim = max(global_pool(:));
%     ax.YLim = [0 ylim];
%     ax.FontSize = 12;
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
%     title('Simulated Metacommunity Dynamics'); hold off
%     
    
% %      plot demes dynamics
%     f2 = figure('units','normalized','outerposition',[0 0 1 1]);
% %     set(f2, 'Visible', 'off');
%     for deme = 1:params.num_demes; hold on
%         subplot(3, 3, deme);
%         yyaxis left
%         p1 = plot(0:params.T-1, deme_pool(:,:,deme), '-', 'LineWidth', 3);
%         set(p1, {'color'}, num2cell(hsv(Data{2,3}),2));
%         ax1 = gca
%         xlabel('Time (iterations)')
%         ylabel('Species Abundance')
%         ylim = max(deme_pool(:));
%         ax1.YLim = [0 ylim];
% %         ax.FontSize = 14;
% %         ax.YTick = 0 : 5: ylim;
%         yyaxis right
%         p2 = plot(0:params.T-1, resource_pool(:,:,deme), 'LineWidth', 2);
%         ax2 = gca
%         ylabel('Resource Concentration')
%         y2lim = max(resource_pool(:));
%         ax2.YLim = [0 y2lim];
%         p2(1).LineStyle = '--';
%         p2(1).Color = 'r';
%         p2(2).LineStyle = '--';
%         p2(2).Color = 'g';
%         p2(3).LineStyle = '--';
%         p2(3).Color = 'b';
% 
%         
%         %         legend('i = 1', 'i = 2', 'i = 3')
%         ax = gca;
%         ax.XGrid = 'on';
%         ax.YGrid = 'on';
%         ax.GridLineStyle = ':';
%         box off
%         title(['Deme' num2str(deme)]); 
%     end
   
%     hold off
%     
%     plot simplex
    f3 = figure(); hold on
%     figtitle = sgtitle('T = 0');
    axis tight manual % this ensures that getframe() returns a consistent size
    filename = 'simplex.gif';    
    deme_color = summer(params.num_demes); %demes colromap
    species_color = hsv(Data{2,3}); %species colormap
    for iteration = 1:params.T
        for deme= 1:params.num_demes
            hAxis = subplot(3,3,deme);
            su = scatter3(params.supply(1,1,deme)/params.total_supply, params.supply(1,2,deme)/params.total_supply, params.supply(1,3,deme)/params.total_supply, 100, 'd', 'filled', 'k'); hold on
    %         su = scatter3(params.homo_supply(1,1)/params.total_supply, params.homo_supply(1,2)/params.total_supply, params.homo_supply(1,3)/params.total_supply, 100, 'd', 'filled', 'k'); hold on
            f = fill3([1 0 0],[0 1 0],[0 0 1],deme_color(deme,:)); 
            for row = 1:params.size(1)
                for col = 1:params.size(2)
                    for i = 1 : Data{2,3}
                        if deme_pool(iteration,i,deme) > 0.000001
    %                     h = scatter3(params.strategies(1,1,i),params.strategies(1,2,i), params.strategies(1,3,i), 70, 'o', 'filled'); hold on
                        h = scatter3(params.strategies(1,1,i), params.strategies(1,2,i), params.strategies(1,3,i), 100, 'o', 'filled'); hold on
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
            title(['Deme' num2str(deme)]); 
            grid on
            box off
            axis off
            newfigtitle = ['T = ', num2str(i)];
            figtitle.String = newfigtitle;
            drawnow
              frame = getframe(f3); 
              im = frame2im(frame); 
              [imind,cm] = rgb2ind(im,256);
            del = 0;
            % Write to the GIF File 
            if iteration == 1 
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',del);
            else 
              imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',del); 
            end 
        end
    end
   
    hold off
    
    
    
%        

%     %      plot resoure dynamics
%     f3 = figure('units','normalized','outerposition',[0 0 1 1]);
% %     set(f3, 'Visible', 'off');
%     for deme = 1:params.num_demes; hold on
%         subplot(3, 3, deme);
%         h3 = plot(0:params.T-1, resource_pool(:,:,deme), 'LineWidth', 2);
%         xlabel('Time (iterations)');
%         ylabel('Nutrient Conentration');
%         legend('i = 1', 'i = 2', 'i = 3')
%         ax = gca;
%         ax.XLim = [0 params.T-1];
%         ylim = max(resource_pool(:));
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
%     figure()
%     semilogy(rank_plot(:,1,file), rank_plot(:,2,file));hold on
%     xlabel('Species Rank')
%     ylabel('Relative Abundance')
%     grid on
%     ax = gca
%     ax.XLim = [1 27]
%     ax.XTick = 1:27
%     box off
%     title('Rank Abundance Plot with Global Shuffle')
%     hold off
    end
%     
    cd /end/home/student2/Desktop/as6616/matlab/plankton-simulator
    
    
%    
% %     %Rank Abundance plot   
%    std_1 = [];
%    mean_1 = [];
%    
%    std_2 = [];
%    mean_2 = [];
   
%    std_3 = [];
%    mean_3 = [];
   
%    std_4 = [];
%    mean_4 = [];
%    
%    std_5 = [];
%    mean_5 = [];
%    
%    for rank = 1 : N
%            std_1(rank) = std(rank_plot(rank,2,1:20))
%            mean_1(rank) = mean(rank_plot(rank,2,1:20))
%            
%            std_2(rank) = std(rank_plot(rank,2,21:40))
%            mean_2(rank) = mean(rank_plot(rank,2,21:40));
% %            
% %            std_3(rank) = std(rank_plot(rank,2,21:30))
% %            mean_3(rank) = mean(rank_plot(rank,2,21:30));
% %            
% %            std_4(rank) = std(rank_plot(rank,2,61:80))
% %            mean_4(rank) = mean(rank_plot(rank,2,61:80));
%            
% %            
% %            std_5(rank) = std(rank_plot(rank,2,81:100))
% %            mean_5(rank) = mean(rank_plot(rank,2,81:100));
% end
%    
%    mean_global = log(mean_global)
%    mean_local = log(mean_local)
%    std_global = log(std_global)
%    std_local = log(std_local)

%    
% 
%     mean_plot = []
%     std_plot = []
%     for i = 1:length(mean_1)
%         if mean_1(i) > 0.0015
%             mean_plot(i) = mean_1(i)
%             std_plot(i) = std_1(i)
%         end
%     end
% 

%       negbars1 = []
%       negbars2 = []
%       negbars3 = []
%       negbars4 = []
%       negbars5 = []
      
%       
%     filename = 'margata.txt';
%     A = importdata(filename)
%     B = A';
%    
%     rel_ab = B(1:3:length(B));
%     
%     Tot = sum(rel_ab)
%     
%     for i = 1:length(rel_ab)
%         rel_ab(i) = rel_ab(i)/Tot;
%     end
%     

%     Data
% 
%    figure()
%    e1 = errorbar(mean_1,std_1, '-o','MarkerSize',7, 'MarkerEdgeColor','blue','MarkerFaceColor','blue', 'LineWidth', 4);hold on
%    for i = 1 : length(std_1)
%        negbars1(i) = std_1(i);
%        if mean_1(i) - negbars1(i) <=0
%            negbars1(i) = mean_1(i) - 10^-6;
%            e1.CapSize = 0;
%        end
%    end
%    e1.YNegativeDelta = negbars1;
%    e1.Bar.LineStyle = 'dotted';
%    e1.Bar.LineWidth = 2.5;
% %    
% % %    s1 = semilogy(1:length(rel_ab),rel_ab, 'r', 'LineWidth', 3);hold on
% % %    s1.LineStyle = '--'
% % % 
%    e2 = errorbar(mean_2,std_2, '-s','MarkerSize',7, 'MarkerEdgeColor','red','MarkerFaceColor','red', 'LineWidth', 4);
%    for i = 1 : length(std_2)
%       negbars2(i) = std_2(i);
%        if mean_2(i) - negbars2(i) <=0
%            negbars2(i) = mean_2(i) - 10^-6;
%        end
%    end
%    e2.YNegativeDelta = negbars2;
%    e2.Bar.LineStyle = 'dotted';
%    e2.Bar.LineWidth = 2.5;
%    e3 = errorbar(mean_3,std_3, '-d','MarkerSize',7, 'MarkerEdgeColor','g','MarkerFaceColor','g', 'LineWidth', 3, 'CapSize',0);
%    e3.Color = 'g';
%    for i = 1 : length(std_3)
%       negbars3(i) = std_3(i);
%        if mean_3(i) - negbars3(i) <=0
%            negbars3(i) = mean_3(i) - 10^-6;
%        end
%    end
%    e3.YNegativeDelta = negbars3;
%    e3.Bar.LineStyle = 'dotted';
%    e3.Bar.LineWidth = 2;
%    e4 = errorbar(mean_4,std_4, '-p','MarkerSize',7, 'MarkerEdgeColor','m','MarkerFaceColor','m', 'LineWidth', 3, 'CapSize',0);
%    e4.Color = 'm';
%    for i = 1 : length(std_4)
%       negbars4(i) = std_4(i);
%        if mean_4(i) - negbars4(i) <=0
%            negbars4(i) = mean_4(i) - 10^-6;
%        end
%    end
%    e4.YNegativeDelta = negbars4;
%    e4.Bar.LineStyle = 'dotted';
%    e4.Bar.LineWidth = 2;
%    e5 = errorbar(mean_5,std_5, '-h','MarkerSize',7, 'MarkerEdgeColor',[255/255,165/255,0],'MarkerFaceColor',[255/255,165/255,0], 'LineWidth', 3, 'CapSize',0);
%    e5.Color = [255/255,165/255,0]; 
%    for i = 1 : length(std_5)
%       negbars5(i) = std_5(i);
%        if mean_5(i) - negbars5(i) <=0
%            negbars5(i) = mean_5(i) - 10^-6;
%        end
%    end
%    e5.YNegativeDelta = negbars5;
%    e5.Bar.LineStyle = 'dotted';
% %    e5.Bar.LineWidth = 2;
%    legend('reference simulation', 'nonfixed enzyme budget');
%    title('Simulated Rank Abundance')
%    xlabel('Species Rank')
%    ylabel('Relative Mean Abundance (n_{\sigma}/N)')
%    ax1 = gca
%    ax1.YScale = 'log';
%    ax1.FontSize = 18;
%    ax1.XLim = [1 11];
%    ax1.XTick = 1:11;
%    ax1.YLim = [10^-5 1];
%    box off
%    grid on
%    hold off
% % %    
% l = length(mean_3)
% % 
% A = [rel_ab(1:45); mean_3(1:45)]
% 
% R = corrcoef(rel_ab(1:45), mean_3(1:45))
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

% 
% % % 
% Entropy = zeros(20,5);
% Surviving = zeros(20,5);

%     for col = 1:cases
%         Entropy(:,col) = H(col*replicates - replicates +1 : col*replicates);
%     end
%     
% %     
%     Entropy(:,1) = H(1:20);
%     Entropy(:,2) = H(21:40);
%     Entropy(:,3) = H(41:60);
%     Entropy(:,4) = H(61:80);
%     Entropy(:,5) = H(81:100);
% %     
% %     
%     Surviving(:,1) = surviving(1:20);
%     Surviving(:,2) = surviving(21:40);
%     Surviving(:,3) = surviving(41:60);
%     Surviving(:,4) = surviving(61:80);
%     Surviving(:,5) = surviving(81:100);
%     
    
%     for col = 1:cases
%         Surviving(:,col) = surviving(col*replicates - replicates +1 : col*replicates);
%     end
    

%      fit_richness = [];
%      for a = 1: 25
%          fit_richness(a) = 5 * a^0.8;
%      end
%          
%   
% % % 

%     f2 = figure()
% %     subplot(1,2,1); hold on
%     X = [1 2]
% %     b1 = boxplot(Entropy, 'positions', X, 'labels', X );hold on
% %     set(b1,{'linew'},{2})
% %     xlabel('Number of Resources')
% %     ylabel('Entropy')
% %     title('Simulated Biodiversity at Different Number of Resources') 
% %     ax = gca
% %     ax.FontSize = 18;
% %     grid on
% %     box off
% %     subplot(1,2,2)  
% 
%     b2 = boxplot(Surviving, 'position', X, 'labels', {'reference simulation', 'non-fixed enzyme budget'} ); hold on
% %     p1 = plot(1:25, fit_richness, '--')
%     set(b2,{'linew'},{2})
% %     xlabel('Number of Resources')
%     ylabel('Species Richness')
%     title('Species Richness - Reference Simulation vs Non-fixed Enzyme Budget') 
%     ax = gca;
%     ax.FontSize = 16;
%     box off
%     grid on; hold off
% %     
% cd /end/home/student2/Desktop/as6616/matlab/plankton-simulator %change diretory back to current folder

end
