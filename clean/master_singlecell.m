%
%% Data loading and normalization
% note that immgen data appear to be standardized to have same
% distributions across cell lines, and are in "linear" expression units,
% and have > ~120 is 'on', <45 is 'off', according to: 
% https://www.immgen.org/Protocols/ImmGen%20QC%20Documentation_ALL-DataGeneration_0612.pdf
clear all;
k=4 ;


[A,B] = xlsread('shalek.xlsx');
gnames = strtok(B(2:end,2)); % strtok gets rid of spaces
cnames = strtok(B(1,4:end)); 
dat = A(1:end,4:end);



%SHAVE OFF THE ZERO COLUMNS 
dat=dat(:,any(dat));
dat = zscore(dat); %Uncomment If you want to use zscores instead..
dat  = dat' ;
strong_dat_actual = dat;

%dat = normr(dat); %Comment out if you want to use zscores
strong_dat=dat;




%%%%%%%%%%
%ClUSTERING SECTION
%EM
%idx = emgm(dat',k);
%idx = idx';

%%%%%%% t-SNE 
%score = tsne(dat, [], 3, 12,30);
%dat = score;
%%%%%%%


%%%%%% PCA
%[pc,score,latent,tsquare] = princomp(dat);
%dat = score ;
%%%%%%%


%K Medoids
[idx,C] = kmedoids(dat,k);



%%%%%%%%%% SORT DATAPOINTS ACCORDING TO CLUSTERS
dat = strong_dat;
%CODE FOR SORTING ACCORDING TO ID -
datforplot=[idx dat];
dat = sortrows(datforplot);
idx = dat(:,1);
dat = dat(:,2:end);
strong_dat_sorted = dat; 
%%%%%%%%%




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Heat Maps
%hm=HeatMap(dat','Colormap','summer','RowLabels',B(2:13,2),'ColumnLabels',yStrings');
%plot(hm)

% Rearrange gene labels for plotting
if 0
temp_dat = dat ;
dat(:,1) = temp_dat(:,10) ;
dat(:,2) = temp_dat(:,11) ;
dat(:,3) = temp_dat(:,8) ;
dat(:,4) = temp_dat(:,9) ;
dat(:,5) = temp_dat(:,6) ;
dat(:,6) = temp_dat(:,7) ;
dat(:,7) = temp_dat(:,12) ;
dat(:,8) = temp_dat(:,2) ;
dat(:,9) = temp_dat(:,3) ;
dat(:,10) = temp_dat(:,5) ;
dat(:,11) = temp_dat(:,1) ;
dat(:,12) = temp_dat(:,4) ;

temp_B = B;
B(2,:) = temp_B(11,:) ;
B(3,:) = temp_B(12,:) ;
B(4,:) = temp_B(9,:) ;
B(5,:) = temp_B(10,:) ;
B(6,:) = temp_B(7,:) ;
B(7,:) = temp_B(8,:) ;
B(8,:) = temp_B(13,:) ;
B(9,:) = temp_B(3,:) ;
B(10,:) = temp_B(4,:) ;
B(11,:) = temp_B(6,:) ;
B(12,:) = temp_B(2,:) ;
B(13,:) = temp_B(5,:) ;
end
%


imagesc(dat')
colormap('summer')
cb = colorbar
title(cb,'Z-Scores')
xlabel('Cells sorted according to the cluster assigned')
title('Heat Map - 1700 DC Cells - KMedoids Clustering')
ax = gca 

ax.YTick = 1:14
ax.YTickLabel = B(2:15,1)'
ax.XTick = [];



%for finding the change point
indices_change_point = [];
for i = 1:k
    indices_change_point = [ indices_change_point find(idx==i,1)];
end
indices_change_point(i+1) = size(dat,1)+1; % for rep_dat loop to be easy




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% My 3D scatter representation

if 1
%%%%% TAKE 2
[pc,score,latent,tsquare] = princomp(dat);
gscatter(score(:,1), score(:,2), idx);

%scatter3(score(:,1), score(:,2),score(:,3),100*ones(1642,1), idx)
rotate3d on
colormap(lines(4))
lcolorbar({'1','2','3','4'},'fontweight','bold')
title('PCA then Cluster - Visualisation of Clusters');
xlabel('PC 1');
ylabel('PC 2');
zlabel('PC 3');
%%%%%
end



if 0

%%%%% TAKE 3
score = tsne(dat, idx, 2, 14,1200);
%score = fast_tsne(dat, 2 , 14, 50,0.5)
%scatter3(score(:,1), score(:,2),score(:,3),10*ones(1642,1), idx)
%score = tsne(score, [], 3, 3, 30);
gscatter(score(:,1), score(:,2), idx);
title('(Cluster actual data) -> (t-sne Dimensionality Reduction) -> (Visualisation)');
xlabel('t-SNE Component 1');
ylabel('t-SNE Component 2');
zlabel('Component 3');
%dat = score;

%%%%%
end


























%%
dat = C;
%%Incorporating Changes in my representation as suggested by Michael
if 1
biased_spaced_dat=dat;

dat = zeros(size(dat,1),36);

dat(:,1) = biased_spaced_dat(:,4);
dat(:,2) = biased_spaced_dat(:,1);

dat(:,9) = biased_spaced_dat(:,2);
dat(:,10) = biased_spaced_dat(:,3);
dat(:,11) = biased_spaced_dat(:,5);

dat(:,18) = biased_spaced_dat(:,10);
dat(:,19) = biased_spaced_dat(:,11);
dat(:,20) = biased_spaced_dat(:,9);
dat(:,21) = biased_spaced_dat(:,8);

dat(:,28) = biased_spaced_dat(:,6);
dat(:,29) = biased_spaced_dat(:,12);
dat(:,30) = biased_spaced_dat(:,7);
end



if 0% averaging for notchimmgen data
%dat = strong_dat_sorted' ;
average_dat(1,:) = (dat(1,:) + dat(4,:)) /2 ;
average_dat(2,:) = (dat(2,:) + dat(3,:) + dat(5,:)) /3 ;
average_dat(3,:) = (dat(6,:)) /1 ;
average_dat(4,:) = (dat(8,:) + dat(9,:) + dat(10,:) + dat(11,:)) /4 ;
average_dat(5,:) = (dat(7,:) + dat(12,:)) /2 ;
dat = average_dat';
end



if 0
% GENERATE REPRESENTATIVE FIGURE FOR EACH CLUSTER
indices_change_point = [];
for i = 1:k
    indices_change_point = [ indices_change_point find(idx==i,1)];
end
indices_change_point(i+1) = size(dat,1)+1; % for rep_dat loop to be easy

%
rep_dat = zeros(k, size(dat,2));
for j = 1:k
    rep_dat(j,:) = mean(dat(indices_change_point(j):indices_change_point(j+1)-1,:));
end
% GLYPH PLOT
dat = rep_dat';

end



%%%%%%%%%%%%%%%%%%%%%%% GLYPH PLOT
%PLOT GLYPHPLOT
%glyphplot(dat,'obslabels',cellstr(num2str(idx))) ;

%glyphplot((dat),'obslabels',cellstr(num2str(idx)),'standardize','matrix') ;

%glyphplot((rep_dat),'obslabels',cellstr(num2str((1:k)')),'standardize','off') ;

%glyphplot((dat),'obslabels',cellstr(num2str((1:k)')),'standardize','off') ;

%glyphplot(rep_dat,'obslabels',cellstr(num2str((1:k)')),'standardize','off') ;

%%%%%%%%%%%%%%%%%%%%%%%%










%% EVAL Cluster tests
for j = 1:10
clust = zeros(size(strong_dat,1),10);
for i=1:10
clust(:,i) = kmedoids(strong_dat,i);
end
eva = evalclusters(strong_dat,clust,'silhouette') ;

CriterionValues(j,:) = eva.CriterionValues ; 
OptimalK(j) = eva.OptimalK ; 
end

%%

%%%%%%%%%%%% PLOT THE RESULTS for Evalclusters
figure

legends = {'k=1','k=2','k=3','k=4','k=5','k=6','k=7','k=8','k=9','k=10','Optimal Eval Criteria'};

subplot(2,1,1)
y= CriterionValues(1:10,:);
x= repmat(1:10,10,1)' ;
plot(x,y, '--o');   % plot with 'DisplayName' property

%hold on 
%plot(x(:,4),y(:,4), '--o' , 'LineWidth' , 2 , 'color','black');   % plot with 'DisplayName' property
title('Cluster Evaluation Criterion is stable over multiple repeated experiments')
legend(legends,'Location','southoutside','Orientation','horizontal') ;


xlabel('Repeated Experiment Number');
ylabel('Cluster Evaluation Criterion Value');
gca.Xtick = 1:1:10;
%set(gca,'XGrid','on')

subplot(2,1,2)
scatter(1:10 , OptimalK(1:10) , 'x' ,'LineWidth',2,'markeredgecolor','black' )
ylim([1 10])
title('Optimal Values of k(number of clusters) is stable over multiple repeated experiments')
xlabel('Repeated Experiment Number');
ylabel('Optimal Number of Clusters');
legend('Optimal Number of Clusters for given repeated experiment','Location','southoutside','Orientation','horizontal') ;

%%%%%%%%%%%%







%% SHAKE 

shaken_dat = shake(dat,1);

AIC = zeros(1,20);
GMModels = cell(1,20);
options = statset('MaxIter',1000);
for k = 1:20
    GMModels{k} = fitgmdist(dat,k,'Options',options,'CovarianceType','full','RegularizationValue',0.001);
    AIC(k)= GMModels{k}.AIC;
end
[minAIC,numComponents] = min(AIC);
AIC(AIC==0) = nan;


shaken_AIC = zeros(1,20);
shaken_GMModels = cell(1,20);
options = statset('MaxIter',1000);
for k = 1:20
    shaken_GMModels{k} = fitgmdist(shaken_dat,k,'Options',options,'CovarianceType','full','RegularizationValue',0.001);
    shaken_AIC(k)= shaken_GMModels{k}.AIC;
end

[shaken_minAIC,shaken_numComponents] = min(shaken_AIC);
shaken_AIC(shaken_AIC==0) = nan;


figure
xlabel('k-The number of clusters')
ylabel('AIC')
title('1700DC cells')
hold on

plot(1:size(AIC,2),AIC,'-',1:size(shaken_AIC,2),shaken_AIC,'--')
legend('real data)','scrambled data')

























































%% OPTICS Clustering


[RD,CD,order]=optics(dat,50)


%% DBSCAN Clustering

%[labs labscore] = dbscan(dat,1,1)
%%




