%
%% Data loading and normalization
% note that immgen data appear to be standardized to have same
% distributions across cell lines, and are in "linear" expression units,
% and have > ~120 is 'on', <45 is 'off', according to: 
% https://www.immgen.org/Protocols/ImmGen%20QC%20Documentation_ALL-DataGeneration_0612.pdf

%% IF YOU WANT TO IMPORT ALL THE GENE DATA
%b = importdata('shalek.txt','\t',1);
%dat = b.data;



%% OLD IS GOLD

clear all;
k=4 ;

[A,B] = xlsread('shalek.xlsx');
gnames = strtok(B(2:end,2)); % strtok gets rid of spaces
cnames = strtok(B(1,4:end)); 
dat = A(1:end,4:end);
ligs = 1:5;
fngs = [6 7 12];
notches = [8:11];


%%
%SHAVE OFF THE ZERO COLUMNS , the columns are ccells still; I havent
%transposed the matrix.
dat=dat(:,any(dat));

%dat(dat<120)=1;


dat  = dat' ;
strong_dat_actual = dat;
dat = normr(dat);
strong_dat=dat;




%CODE FOR SAMPLING AND AVERAGING. CLUSTER OF CLUSTERS IS BEING ATTeMPTED

sample_size = 1000;
sampled_data = zeros(sample_size,size(strong_dat,2));
for q = 1:sample_size 
    q
dat = strong_dat;



%EM
idx = emgm(dat',k);
idx = idx';
%idx = kmeans(dat,k);

%AVERAGE CELLS IN THESE CLUSTERS.


%CODE FOR SORTING ACCORDING TO ID - (WHICH MIGHT EVEN BE THE STATE!)
datforplot=[idx dat];
dat = sortrows(datforplot);
%dat=datforplot;
idx = dat(:,1);
dat = dat(:,2:end);




%Finding the mean
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



sampled_data((4*q-3):4*q,:) = rep_dat;

end%end sampling

%%












%%Incorporating Changes in my representation as suggested by Michael
if 0
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

%ASSIGNE THE SAMPLED DATA TO THE DAT VARIABLE
dat = sampled_data;

%EM AGAIN FOR THE NEWLY SAMPLED DATA
idx = emgm(dat',k);
idx = idx';


%CODE FOR SORTING ACCORDING TO ID - (WHICH MIGHT EVEN BE THE STATE!)
datforplot=[idx dat];
dat = sortrows(datforplot);
%dat=datforplot;
idx = dat(:,1);
dat = dat(:,2:end);


if 1
% FOR GENERATING THE REPRESENTATION
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


if 1 % averaging for shalek data
average_dat(1,:) = (dat(7,:) + dat(13,:)) /2 ;
average_dat(2,:) = (dat(1,:) + dat(10,:) + dat(6,:)) /3 ;
average_dat(3,:) = (dat(11,:)) /1 ;
average_dat(4,:) = (dat(2,:) + dat(3,:) + dat(5,:) + dat(8,:)) /4 ;
average_dat(5,:) = (dat(14,:) + dat(11,:)) /2 ;
rep_dat = average_dat';
end





%THIS IS WORKING !

%glyphplot(dat,'obslabels',cellstr(num2str(idx))) ;

%glyphplot((dat),'obslabels',cellstr(num2str(idx)),'standardize','off') ;


glyphplot((rep_dat),'obslabels',cellstr(num2str((1:k)')),'standardize','off') ;

%glyphplot(rep_dat,'obslabels',cellstr(num2str((1:k)')),'standardize','off') ;


%glyphplot(dat) 






%% AIC 
dat=sampled_data;
shaken_dat = shake(dat,2);

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


%[RD,CD,order]=optics(dat,50)


%% DBSCAN Clustering

%[labs labscore] = dbscan(dat,1,1)
%%





%% SCATTER PLOTS

scatter (idx',dat(:,10),'r.')
hold on;
scatter (idx',dat(:,2),'go')
hold on;
scatter (idx',dat(:,5),'bx')
hold on;
scatter (idx',dat(:,12),'bo')
hold on;
scatter (idx',dat(:,4),'rx')


%%
%% EM
obj = fitgmdist(dat,2);



mindat = min(min(dat))-1;

[Jag2,Dll1,Dll4,Jag1,Dll3,Rfng,Mfng,Notch4,Notch3,Notch1,Notch2,Lfng]=deal(1,2,3,4,5,6,7,8,9,10,11,12);

cgo_all = clustergram(dat,'Standardize','column','RowLabels',gnames,'DisplayRange',1.5,'Symmetric','false','Colormap',redbluecmap);

cgo_all = clustergram(log(dat),'Standardize','column','RowLabels',gnames,'DisplayRange',1.5,'Symmetric','false','Colormap',redbluecmap);

cgo_receptors = clustergram(dat(notches,:),'Standardize','column','RowLabels',gnames(notches),'DisplayRange',1.5,'Symmetric','false','Colormap',redbluecmap);
cgo_ligs = clustergram(log(dat(ligs,:)-mindat),'Standardize','column','RowLabels',gnames(ligs),'DisplayRange',1.5,'Symmetric','false','Colormap',redbluecmap);

%% Plot the relationship between Notch1 and Notch2.

figure(101);
loglog(dat(Notch1,:),dat(Notch2,:),'.',dat(Notch1,:),dat(Notch3,:),'o',dat(Notch1,:),dat(Notch4,:),'.');
a = axis;
hl(1)=line([120 120],[a(3) a(4)]);
hl(2)=line([a(1) a(2)],[120 120]);

xlabel('Notch1');
legend('Notch2','Notch3','Notch4');


% now let's do histograms of each notch:

bins = logspace(1.5,4,40);
figure(102);
colors = 'bgrmcyk';
for i = 1:4
    [yn{i},xn{i}]=hist(dat(notches(i),:),bins);
    semilogx(xn{i},yn{i},colors(i));hold on
end
legend(gnames(notches))
hold off;

% now let's do histograms of fringes:

figure(103);
for i = 1:3
    [yn{i},xn{i}]=hist(dat(fngs(i),:),bins);
    semilogx(xn{i},yn{i},colors(i));hold on
end
legend(gnames(fngs))
hold off;

figure(104);
for i = 1:5
    [yn{i},xn{i}]=hist(dat(ligs(i),:),bins);
    semilogx(xn{i},yn{i},colors(i));hold on
end
legend(gnames(ligs))
hold off;




%%




