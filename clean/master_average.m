%SUPERPOSE DATA
%% Data loading and normalization
% note that immgen data appear to be standardized to have same
% distributions across cell lines, and are in "linear" expression units,
% and have > ~120 is 'on', <45 is 'off', according to: 
% https://www.immgen.org/Protocols/ImmGen%20QC%20Documentation_ALL-DataGeneration_0612.pdf
clear all;
k=4 ;
 

[A,B] = xlsread('notchimmgen.xlsx');
gnames = strtok(B(2:end,2)); % strtok gets rid of spaces
cnames = strtok(B(1,4:end)); 
dat = A(1:end,4:end);
ligs = 1:5;
fngs = [6 7 12];
notches = [8:11];

%SHAVE OFF THE ZERO COLUMNS 
dat=dat(:,any(dat));

dat(dat<120)=1;

average_dat=dat;

if 0 % averaging for notchimmgen data
average_dat(1,:) = (dat(1,:) + dat(4,:)) /2 ;
average_dat(2,:) = (dat(2,:) + dat(3,:) + dat(5,:)) /3 ;
average_dat(3,:) = (dat(6,:)) /1 ;
average_dat(4,:) = (dat(8,:) + dat(9,:) + dat(10,:) + dat(11,:)) /4 ;
average_dat(5,:) = (dat(7,:) + dat(12,:)) /2 ;

end


%

%dat(dat<120)=1;


if 0 % averaging for shalek data
average_dat(1,:) = (dat(7,:) + dat(13,:)) /2 ;
average_dat(2,:) = (dat(1,:) + dat(10,:) + dat(6,:)) /3 ;
average_dat(3,:) = (dat(11,:)) /1 ;
average_dat(4,:) = (dat(2,:) + dat(3,:) + dat(5,:) + dat(8,:)) /4 ;
average_dat(5,:) = (dat(14,:) + dat(11,:)) /2 ;

end

dat  = dat' ;
dat = normr(dat);
strong_dat=dat;

average_dat = average_dat';
average_dat = normr(average_dat);
dat=average_dat;





%EM
idx = emgm(strong_dat',k);
idx = idx';
%idx = kmeans(dat,k);

%

% PLOT 

%parallelcoords(dat)

%andrewsplot(dat,'group',idx)

%CODE FOR SORTING ACCORDING TO ID - (WHICH MIGHT EVEN BE THE STATE!)
datforplot=[idx dat];
dat = sortrows(datforplot);
%dat=datforplot;
idx = dat(:,1);
dat = dat(:,2:end);

%dat=dat(:,any(~dat)); 

if 1

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
end



%%Incorporating Changes in my representation as suggested by Michael

if 0 %%A Little LESS SPACE

dat = zeros(size(dat,1),24);

dat(:,1) = strong_dat(:,4);
dat(:,2) = strong_dat(:,1);

dat(:,6) = strong_dat(:,2);
dat(:,7) = strong_dat(:,3);
dat(:,8) = strong_dat(:,5);

dat(:,12) = strong_dat(:,10);
dat(:,13) = strong_dat(:,11);
dat(:,14) = strong_dat(:,9);
dat(:,15) = strong_dat(:,8);

dat(:,19) = strong_dat(:,6);
dat(:,20) = strong_dat(:,12);
dat(:,21) = strong_dat(:,7);
end




if 0 %%LOTS OF SPACE

dat = zeros(size(dat,1),36);

dat(:,1) = strong_dat(:,4);
dat(:,2) = strong_dat(:,1);

dat(:,9) = strong_dat(:,2);
dat(:,10) = strong_dat(:,3);
dat(:,11) = strong_dat(:,5);

dat(:,18) = strong_dat(:,10);
dat(:,19) = strong_dat(:,11);
dat(:,20) = strong_dat(:,9);
dat(:,21) = strong_dat(:,8);

dat(:,28) = strong_dat(:,6);
dat(:,29) = strong_dat(:,12);
dat(:,30) = strong_dat(:,7);
end
%THIS IS WORKING !

%glyphplot(dat,'obslabels',cellstr(num2str(idx))) ;

%glyphplot((dat),'standardize','off') ;

%glyphplot((rep_dat),'obslabels',cellstr(num2str((1:k)')),'standardize','off') ;

%glyphplot(rep_dat,'obslabels',cellstr(num2str((1:k)')),'standardize','off') ;


%glyphplot(dat) 






%% AIC 
dat=strong_dat;
AIC = zeros(1,20);
GMModels = cell(1,20);
options = statset('MaxIter',1000);
for k = 1:20
    GMModels{k} = fitgmdist(dat,k,'Options',options,'CovarianceType','full','RegularizationValue',0.01);
    AIC(k)= GMModels{k}.AIC;
end

[minAIC,numComponents] = min(AIC);
numComponents

AIC(AIC==0) = nan;
plot(1:size(AIC,2),AIC)
xlabel('k-The number of clusters')
ylabel('AIC')
title('The value of AIC for different number of clusters')



%% K Means Clustering

idx = kmeans(dat,4);



%% OPTICS Clustering


[RD,CD,order]=optics(dat,50)


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




