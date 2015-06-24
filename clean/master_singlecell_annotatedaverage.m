
%% DATA IMPORT : OLD IS GOLD

clear all;
k=4 ;

[A,B] = xlsread('shalek.xlsx');
gnames = strtok(B(2:end,2)); % strtok gets rid of spaces
cnames = strtok(B(1,4:end)); 
dat = A;
ligs = 1:5;
fngs = [6 7 12];
notches = [8:11];

%DO NOT SHAVE OFF THE ZERO COLUMNS HERE , the columns are ccells still; I havent
%transposed the matrix.
%dat=dat(:,any(dat));
%dat(dat<120)=1;
%dat = normr(dat);


%% AVERAGE OVER ANNOTATIONS

    %EXAMPLE
matchStr = regexp(B,'Ifnar1_KO_LPS_4h\w*');
find(~cellfun(@isempty,matchStr));
    %
    

cell_type = {'On_Chip_Stimulation_LPS';'Ifnar1_KO_LPS_4h';...
'IFNB_2h';...
'Unstimulated_Replicate';...
'LPS_1h';...
'LPS_2h';...
'(?m)^LPS_4h_S';...
'LPS_6h';...
'Unstimulated_S';...
'LPS_4h_Replicate';...
'PAM_1h';...
'PAM_2h';...
'PAM_4h';...
'PAM_6h';...
'PIC_1h';...
'PIC_2h';...
'PIC_4h';...
'PIC_6h';...
'Stat1_KO_LPS_4h';...
'Tnfr_KO_LPS_4h' };


%just collect the indices of each cell type in a cell array :(..I hate cell
%arrays but I must use them here :( :( :(
for typ = 1:length(cell_type)
    type_string = regexp(B,strcat(cell_type{typ},'\w*'));
    temp_store = find(~cellfun(@isempty,type_string));
    type_indices_cell{typ} = temp_store;
end



%this loop is complicated. It randomly chooses 20 cells from each of the
%~90 cells belonging to each cell type and then averages their expression values for the notch relevant
%genes. This way, I will generate 200 cells for each of the 20 cell types.
data_point = [];
type_indices_mat = [];
for typ = 1:length(cell_type)
    type_indices_mat = type_indices_cell{typ};
    
    %%select twenty cells for a given cell type and average them to get one
    %%data point. Do this 200 times to get 200 data points.
    for j = 1:200
        
      
        %Randomly select 20 cells for a given cell type
        msize = numel(type_indices_mat);
        randomly_chosen = type_indices_mat(randperm(msize, 20));
        
        %CRAZY STEP : B AND A HAVE Different dimensions. So, I convert the
        %numerical dimension of B to the appropriate column number of A. a
        %week later, I will have no idea how this works. But it works. abs
        %is for the first value which will be -1 but should be 1
        randomly_chosen  = (((randomly_chosen - 1 )/15) + 1)';
        
        %average the 20 cells to get 1 data point outta 200 for this cell type
        data_point = [data_point mean(dat(:,randomly_chosen),2)] ;
              
        %
    end
    
end




%% ASSIGNd THE SAMPLED DATA TO THE DAT VARIABLE
dat = data_point';

%EM AGAIN FOR THE NEWLY SAMPLED DATA
%idx = emgm(dat',k);
%idx = idx';

%k-medoids
[idx,C]  = kmedoids(dat,k);

%CODE FOR SORTING ACCORDING TO ID - (WHICH MIGHT EVEN BE THE STATE!)
datforplot=[idx dat];
dat = sortrows(datforplot);
%dat=datforplot;
idx = dat(:,1);
dat = dat(:,2:end);

%%

if 1%DONT NEED THIS FOR AIC PLOTS
if 1
    dat = normr(dat)

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

%glyphplot(dat,'obslabels',cellstr(num2str(idx))) ;

%glyphplot((dat),'obslabels',cellstr(num2str(idx)),'standardize','off') ;

glyphplot((rep_dat),'obslabels',cellstr(num2str((1:k)')),'standardize','off') ;

%glyphplot(rep_dat,'obslabels',cellstr(num2str((1:k)')),'standardize','off') ;


%glyphplot(dat) 
end %DONT NEED THIS FOR AIC PLOTS



%% AIC 
dat=data_point';
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





