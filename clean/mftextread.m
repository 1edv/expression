function [colHeaders,rowHeaders,mydata] = mftextread(fileName)
%  Syntax to be used: [columnHeaders, rowHeaders, numericData] =  mytextread('filename.txt')
%  where 
%      columnHeaders will be a cell array
%      rowHeaders will be a cell array
%      numericData will be a single-precision array
% 
%  The contents of the input text-file should be formatted as follows:
% 
%      Title A  B  C
%        X  1  2  3
%        Y  3  4  5
%        Z  4  5  6
% 
%  Manu Raghavan,
%  June 25, 2008

% Open the file
fid=fopen(fileName);

% Start reading the data
tline = fgetl(fid); % Read in the first line only (text headers)
% tline = tline(tline~=' '); %Get rid of any spaces
tabLocs=findstr(char(9),tline); % find the tabs
colHeaders=cell(1,(length(tabLocs)+1));
start=1;

for col=1:length(tabLocs)
    colHeaders{col}=tline(start:tabLocs(col)-1);
    start=tabLocs(col)+1;
end
colHeaders{col+1}=tline(start:end);

row = 1;
tline = fgetl(fid) % Get second row (first row of data)
while(1)
    tabLocs=findstr(char(9),tline); % find the tabs
    c = textscan(tline,'%s%f32%f32%f32%f32%f32%f32%f32%f32%f32%f32%f32%f32%f32%f32%f32%f32%f32%f32%f32%f32%f32%f32','Delimiter',char(9))
    rowHeaders{row} = c{1}; % Get column header
    for i=2:length(c)-1
        mydata(row,i-1) = c{i}; % Get numeric data
    end
    tline = fgetl(fid); % Go to next line in text file
    if(length(tline)==1)
        if(tline==-1) % Reached end of file, terminate
            break
        end
    else
        row = row+1;
    end        
end