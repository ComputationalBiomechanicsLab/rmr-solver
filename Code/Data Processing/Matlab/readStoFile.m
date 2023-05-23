function [labels, data] = readStoFile(filename)

% Read in motion file
fid = fopen(filename);

disp(['Loading file...' filename] );

name = fgetl(fid);
versionLine = fgetl(fid);
version = str2num(versionLine(length('version= '):length(versionLine)));
rowLine = fgetl(fid);
row = str2num(rowLine(length('nRows= '):length(rowLine)));
colLine = fgetl(fid);
col = str2num(colLine(length('nColumns= '):length(colLine)));

line = ' ';
while ~strcmp(line,'endheader')
    line = fgetl(fid);
end

line = fgetl(fid);
while length(line) == 0;
    line = fgetl(fid);
end
token = ' ';
remainder = line;
count = 1;
labels = {};
while length(remainder)~=0
%while ~strcmp(token,remainder)
    [token,remainder] = strtok(remainder);
    labels(count) = cellstr(char(token));
    count = count + 1;
end

if length(labels) > col
    labels = labels(1:length(labels)-1); %Get rid of extra blank cell
elseif length(labels) < col
    error('Problem reading in motion file headers');
else
    labels = labels;
end

data = zeros(row, col);
for i=1:row
    line = fgetl(fid);
    temp = sscanf(line,'%f');
    if(contains(line, 'nan'))
        temp(2:col) = NaN;
    end
    data(i,:) = temp';
end

fclose(fid);

return
