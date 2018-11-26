function convert_database(folder,filename,headers)
% function convert_database(folder,filename,headers)
% This function loads all files from the folder and saves them in the
% database format in a csv with headers a cell array of strings
% 
% For the date this always saves +01:00 i.e. German time, as I cannot
% access the timezone from matlab
%
% Also we assume that the date is the last column of the data
% 
% Finally this adds two columns: 
%  subject -> first column, Name of the Subject extracted from the filename
%  session -> increasing session number for each subject

headers = [{'subject','session'},headers,{'datetime'}];

files = dir(folder);
files = files(~[files.isdir]);

% sort fils in order of creation
[~,idx] = sort([files.datenum]);
files = files(idx);

f = fopen(filename,'wt');

fprintf(f,'%s;', headers{1:(end-1)});
fprintf(f,'%s\n', headers{end});
subs = cell(0);
split = strsplit(files(1).name,'_');
subj = split{1};
subs{1,1} = subj;
subs{1,2} = 0;
for iFile = 1:length(files)
    data = load(fullfile(folder,files(iFile).name));
    split = strsplit(files(iFile).name,'_');
    subj = split{1};
    subfit = strcmp(subj,subs(:,1));
    if any(subfit)
        subidx = find(strcmp(subj,subs(:,1)));
        subs{subidx,2} = subs{subidx,2} + 1;
        session = subs{subidx,2};
    else
        subs{size(subs,1)+1,1} = subj;
        subs{size(subs,1)+1,2} = 1;
        session = 1;
    end
    for iLine = 1:size(data.results,1)
        fprintf(f,'%s;',subj);
        fprintf(f,'%d;',session);
        fprintf(f,'%d;', data.results(iLine,1:(end-1)));
        fprintf(f,'%s;',files(iFile).name);
        fprintf(f,'%s\n', datestr(data.results(iLine,end),'YYYY-mm-DDTHH:MM:ss+01:00'));
    end
end

fclose(f);