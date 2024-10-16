% This script extracts baseline values of different quantities from
% simDataPSAextract and writes in a file

% number of quantities to plot 
nquantities=length(simDataPSAextract);
index = params_out.iPatient;

dataset=table();

for i=1:nquantities
    data=table(simDataPSAextract(i).simData.Data(1,:)');
    dataset=[splitvars(dataset) data];
    dataset.Properties.VariableNames(i)=simDataPSAextract(i).simData.DataNames{1,1};
end


response_status=params_out.RECIST_mets(index);
response_status=table(response_status);
dataset=[splitvars(dataset) response_status];
dataset.Properties.VariableNames(nquantities+1) = "Response";

writetable(dataset,'baseline_biomarkers.csv')

names=dataset.Properties.VariableNames;
fileID = fopen(fullfile(pwd, 'biomarker_names_units.txt'), 'wt');
if fileID == -1
  error('Cannot open the file');
end
%fprintf(fileID,'i NOT A\n');
for row = 1:numel(names)
  fprintf(fileID, '%s\n', names{row});
end
fclose(fileID);
