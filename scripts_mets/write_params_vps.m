% This script writes parameter values of all virtual patients and response
% status in a file

% patients in rows and parameters in columns
data=params_in.all;
data=table(data);

response_status=params_out.RECIST_mets;
response_status=table(response_status);

% [len,~]=size(params_in.names);
% for i=1:len
%     data.Properties.VariableNames(i) = params_in.names(i);
% end

dataset=[splitvars(data) response_status];
dataset.Properties.VariableNames(:) = [params_in.names(:);"Response"];

writetable(dataset,'params_vps.csv')

