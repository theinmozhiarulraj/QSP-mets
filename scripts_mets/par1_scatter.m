% patients in rows and parameters in columns
data=params_in.all;

response_status=params_out.RECIST_mets;

[xsize,ysize]=size(data);
gscatter(repelem(1,xsize),data(:,1),response_status)

response_status=table(response_status);
data=table(data);
dataset=[data response_status];
T.Properties.VariableNames = params_in.names;


writetable(dataset,'params_vps.csv')

