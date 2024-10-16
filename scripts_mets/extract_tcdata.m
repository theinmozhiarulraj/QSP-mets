% thein: This script extracts timecourse data from simulation results and
% also calculates some derived quantities for biomarker analysis

% total number of patients
n_PSA = length(params_out.iPatient);

% index of patients
index = params_out.iPatient;

% initialize the structure
simDataPSAextract = struct;

% number of timepoints
ntimepoints=length(simDataPSA(index(1)).simData.Time);

%sfiltered=simDataPSA;

%ncols=size(simDataPSA_filtered(index(i)).simData.Data,2);

% % set all -ve values to 0
% for i=1:length(index)
%     sfiltered(index(i)).simData.Data(simDataPSA(index(i)).simData.Data<0)=0;
% end

% for i=1:length(index)
%     for j=1:ncols
%         simDataPSA_filtered(index(i)).simData.Data(simDataPSA(index(i)).simData.Data(:,j)<0)=0;
%     end
% end

% extract quantities of central or peripheral compartments
simDataPSAextract=extract_tcdata_central_periphery(n_PSA,n_T_specs,ntimepoints,index,simDataPSA,simDataPSAextract);

% quantities from tumor compartments - average of all tumors
simDataPSAextract=extract_tcdata_tumor(n_PSA,params_in,n_T_specs,cancer_clones,TumorArray,ntimepoints,index,simDataPSA,simDataPSApost,simDataPSAextract);

% quantities from synapse compartments such as PDL1 expressions
simDataPSAextract=extract_tcdata_synapse(n_PSA,n_T_specs,TumorArray,LNArray,ntimepoints,index,simDataPSA,simDataPSApost,simDataPSAextract);

% quantities from LN compartments
simDataPSAextract=extract_tcdata_lymphnodes(n_PSA,n_T_specs,TumorArray,LNArray,ntimepoints,index,simDataPSA,simDataPSApost,simDataPSAextract);

%save(simDataPSAextract,'simDataPSAextract_treatment.mat')

save('simDataPSAextract_treatment.mat','simDataPSAextract')

nquantities=length(simDataPSAextract);

for nq=1:nquantities

    test_real=isreal(simDataPSAextract(nq).simData.Data);
    if(test_real==0)
        disp(nq)
    end

%     test_nan=isnan(simDataPSAextract(nq).simData.Data);
%     if(test_nan==1)
%         disp(nq)
%     end
% 
%     test_inf=isinf(simDataPSAextract(nq).simData.Data);
%     if(test_inf==1)
%         disp(nq)
%     end
end

