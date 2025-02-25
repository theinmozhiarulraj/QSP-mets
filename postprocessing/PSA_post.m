% Function to postprocess the batch simulation data and extract addtional
% outputs of interest like clonality and precentage changes and adding them
% to individual simulations as parameters
%
% Inputs: simDataPSA   -- Object containing the whole simbiology model
%                         outputs for all batch simulations
%         params_in    -- object containing model inputs to be organized
%                         for future sensitivity analysis
%         params_out   -- object containing model outputs to be organized
%                         for future sensitivity analysis
%
% Outputs: simDataPSApost  -- Object containing postprocessed outputs


function simDataPSApost = PSA_post(simDataPSA,params_in,params_out,TumorArray)

n_PSA = length(params_out.iPatient);
index = params_out.iPatient;
ntumors=length(TumorArray); %thein

% initialize a strutrue for postprocess and adds the first postprocessed
% data as simple tumor volume
simDataPSApost = struct;
for i = 1:n_PSA
    [~,temp,~] = selectbyname(simDataPSA(index(i)).simData, 'V_T');
    simDataPSApost(index(i)).simData.DataNames = {'Voltum'};
    simDataPSApost(index(i)).simData.Data      = temp;
end
% calculates T cell clonality and total CD8 T cell in blood
simDataPSApost = clonality(simDataPSA,simDataPSApost,params_out);

% Add percentage change in tumor size
simDataPSApost = tumSizePerc(simDataPSA,simDataPSApost,params_out,TumorArray);

% Add T cell subset ratio
simDataPSApost = calculate_ratio(simDataPSA,simDataPSApost,params_out);

disp('PSA data post-processed.');
