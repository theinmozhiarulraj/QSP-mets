% Function to calculate and add percentage change in tumor size to
% postprocessed parameters
%
% Inputs: simDataPSA     -- Object containing the whole simbiology model
%                           outputs for all batch simulations
%         simDataPSApost -- Object containing the postprocessed simbiology
%                           model outputs for all batch simulations
%         params_out     -- object containing model outputs to be organized
%                           for future sensitivity analysis
%
% Outputs: simDataPSApost  -- Updated object containing postprocessed
%                             outputs

function simDataPSAout = tumSizePerc(simDataPSA,simDataPSApost,params_out,TumorArray)

n_PSA = length(params_out.iPatient);
index = params_out.iPatient;
simDataPSAout = simDataPSApost;
ntumors=length(TumorArray); 

for i = 1:n_PSA

    % thein: calculate sum of longest diameter of tumors
    sum_longest_dia=[];

    for j=1:ntumors
        % Access tumor volume (microliter = 1e-9*m^3)
        [~,V_T_temp,~] = selectbyname(simDataPSA(index(i)).simData, ['V_T' TumorArray{j}]);
        % Calculate tumor size (diameter) assuming a single sphere nodule (comes out as mm if Tum volume was at microliter)
        D_T_temp = 2*(3/(4*pi)*V_T_temp).^(1/3);
        % Calculates percentage change of tumor size 
        D_T_perc_temp = (D_T_temp - D_T_temp(1))/D_T_temp(1)*100;
        
        % thein: add only if initial tumor size>=0.2 cm
        if(D_T_temp(1)>=0.2)
            [xsize,~]=size(sum_longest_dia);
            if(xsize==0)
                sum_longest_dia=D_T_temp;
            else
                sum_longest_dia=sum_longest_dia+D_T_temp;
            end
        end

        % Add calculated tumor size to postprocess structure
        simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {['D_T' TumorArray{j} '_perc']}];
        simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , D_T_perc_temp];
        simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {['D_T' TumorArray{j}]}];
        simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , D_T_temp];
    end

    D_T_all_perc=(sum_longest_dia - sum_longest_dia(1))/sum_longest_dia(1)*100;
    D_T_all=sum_longest_dia;

    simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {'D_T_all_perc'}];
    simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , D_T_all_perc];

    simDataPSAout(index(i)).simData.DataNames = [simDataPSAout(index(i)).simData.DataNames; {'D_T_all'}];
    simDataPSAout(index(i)).simData.Data      = [simDataPSAout(index(i)).simData.Data     , D_T_all];

end
