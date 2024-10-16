% this script counts the number of patients with lung metastases

% total number of patients
total_patients=length(params_out.iPatient);

lung_mets_suffix=["_Ln1","_Ln2"];

npatients_with_lungmets=0;

for i = 1:total_patients

    count=0;


    for array = 1:length(lung_mets_suffix)

        pos = find(ismember(simDataPSApost(params_out.iPatient(i)).simData.DataNames,['D_T' lung_mets_suffix{array}]));
        D_T_temp = simDataPSApost(params_out.iPatient(i)).simData.Data(:,pos);

        % if the tumor dia is > 0.2 cm, then the tumor is present
        if(D_T_temp(1)>=0.2)
            
            count=count+1;

        end

    end


    if(count>0)

        npatients_with_lungmets=npatients_with_lungmets+1;

    end

end

frac_lungmets=npatients_with_lungmets/total_patients;

disp(frac_lungmets)