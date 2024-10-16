% overall readouts of all tumors together

% total number of virtual patients generated
% includes the ones that did not reach initial conditions

total_patients=length(simDataPSA);

% id for each patient
patient_id=transpose(1:total_patients);

% recist status
recist=params_out.RECIST_mets;

% shannon index of neoag
shannon_index_neoag=zeros(total_patients,1);
total_neoag =zeros(total_patients,n_T_specs);
proportions_neoag=zeros(total_patients,n_T_specs);

% get the positions using 1 virtual patient
positions_cancer_clones=zeros(length(cancer_clones),length(TumorArray));
positions_neoag=zeros(length(cancer_clones),n_T_specs);

done=false;
i=1;

while(done==false)

    if(params_out.patient(i)==1)

        done=true;

        for array=1:length(TumorArray)

            for j=1:length(cancer_clones)
                positions_cancer_clones(j,array)=find_index(simDataPSA,['C' num2str(j)],convertCharsToStrings(['V_T' TumorArray{array}]),i);
            end

        end
    
    end

    i=i+1;

end

for i = 1:total_patients

    % if the patient reached initial tumor size
    if(params_out.patient(i)==1)

        for array = 1:length(TumorArray)

            pos = find(ismember(simDataPSApost(i).simData.DataNames,['D_T' TumorArray{j}]));
            D_T_temp = simDataPSApost(i).simData.Data(:,pos);

            if(D_T_temp(1)>=0.2)

                for k=1:n_T_specs
    
                    for j = 1:length(cancer_clones)

                        neoag_conc_perclone=params_in.(['agconc_C' num2str(j) '_' num2str(k)]).LHS(i);
                        
                        %pos = find_index(simDataPSA,['C' num2str(j)],convertCharsToStrings(['V_T' TumorArray{array}]),i);
                        total_neoag(i,k) = total_neoag(i,k)+(simDataPSA(i).simData.Data(1,positions_cancer_clones(j,array))*neoag_conc_perclone);
                        
                    end

                end

            end

        end
 
    end

end

for i=1:total_patients

    for k=1:n_T_specs
    
        proportions_neoag(i,k)=total_neoag(i,k)/sum(total_neoag(i,:));
        shannon_index_neoag(i,1)=shannon_index_neoag(i,1)-(proportions_neoag(i,k)*log(proportions_neoag(i,k)));
    
    end

end

% calculate evenness 
evenness_neoag=shannon_index_neoag/log(n_T_specs);  

dataset_diversity_mets=table(patient_id,recist,shannon_index_neoag,evenness_neoag);

writetable(dataset_diversity_mets,'test_diversity_mets_neoag.csv')