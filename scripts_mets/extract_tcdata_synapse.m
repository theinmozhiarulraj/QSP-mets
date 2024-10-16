function simDataPSAextract=extract_tcdata_synapse(n_PSA,n_T_specs,TumorArray,LNArray,ntimepoints,index,simDataPSA,simDataPSApost,simDataPSAextract)


    % species of interest, names and compartments
    species_name=["PDL1_total","PDL1_total"];
    varname=["PDL1 in tumor","PDL1 on APCs"];
    compartment_name=["syn_T_C1","syn_T_APC"];
    
    species_name_others=["TPDL1","TPDL1_aPDL1","TPDL1_aPDL1_TPDL1","TPDL1","TPDL1_aPDL1","TPDL1_aPDL1_TPDL1"];
    compartment_name_others=["syn_T_C1","syn_T_C1","syn_T_C1","syn_T_APC","syn_T_APC","syn_T_APC"];
    species_name_derived=["PDL1 on T cells in tumor","PDL1 on T cells in LNs"];
    
    positions=zeros(length(species_name),length(TumorArray));
    positions_others=zeros(length(species_name_others),length(TumorArray));
    pos_tumorvol=zeros(1,length(TumorArray));

    ntimepoints=length(simDataPSA(index(1)).simData.Time);
    
    % get the index using 1 virtual patient
    for array=1:length(TumorArray)
    
        pos_tumorvol(1,array) = find(ismember(simDataPSApost(index(1)).simData.DataNames,['D_T' TumorArray{array}]));
        positions(1,array)=find_index(simDataPSA,species_name(1),convertCharsToStrings(append(compartment_name(1),TumorArray{array})),index(1));
        positions(2,array)=find_index(simDataPSA,species_name(2),convertCharsToStrings(append(compartment_name(2),LNArray{array})),index(1));
    
    end
    
    for j=1:length(species_name_others)
    
        for array=1:length(TumorArray)
            if(j<=3)
                positions_others(j,array)=find_index(simDataPSA,species_name_others(j),convertCharsToStrings(append(compartment_name_others(j),TumorArray{array})),index(1));
            else
                positions_others(j,array)=find_index(simDataPSA,species_name_others(j),convertCharsToStrings(append(compartment_name_others(j),LNArray{array})),index(1));
            end
        end
    
    end
        
    for i = 1:n_PSA
    
        if(i==1)
            len_struct=length(simDataPSAextract);
        end
    
        % number of tumors and LNs
        ntumors_temp=0;
        nLNs_temp=0;
        
        % sum of each quantity from all tumors
        temp_sum=zeros(ntimepoints,length(species_name));
        temp_derived_sum=zeros(ntimepoints,length(species_name_derived));
        
        added=[];
        added_others=[];
    
        for array=1:length(TumorArray)
    
            temp=zeros(ntimepoints,length(species_name));
            temp_others=zeros(ntimepoints,length(species_name_others));
    
            D_T_temp = simDataPSApost(index(i)).simData.Data(:,pos_tumorvol(1,array));
    
            if(D_T_temp(1)>=0.2)
                ntumors_temp=ntumors_temp+1;
    
                for k=1:length(species_name)
    
                    if(k==2)
                        if(~sum(any(strcmp(added,LNArray{array}))))
                            temp(:,k) = simDataPSA(index(i)).simData.Data(:,positions(k,array));
                            temp_sum(:,k)=temp_sum(:,k)+temp(:,k);
                            added=[added;LNArray{array}];
                            nLNs_temp=nLNs_temp+1;
                        end
                    else
                         temp(:,k) = simDataPSA(index(i)).simData.Data(:,positions(k,array));
                         temp_sum(:,k)=temp_sum(:,k)+temp(:,k);
                    end
    
                end
    
               % others: PDL1 expression on T cells
                for q=1:length(species_name_others)
    
                    if(q>3)
                        if(~sum(any(strcmp(added_others,LNArray{array}))))
                            temp_others(:,q) = simDataPSA(index(i)).simData.Data(:,positions_others(q,array));
                            added_others=[added_others;LNArray{array}];
                        end
                    else
                         temp_others(:,q) = simDataPSA(index(i)).simData.Data(:,positions_others(q,array));
                    end
    
                end
    
                % calculate derived quantities
                % PDL1 on T cells in tumor
                PDL1_Tcell_tumor=temp_others(:,1)+temp_others(:,2)+temp_others(:,3);
    
                % PDL1 on T cells in LN
                PDL1_Tcell_LN=temp_others(:,4)+temp_others(:,5)+temp_others(:,6);
    
                temp_derived=[PDL1_Tcell_tumor,PDL1_Tcell_LN];
                temp_derived_sum=temp_derived_sum+temp_derived;
    
            end
    
        end
    
    
        for k_new=1:length(species_name)
    
            index_struct=len_struct+k_new;
    
            if(i==1)
                simDataPSAextract(index_struct).simData.DataNames = {varname(k_new)};
                if(k_new==2)
                    simDataPSAextract(index_struct).simData.Data      = temp_sum(:,k_new)/nLNs_temp;
                else
                    simDataPSAextract(index_struct).simData.Data      = temp_sum(:,k_new)/ntumors_temp;
                end
            else
                if(k_new==2)
                    simDataPSAextract(index_struct).simData.Data      = [simDataPSAextract(index_struct).simData.Data     , temp_sum(:,k_new)/nLNs_temp];
                else
                    simDataPSAextract(index_struct).simData.Data      = [simDataPSAextract(index_struct).simData.Data     , temp_sum(:,k_new)/ntumors_temp];
                end
            end
    
        end
    
        for k_new=1:length(species_name_derived)
    
            index_struct=len_struct+length(species_name)+k_new;
    
            if(i==1)
                simDataPSAextract(index_struct).simData.DataNames = {species_name_derived(k_new)};
                if(k_new==2)
                    simDataPSAextract(index_struct).simData.Data      = temp_derived_sum(:,k_new)/nLNs_temp;                
                else
                    simDataPSAextract(index_struct).simData.Data      = temp_derived_sum(:,k_new)/ntumors_temp;
                end
            else
                if(k_new==2)
                    simDataPSAextract(index_struct).simData.Data      = [simDataPSAextract(index_struct).simData.Data     , temp_derived_sum(:,k_new)/nLNs_temp];
                else
                    simDataPSAextract(index_struct).simData.Data      = [simDataPSAextract(index_struct).simData.Data     , temp_derived_sum(:,k_new)/ntumors_temp];
                end
            end
    
        end
    
    end

end