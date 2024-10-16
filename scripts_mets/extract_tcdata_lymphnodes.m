function simDataPSAextract=extract_tcdata_lymphnodes(n_PSA,n_T_specs,TumorArray,LNArray,ntimepoints,index,simDataPSA,simDataPSApost,simDataPSAextract)

    volume_lymphnode=1112.65; % mm3

    % species of interest, names
    species_name=["nT0","aT0","T0","IL2","nT1","Th","aTh","aPD1"];
    varname=["Naive CD4+ T cells in LNs","Act. CD4+ T cells in LNs","Tregs in LNs","IL2","Naive CD8+ T cells in LNs","Th in LNs","Act. Th in LNs","Antibody in LNs"];
    
    positions=zeros(length(species_name),length(TumorArray));
    pos_tumorvol=zeros(1,length(TumorArray));
    
    % cytotoxic T cells in tumor
    species_name_Tcyt_LN=[];
    % mAPCs in LN
    species_name_APC_LN=[];
    
    for nt=1:n_T_specs
        species_name_Tcyt_LN=[species_name_Tcyt_LN,append("T",num2str(nt))];
    end
    
    for nt=1:length(TumorArray)
        species_name_APC_LN=[species_name_APC_LN,append("mAPC",TumorArray{nt})];
    end
    
    positions_Tcyt_LN=zeros(length(species_name_Tcyt_LN),length(LNArray));
    positions_APC_LN=zeros(1,length(TumorArray));
    
    
    % get the index using 1 virtual patient
    for j=1:length(species_name_Tcyt_LN)
        for array=1:length(LNArray)
            positions_Tcyt_LN(j,array)=find_index(simDataPSA,species_name_Tcyt_LN(j),convertCharsToStrings(['V_LN' LNArray{array}]),index(1));
        end
    end
    
    for array=1:length(TumorArray)
        positions_APC_LN(1,array)=find_index(simDataPSA,species_name_APC_LN(array),convertCharsToStrings(['V_LN' LNArray{array}]),index(1));
    end
    
    species_name_derived=["Tcyt in LNs","Treg/Tcyt in LNs","Th/Tcyt in LNs","Treg/Th in LNs",...
       "Th+Tcyt in LNs","T cells in LNs","(Th+Tcyt)/T cells in LNs","Th/T cells in LNs","Tcyt/T cells in LNs","Treg/T cells in LNs","APCs in LNs",...
       "Tcyt richness in LNs","Tcyt diversity in LNs","Tcyt evenness in LNs"];
    
    ntimepoints=length(simDataPSA(index(1)).simData.Time);
    
    % get the index using 1 virtual patient
    for array=1:length(TumorArray)
    
        pos_tumorvol(1,array) = find(ismember(simDataPSApost(index(1)).simData.DataNames,['D_T' TumorArray{array}]));
        for j=1:length(species_name)
            positions(j,array)=find_index(simDataPSA,species_name(j),convertCharsToStrings(['V_LN' LNArray{array}]),index(1));
        end
    
    end
        
    for i = 1:n_PSA
    
        if(i==1)
            len_struct=length(simDataPSAextract);
        end
    
        nLNs_temp=0;
        added=[];
        temp_sum=zeros(ntimepoints,length(species_name));
        temp_derived_sum=zeros(ntimepoints,length(species_name_derived));
    
        for array=1:length(TumorArray)
    
            temp=zeros(ntimepoints,length(species_name));
            temp_tcyt=zeros(ntimepoints,n_T_specs);
    
            D_T_temp = simDataPSApost(index(i)).simData.Data(:,pos_tumorvol(1,array));
    
            temp_apc=zeros(ntimepoints,1);

            if(D_T_temp(1)>=0.2) 
    
                % mAPCs in LNs
                temp_apc = simDataPSA(index(i)).simData.Data(:,positions_APC_LN(1,array));
    
                if(~sum(any(strcmp(added,LNArray{array}))))
    
                    added=[added;LNArray{array}];
                    nLNs_temp=nLNs_temp+1;
    
                    for k=1:length(species_name)
                        temp(:,k) = simDataPSA(index(i)).simData.Data(:,positions(k,array));

                        if(k==1 || k==2 || k==3 || k==5 || k==6 || k==7)
                            
                            temp_sum(:,k) = temp_sum(:,k)+(temp(:,k)/volume_lymphnode);

                        else

                            temp_sum(:,k) = temp_sum(:,k)+temp(:,k);

                        end

                    end
        
                    for k=1:n_T_specs
                        temp_tcyt(:,k) = simDataPSA(index(i)).simData.Data(:,positions_Tcyt_LN(k,array));
                    end
    
                    % cytotoxic t cells 
                    %Tcyt_LN=sum(temp_tcyt,2);
    
                    % richness, diversity and evenness of Cytotoxic T cells
                    richness_tcyt_LN=zeros(ntimepoints,1);
                    Tcyt_LN=zeros(ntimepoints,1);
                    shannon_index_tcyt_LN=zeros(ntimepoints,1);
                
                    for g=1:n_T_specs
                        idx=find(temp_tcyt(:,g)>=1);
                        if(length(idx)>=1)
                            richness_tcyt_LN(idx,1)=richness_tcyt_LN(idx,1)+1;
                            Tcyt_LN(idx,1)=Tcyt_LN(idx,1)+temp_tcyt(idx,g);
                        end
                    end
                
                    for g=1:n_T_specs
                        idx=find(temp_tcyt(:,g)>=1);
                        if(length(idx)>=1)
                            proportion=temp_tcyt(idx,g)./Tcyt_LN(idx,1);
                            shannon_index_tcyt_LN(idx,1)=shannon_index_tcyt_LN(idx,1)-(proportion.*log(proportion));
                        end
                    end
    
                    % evenness of T cyt in central comp.
                    evenness_tcyt_LN=shannon_index_tcyt_LN./log(richness_tcyt_LN);
    
                    % derived quantities
                    % Treg/Tcyt ratio
                    Treg_Tcyt_ratio_LN=temp(:,3)./Tcyt_LN;
                    
                    % Th/Tcyt ratio
                    Th_Tcyt_ratio_LN=temp(:,6)./Tcyt_LN;
                    
                    % Treg/Th ratio
                    Treg_Th_ratio_LN=temp(:,3)./temp(:,6);
                    
                    % Th+Tcyt 
                    Th_Tcyt_sum_LN=temp(:,6)+Tcyt_LN;
                    
                    % Th+Tcyt+Treg
                    Th_Treg_Tcyt_sum_LN=temp(:,6)+temp(:,3)+Tcyt_LN;
                    
                    % (Th+Tcyt)/T total
                    Tact_Ttotal_ratio=(temp(:,6)+Tcyt_LN)./Th_Treg_Tcyt_sum_LN;
                    
                    % Th/Ttotal 
                    Th_Ttotal_ratio=temp(:,6)./Th_Treg_Tcyt_sum_LN;
                    
                    % Tcyt/Ttotal
                    Tcyt_Ttotal_ratio=Tcyt_LN./Th_Treg_Tcyt_sum_LN;
                    
                    % Treg/Ttotal
                    Treg_Ttotal_ratio=temp(:,3)./Th_Treg_Tcyt_sum_LN;
    
                   temp_derived=[Tcyt_LN/volume_lymphnode,Treg_Tcyt_ratio_LN,Th_Tcyt_ratio_LN,Treg_Th_ratio_LN,...
                   Th_Tcyt_sum_LN/volume_lymphnode,Th_Treg_Tcyt_sum_LN/volume_lymphnode,Tact_Ttotal_ratio,Th_Ttotal_ratio,Tcyt_Ttotal_ratio,Treg_Ttotal_ratio,temp_apc/volume_lymphnode,...
                   richness_tcyt_LN,shannon_index_tcyt_LN,evenness_tcyt_LN];
    
                   temp_derived_sum=temp_derived_sum+temp_derived;
    
                end
                
            end
            
        end
    
    
         for k_new=1:length(species_name)
    
            index_struct=len_struct+k_new;
    
            if(i==1)
                simDataPSAextract(index_struct).simData.DataNames = {varname(k_new)};
                simDataPSAextract(index_struct).simData.Data      = temp_sum(:,k_new)/nLNs_temp;
            else
                simDataPSAextract(index_struct).simData.Data      = [simDataPSAextract(index_struct).simData.Data     , temp_sum(:,k_new)/nLNs_temp];
            end
    
         end
    
         % derived quantities
    
         for k_new=1:length(species_name_derived)
    
            index_struct=len_struct+k_new+length(species_name);
    
            if(i==1)
                simDataPSAextract(index_struct).simData.DataNames = {species_name_derived(k_new)};
                simDataPSAextract(index_struct).simData.Data      = temp_derived_sum(:,k_new)/nLNs_temp;
            else
                simDataPSAextract(index_struct).simData.Data      = [simDataPSAextract(index_struct).simData.Data     , temp_derived_sum(:,k_new)/nLNs_temp];
            end
    
         end
    
    end

end
