function simDataPSAextract=extract_tcdata_central_periphery(n_PSA,n_T_specs,ntimepoints,index,simDataPSA,simDataPSAextract)

    volume_central=5*1000; % 5000 ml = 5 l

    % species of interest, names and compartments
    species_name=["nT0","nT1","T0","Th","aPD1","nT0","nT1","T0","Th","aPD1"];

    varname=["Naive CD4+ T cells in central comp.","Naive CD8+ T cells in central comp.",...
        "Tregs in central comp.","Th in central comp.","Antibody in central comp.",...
        "Naive CD4+ T cells in periphery","Naive CD8+ T cells in periphery","Tregs in periphery",...
        "Th in periphery","Antibody in periphery"];

    compartment_name=["V_C","V_C","V_C","V_C","V_C","V_P","V_P","V_P","V_P","V_P"];
    
    % get the index using 1 virtual patient
    positions=zeros(length(species_name),1);
    for j=1:length(species_name)
        positions(j,1)=find_index(simDataPSA,species_name(j),compartment_name(j),index(1));
    end
    
    % cytotoxic T cells
    species_name_Tcyt_central=[];
    for nt=1:n_T_specs
        species_name_Tcyt_central=[species_name_Tcyt_central,append("T",num2str(nt))];
    end
    
    % get the index using 1 virtual patient
    positions_Tcyt_central=zeros(length(species_name_Tcyt_central),1);
    for j=1:length(species_name_Tcyt_central)
        positions_Tcyt_central(j,1)=find_index(simDataPSA,species_name_Tcyt_central(j),"V_C",index(1));
    end
    
        
    for i = 1:n_PSA
    
        temp=NaN(ntimepoints,length(species_name));
    
        % store timecourse of species of interest
        for k=1:length(species_name)
                
               temp(:,k) = simDataPSA(index(i)).simData.Data(:,positions(k,1));

               if(k<5)

                   if(i==1)
                       simDataPSAextract(k).simData.DataNames = {varname(k)};
                       simDataPSAextract(k).simData.Data      = temp(:,k)/volume_central;
                   else
                       simDataPSAextract(k).simData.Data      = [simDataPSAextract(k).simData.Data     , temp(:,k)/volume_central];
                   end

               else

                   if(i==1)
                       simDataPSAextract(k).simData.DataNames = {varname(k)};
                       simDataPSAextract(k).simData.Data      = temp(:,k);
                   else
                       simDataPSAextract(k).simData.Data      = [simDataPSAextract(k).simData.Data     , temp(:,k)];
                   end

               end

        
        end
    
        temp_Tcyt_central=NaN(ntimepoints,length(species_name_Tcyt_central));
    
        % extract cytotoxic T cell counts (different antigen specificities)
        for p=1:length(species_name_Tcyt_central)
               temp_Tcyt_central(:,p) = simDataPSA(index(i)).simData.Data(:,positions_Tcyt_central(p,1));
        end
    
        % total cytotoxic T cell count in central comp. (includes T cells with different Ag specificities)
        % total_Tcyt_central=sum(temp_Tcyt_central,2); % sum of rows
    
        % richness, diversity and evenness of Cytotoxic T cells
        richness_tcyt_central=zeros(ntimepoints,1);
        total_Tcyt_central=zeros(ntimepoints,1);
        shannon_index_tcyt_central=zeros(ntimepoints,1);
    
        for g=1:n_T_specs
            idx=find(temp_Tcyt_central(:,g)>=1);
            if(length(idx)>=1)
                richness_tcyt_central(idx,1)=richness_tcyt_central(idx,1)+1;
                total_Tcyt_central(idx,1)=total_Tcyt_central(idx,1)+temp_Tcyt_central(idx,g);
            end
        end
    
        for g=1:n_T_specs
            idx=find(temp_Tcyt_central(:,g)>=1);
            if(length(idx)>=1)
                proportion=temp_Tcyt_central(idx,g)./total_Tcyt_central(idx,1);
                shannon_index_tcyt_central(idx,1)=shannon_index_tcyt_central(idx,1)-(proportion.*log(proportion));
            end
        end
    
        evenness_tcyt_central=shannon_index_tcyt_central./log(richness_tcyt_central);
    
        % calculate ratios
        % Treg/Tcyt ratio
        Treg_Tcyt_ratio_central=temp(:,3)./total_Tcyt_central;
    
        % Th/Tcyt ratio
        Th_Tcyt_ratio_central=temp(:,4)./total_Tcyt_central;
    
        % Treg/Th ratio
        Treg_Th_ratio_central=temp(:,3)./temp(:,4);
    
        % Th+Tcyt 
        Th_Tcyt_sum_central=temp(:,4)+total_Tcyt_central;
    
        % Th+Tcyt+Treg 
        Th_Treg_Tcyt_sum_central=temp(:,3)+temp(:,4)+total_Tcyt_central;
        
        % (Th+Tcyt)/Th+Tcyt+Treg
        Tact_Ttotal_ratio=(temp(:,4)+total_Tcyt_central)./Th_Treg_Tcyt_sum_central;
    
        % Th/Ttotal 
        Th_Ttotal_ratio=temp(:,4)./Th_Treg_Tcyt_sum_central;
    
        % Tcyt/Ttotal
        Tcyt_Ttotal_ratio=total_Tcyt_central./Th_Treg_Tcyt_sum_central;
    
        % Treg/Ttotal
        Treg_Ttotal_ratio=temp(:,3)./Th_Treg_Tcyt_sum_central;
    
        ratios_temp=[total_Tcyt_central/volume_central,Treg_Tcyt_ratio_central,Th_Tcyt_ratio_central,Treg_Th_ratio_central,Th_Tcyt_sum_central/volume_central,Th_Treg_Tcyt_sum_central/volume_central,...
            Tact_Ttotal_ratio,Th_Ttotal_ratio,Tcyt_Ttotal_ratio,Treg_Ttotal_ratio,richness_tcyt_central,shannon_index_tcyt_central,evenness_tcyt_central];
    
        ratios_names=["Tcyt in central comp.","Treg/Tcyt in central comp.","Th/Tcyt in central comp.","Treg/Th in central comp.","Th+Tcyt in central comp.","T cells in central comp.",...
            "(Th+Tcyt)/T cells in central comp.","Th/T cells in central comp.","Tcyt/T cells in central comp.","Treg/T cells in central comp.","Tcyt richness in central comp.","Tcyt diversity in central comp.","Tcyt evenness in central comp."];
    
        [~,len_ratios_temp]=size(ratios_temp);
        % Length should be the same
        if(len_ratios_temp~=length(ratios_names))
        
            disp("ERROR: check the length of variables and names")
        
        end
    
        % add all the ratios and sum of cell counts to the struct
        if(i==1)
            len_struct=length(simDataPSAextract);
        end
    
        for q=1:length(ratios_names)
    
            index_struct=len_struct+q;
        
            if(i==1)
               simDataPSAextract(index_struct).simData.DataNames = {ratios_names(q)};
               simDataPSAextract(index_struct).simData.Data      = ratios_temp(:,q);
            
            else
               simDataPSAextract(index_struct).simData.Data      = [simDataPSAextract(index_struct).simData.Data     , ratios_temp(:,q)];
            end
        
        end
    
    end

end