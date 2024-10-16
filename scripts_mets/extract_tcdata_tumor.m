function simDataPSAextract=extract_tcdata_tumor(n_PSA,params_in,n_T_specs,cancer_clones,TumorArray,ntimepoints,index,simDataPSA,simDataPSApost,simDataPSAextract)

    % species of interest, names
    species_name=["C_x","T1_exh","Th_exh","K",...
        "c_vas","T0","APC","mAPC","c","P0",...
        "aPD1","Th","TGFb","IFNg","MDSC","NO","ArgI","CCL2","Mac_M1","Mac_M2","IL12","IL10"];

    varname=["Dead cancer cell","Exhausted Tcyt","Exhausted Th","Carrying capacity",...
        "Angiogenic factor","Treg","immature APC","mAPC","APC maturation cytokine","Self Ag",...
        "anti-PD-1 antibody","Th","TGFb","IFNg","MDSC","NO","ArgI","CCL2","M1 Mac","M2 Mac","IL12","IL10"];
    
    ntimepoints=length(simDataPSA(index(1)).simData.Time);
    
    % get the index using 1 virtual patient
    positions=zeros(length(species_name),length(TumorArray));
    pos_tumorvol=zeros(1,length(TumorArray)); % tumor dia
    pos_tumorVT=zeros(1,length(TumorArray)); % tumor volume in ml?
   
    for array=1:length(TumorArray)
    
        pos_tumorvol(1,array) = find(ismember(simDataPSApost(index(1)).simData.DataNames,['D_T' TumorArray{array}]));
        pos_tumorVT(1,array) = find(ismember(simDataPSA(index(1)).simData.DataNames,['V_T' TumorArray{array}]));

        for j=1:length(species_name)
            positions(j,array)=find_index(simDataPSA,species_name(j),convertCharsToStrings(['V_T' TumorArray{array}]),index(1));
        end
    
    end
    
    % cytotoxic T cells in tumor
    species_name_Tcyt_tumor=[];
    % antigen in tumor
    species_name_Ag_tumor=[];
    % cancer cells in tumor
    species_name_C_tumor=[];
    
    for nt=1:n_T_specs
        species_name_Tcyt_tumor=[species_name_Tcyt_tumor,append("T",num2str(nt))];
        species_name_Ag_tumor=[species_name_Ag_tumor,append("P",num2str(nt))];
    end
    
    for nc=1:length(cancer_clones)
        species_name_C_tumor=[species_name_C_tumor,append("C",num2str(nc))];
    end
    
    positions_Tcyt_tumor=zeros(length(species_name_Tcyt_tumor),length(TumorArray));
    positions_Ag_tumor=zeros(length(species_name_Ag_tumor),length(TumorArray));
    positions_C_tumor=zeros(length(species_name_C_tumor),length(TumorArray));
    
    
    % get the index using 1 virtual patient
    for j=1:length(species_name_Tcyt_tumor)
        for array=1:length(TumorArray)
            positions_Tcyt_tumor(j,array)=find_index(simDataPSA,species_name_Tcyt_tumor(j),convertCharsToStrings(['V_T' TumorArray{array}]),index(1));
            positions_Ag_tumor(j,array)=find_index(simDataPSA,species_name_Ag_tumor(j),convertCharsToStrings(['V_T' TumorArray{array}]),index(1));
        end
    end
    
    for j=1:length(species_name_C_tumor)
        for array=1:length(TumorArray)
            positions_C_tumor(j,array)=find_index(simDataPSA,species_name_C_tumor(j),convertCharsToStrings(['V_T' TumorArray{array}]),index(1));
        end
    end
    
    species_name_derived=["Tcyt","Ag","Cancer cells","Treg/Tcyt","Th/Tcyt","Treg/Th",...
       "Th+Tcyt","T cells","(Th+Tcyt)/T cells","Th/T cells","Tcyt/T cells","Treg/T cells",...
       "Mac","M2/M1 Mac ratio","M1/(M1+M2) Mac","M2/(M1+M2) Mac","Total number of cells in tumor","TIL fraction","Immune supp. cells","Immune act. cells","Immune supp./act. cells","Immune supp. cell fraction","Immune act. cell fraction","MDSC fraction","M1 Mac fraction","M2 Mac fraction",...
       "Mac fraction","Treg fraction","Tcyt fraction","Th fraction","Tcyt+Th cell fraction","mAPC fraction","APC fraction","Exhausted Tcyt fraction","Exhausted Th fraction","Immune cell fraction",...
       "Cancer clone richness","Tcyt richness","Cancer clone diversity","Tcyt diversity","Cancer clones evenness","Tcyt evenness",...
       "neo-Ag diversity","neo-Ag evenness","Tumor diameter","Neo-Ag specific T cell clones"];
    
    
    for i = 1:n_PSA
    
        if(i==1)
            len_struct=length(simDataPSAextract);
        end
    
        % number of tumors
        ntumors_temp=0;
    
        % sum of each quantity from all tumors
        temp_sum=zeros(ntimepoints,length(species_name));
        temp_derived_sum=zeros(ntimepoints,length(species_name_derived));
    
        for array=1:length(TumorArray)
    
            temp=NaN(ntimepoints,length(species_name));
            temp_tcyt=NaN(ntimepoints,n_T_specs);
            temp_ag=NaN(ntimepoints,n_T_specs);
            temp_cancer=NaN(ntimepoints,length(cancer_clones));
    
            % tumor diameter
            D_T_temp = simDataPSApost(index(i)).simData.Data(:,pos_tumorvol(1,array));
            vol_temp = simDataPSA(index(i)).simData.Data(:,pos_tumorVT(1,array));

            if(D_T_temp(1)>=0.2)
    
                ntumors_temp=ntumors_temp+1;
    
                for k=1:length(species_name)
                    temp(:,k) = simDataPSA(index(i)).simData.Data(:,positions(k,array));
                    %temp_sum(:,k) = temp_sum(:,k)+temp(:,k);
                end
    
                % calculate derived quantities
    
                for k=1:n_T_specs
                    temp_tcyt(:,k) = simDataPSA(index(i)).simData.Data(:,positions_Tcyt_tumor(k,array));
                    temp_ag(:,k) = simDataPSA(index(i)).simData.Data(:,positions_Ag_tumor(k,array));
                end

                % set -ve values in temp_ag to 0
                temp_ag(temp_ag<0)=0;
    
                for k=1:length(cancer_clones)
                    temp_cancer(:,k) = simDataPSA(index(i)).simData.Data(:,positions_C_tumor(k,array));
                end
    
    
                sum_TMB=zeros(ntimepoints,1);
                for k=1:n_T_specs
    
                    % calculate sum TMB
                    sum_TMB(:,1)=sum_TMB(:,1)+repelem(params_in.(['n_T' num2str(k) '_clones']).LHS(index(i)),ntimepoints)';
    
                end
    
    
                % cytotoxic t cells & Ag concentrations 
                %Tcyt_tumor=sum(temp_tcyt,2);
                Ag_tumor=sum(temp_ag,2);
                % cancercells_tumor=sum(temp_cancer,2);
    
                % richness, diversity and evenness of cancer cells
                richness_cancer=zeros(ntimepoints,1);
                cancercells_tumor=zeros(ntimepoints,1);
                shannon_index_cancer=zeros(ntimepoints,1);
    
                for g=1:length(cancer_clones)
                    idx=find(temp_cancer(:,g)>=1);
                    if(length(idx)>=1)
                        richness_cancer(idx,1)=richness_cancer(idx,1)+1;
                        cancercells_tumor(idx,1)=cancercells_tumor(idx,1)+temp_cancer(idx,g);
                    end
                end
    
                for g=1:length(cancer_clones)
                    idx=find(temp_cancer(:,g)>=1);
                    if(length(idx)>=1)
                        proportion=temp_cancer(idx,g)./cancercells_tumor(idx,1);
                        shannon_index_cancer(idx,1)=shannon_index_cancer(idx,1)-(proportion.*log(proportion));
                    end
                end
    
                % evenness of cancer cells
                evenness_cancer=shannon_index_cancer./log(richness_cancer);
    
                % richness, diversity and evenness of Cytotoxic T cells
                richness_tcyt=zeros(ntimepoints,1);
                Tcyt_tumor=zeros(ntimepoints,1);
                shannon_index_tcyt=zeros(ntimepoints,1);
    
                for g=1:n_T_specs
                    idx=find(temp_tcyt(:,g)>=1);
                    if(length(idx)>=1)
                        richness_tcyt(idx,1)=richness_tcyt(idx,1)+1;
                        Tcyt_tumor(idx,1)=Tcyt_tumor(idx,1)+temp_tcyt(idx,g);
                    end
                end
    
                for g=1:n_T_specs
                    idx=find(temp_tcyt(:,g)>=1);
                    if(length(idx)>=1)
                        proportion=temp_tcyt(idx,g)./Tcyt_tumor(idx,1);
                        shannon_index_tcyt(idx,1)=shannon_index_tcyt(idx,1)-(proportion.*log(proportion));
                    end
                end
    
                % evenness of Tcyt cells
                evenness_tcyt=shannon_index_tcyt./log(richness_tcyt);
    
                % diversity of neoag in tumor
                shannon_index_ag=zeros(ntimepoints,1);
                for g=1:n_T_specs
                    proportion=temp_ag(:,g)./Ag_tumor;
                    shannon_index_ag(:,1)=shannon_index_ag(:,1)-(proportion.*log(proportion));
                end
    
                % evenness of neoag
                evenness_ag=shannon_index_ag./log(n_T_specs);
    
                % Treg/Tcyt ratio
                Treg_Tcyt_ratio_tumor=temp(:,6)./Tcyt_tumor;
                
                % Th/Tcyt ratio
                Th_Tcyt_ratio_tumor=temp(:,12)./Tcyt_tumor;
                
                % Treg/Th ratio
                Treg_Th_ratio_tumor=temp(:,6)./temp(:,12);
                
                % Th+Tcyt 
                Th_Tcyt_sum_tumor=temp(:,12)+Tcyt_tumor;
                
                % Th+Tcyt+Treg+exhausted T cells
                Th_Treg_Tcyt_sum_tumor=temp(:,12)+temp(:,6)+Tcyt_tumor+temp(:,2)+temp(:,3);
                
                % (Th+Tcyt)/T total
                Tact_Ttotal_ratio=(temp(:,12)+Tcyt_tumor)./Th_Treg_Tcyt_sum_tumor;
                
                % Th/Ttotal 
                Th_Ttotal_ratio=temp(:,12)./Th_Treg_Tcyt_sum_tumor;
                
                % Tcyt/Ttotal
                Tcyt_Ttotal_ratio=Tcyt_tumor./Th_Treg_Tcyt_sum_tumor;
                
                % Treg/Ttotal
                Treg_Ttotal_ratio=temp(:,6)./Th_Treg_Tcyt_sum_tumor;
    
                % sum of M1 and M2 macrophages
                M1_M2_sum=temp(:,19)+temp(:,20);
    
                % M2/M1 ratio
                M2_M1_ratio=temp(:,20)./temp(:,19);
    
                % Fraction of M1 among all macrophages
                Frac_M1_Mac=temp(:,19)./M1_M2_sum;
    
                % Fraction of M2 among all macrophages
                Frac_M2_Mac=temp(:,20)./M1_M2_sum;
    
                % count all cells in tumor
                total_cells=Tcyt_tumor+cancercells_tumor+temp(:,1)+temp(:,2)+temp(:,3)+temp(:,6)+temp(:,7)+temp(:,8)+...
                    temp(:,12)+temp(:,15)+temp(:,19)+temp(:,20);
    
                % fraction of tumor infiltrating T cells
                TILs=Th_Treg_Tcyt_sum_tumor./total_cells;
    
                % total suppressive cell types MDSCs+M2 macs+Tregs
                supp_cells=temp(:,15)+temp(:,20)+temp(:,6);
    
                % Immune response promoting cell types (APCs+mAPCs+Tcyt+Th+Mac_M1+T1_exh+Th_exh)
                act_cells=temp(:,7)+temp(:,8)+Tcyt_tumor+temp(:,12)+temp(:,19)+temp(:,2)+temp(:,3);
    
                % suppressive/immune response promoting cell types
                supp_act_ratio=supp_cells./act_cells;
    
                % fraction of suppressive cells
                supp_frac=supp_cells./total_cells;
    
                % fraction of activating cells
                act_frac=act_cells./total_cells;
    
                % fraction of MDSCs
                MDSC_frac=temp(:,15)./total_cells;
    
                % fraction of M1 macrophages
                M1_frac=temp(:,19)./total_cells;
    
                % fraction of M2 macrophages
                M2_frac=temp(:,20)./total_cells;
    
                % fraction of macrophages M1+M2
                M1M2_frac=M1_M2_sum./total_cells;
    
                % fraction of Tregs
                Tregs_frac=temp(:,6)./total_cells;
    
                % fraction of Tcyt
                Tcyt_frac=Tcyt_tumor./total_cells;
    
                % fraction of Th
                Th_frac=temp(:,12)./total_cells;
    
                % fraction of (Tcyt+Th)
                Tcyt_Th_frac=Th_Tcyt_sum_tumor./total_cells;
    
                % fraction of mAPCs
                mAPC_frac=temp(:,8)./total_cells;
    
                % fraction of total APCs
                APC_mAPC_frac=(temp(:,8)+temp(:,7))./total_cells;
    
                % fraction of exhausted CD8+T cells
                T1exh_frac=temp(:,2)./total_cells;
    
                % fraction of exhausted Th cells
                Thexh_frac=temp(:,3)./total_cells;
    
                %fraction of immune cells
                immune_frac=1-(cancercells_tumor./total_cells);
    
                temp_derived=[Tcyt_tumor./vol_temp,Ag_tumor,cancercells_tumor./vol_temp,Treg_Tcyt_ratio_tumor,Th_Tcyt_ratio_tumor,Treg_Th_ratio_tumor,...
                   Th_Tcyt_sum_tumor./vol_temp,Th_Treg_Tcyt_sum_tumor./vol_temp,Tact_Ttotal_ratio,Th_Ttotal_ratio,Tcyt_Ttotal_ratio,Treg_Ttotal_ratio,...
                   M1_M2_sum./vol_temp,M2_M1_ratio,Frac_M1_Mac,Frac_M2_Mac,total_cells,TILs,supp_cells./vol_temp,act_cells./vol_temp,supp_act_ratio,supp_frac,act_frac,MDSC_frac,M1_frac,M2_frac,...
                   M1M2_frac,Tregs_frac,Tcyt_frac,Th_frac,Tcyt_Th_frac,mAPC_frac,APC_mAPC_frac,T1exh_frac,Thexh_frac,immune_frac,...
                   richness_cancer,richness_tcyt,shannon_index_cancer,shannon_index_tcyt,evenness_cancer,evenness_tcyt,...
                   shannon_index_ag,evenness_ag,D_T_temp,sum_TMB];
    
                temp_derived_sum=temp_derived_sum+temp_derived;

               for k=1:length(species_name)

                   if(k==1 || k==2 || k==3 || k==6 || k==7 || k==8 || k==12 || k==15 || k==19 || k==20)

                        temp_sum(:,k) = temp_sum(:,k)+(temp(:,k)./vol_temp);

                   else

                        temp_sum(:,k) = temp_sum(:,k)+temp(:,k);

                   end

               end
    
            end
    
        end
    
        for k_new=1:length(species_name)
    
            index_struct=len_struct+k_new;
    
            if(i==1)
                simDataPSAextract(index_struct).simData.DataNames = {varname(k_new)};
                simDataPSAextract(index_struct).simData.Data      = temp_sum(:,k_new)/ntumors_temp;
            else
                simDataPSAextract(index_struct).simData.Data      = [simDataPSAextract(index_struct).simData.Data     , temp_sum(:,k_new)/ntumors_temp];
            end
    
        end
    
        for k_new=1:length(species_name_derived)
    
            index_struct=len_struct+length(species_name)+k_new;
    
            if(i==1)
                simDataPSAextract(index_struct).simData.DataNames = {species_name_derived(k_new)};
                simDataPSAextract(index_struct).simData.Data      = temp_derived_sum(:,k_new)/ntumors_temp;
            else
                simDataPSAextract(index_struct).simData.Data      = [simDataPSAextract(index_struct).simData.Data     , temp_derived_sum(:,k_new)/ntumors_temp];
            end
    
        end
    
    end

end
