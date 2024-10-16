organs=["lung1","lung2","others"];
suffix=["_Ln1","_Ln2","_other"];

for j=1:length(organs)

    npatients = length(params_out.iPatient);
    
    % find the number of timepoints in data
    ntimepoints=length(simDataPSA(params_out.iPatient(1)).simData.Data(:,1));
    
    t_time=transpose(time);
    
    MDSC_ratio=t_time;
    Teff_ratio=t_time;
    Treg_ratio=t_time;
    Th_ratio=t_time;
    Mac_ratio=t_time;
    MacM1_ratio=t_time;
    MacM2_ratio=t_time;
    mAPC_ratio=t_time;
    APC_ratio=t_time;
    allAPC_ratio=t_time;
    total_tcells_ratio=t_time;
    total_teff_ratio=t_time;
    total_th_ratio=t_time;
    total_cd4_ratio=t_time;

    for i = 1:npatients
    
        index = params_out.iPatient(i);
    
        if(i==1)
    
            pos_MDSC = find_index(simDataPSA,"MDSC","V_T",index);
            pos_Teff=find(ismember(simDataPSA(index).simData.DataNames,'Tcyt_total'));
            pos_Th=find_index(simDataPSA,"Th","V_T",index);
            pos_Teff_exh=find_index(simDataPSA,"T1_exh","V_T",index);
            pos_Th_exh=find_index(simDataPSA,"Th_exh","V_T",index);
            pos_Ctotal=find(ismember(simDataPSA(index).simData.DataNames,'C_total'));
            pos_Cx=find_index(simDataPSA,"C_x","V_T",index);
            pos_Treg=find_index(simDataPSA,"T0","V_T",index);
            pos_MacM1=find_index(simDataPSA,"Mac_M1","V_T",index);
            pos_MacM2=find_index(simDataPSA,"Mac_M2","V_T",index);
            pos_mAPC=find_index(simDataPSA,"mAPC","V_T",index);
            pos_APC=find_index(simDataPSA,"APC","V_T",index);
    
            pos_MDSC_met = find_index(simDataPSA,"MDSC",append('V_T',suffix(j)),index);
            pos_Teff_met=find(ismember(simDataPSA(index).simData.DataNames,append('Tcyt_total',suffix(j))));
            pos_Th_met=find_index(simDataPSA,"Th",append('V_T',suffix(j)),index);
            pos_Teff_exh_met=find_index(simDataPSA,"T1_exh",append('V_T',suffix(j)),index);
            pos_Th_exh_met=find_index(simDataPSA,"Th_exh",append('V_T',suffix(j)),index);
            pos_Ctotal_met=find(ismember(simDataPSA(index).simData.DataNames,append('C_total',suffix(j))));
            pos_Cx_met=find_index(simDataPSA,"C_x",append('V_T',suffix(j)),index);
            pos_Treg_met=find_index(simDataPSA,"T0",append('V_T',suffix(j)),index);
            pos_MacM1_met=find_index(simDataPSA,"Mac_M1",append('V_T',suffix(j)),index);
            pos_MacM2_met=find_index(simDataPSA,"Mac_M2",append('V_T',suffix(j)),index);
            pos_mAPC_met=find_index(simDataPSA,"mAPC",append('V_T',suffix(j)),index);
            pos_APC_met=find_index(simDataPSA,"APC",append('V_T',suffix(j)),index);
    
        end
    
        %find(ismember(simDataPSApost(i).simData.DataNames,['D_T' TumorArray{j}]));
    
        pos_dia = find(ismember(simDataPSApost(index).simData.DataNames,append('D_T',suffix(j))));
        D_T_temp = simDataPSApost(index).simData.Data(:,pos_dia);

        if(D_T_temp(1)>=0.2)

            MDSC_primary=simDataPSA(index).simData.Data(:,pos_MDSC);
            Teff_primary=simDataPSA(index).simData.Data(:,pos_Teff);
            Th_primary=simDataPSA(index).simData.Data(:,pos_Th);
            Teff_exh_primary=simDataPSA(index).simData.Data(:,pos_Teff_exh);
            Th_exh_primary=simDataPSA(index).simData.Data(:,pos_Th_exh);
            Ctotal_primary=simDataPSA(index).simData.Data(:,pos_Ctotal);
            Cx_primary=simDataPSA(index).simData.Data(:,pos_Cx);
            Treg_primary=simDataPSA(index).simData.Data(:,pos_Treg);
            MacM1_primary=simDataPSA(index).simData.Data(:,pos_MacM1);
            MacM2_primary=simDataPSA(index).simData.Data(:,pos_MacM2);
            mAPC_primary=simDataPSA(index).simData.Data(:,pos_mAPC);
            APC_primary=simDataPSA(index).simData.Data(:,pos_APC);
        
            total_teff_primary=Teff_primary+Teff_exh_primary;
            total_th_primary=Th_primary+Th_exh_primary;
            total_tcells_primary=Teff_primary+Th_primary+Treg_primary+Teff_exh_primary+Th_exh_primary;
            total_primary=MDSC_primary+Teff_primary+Th_primary+Teff_exh_primary+Th_exh_primary+Ctotal_primary+Cx_primary...
            +Treg_primary+MacM1_primary+MacM2_primary+mAPC_primary+APC_primary;
            total_cd4_primary=total_th_primary+Treg_primary;
        
            MDSC_met=simDataPSA(index).simData.Data(:,pos_MDSC_met);
            Teff_met=simDataPSA(index).simData.Data(:,pos_Teff_met);
            Th_met=simDataPSA(index).simData.Data(:,pos_Th_met);
            Teff_exh_met=simDataPSA(index).simData.Data(:,pos_Teff_exh_met);
            Th_exh_met=simDataPSA(index).simData.Data(:,pos_Th_exh_met);
            Ctotal_met=simDataPSA(index).simData.Data(:,pos_Ctotal_met);
            Cx_met=simDataPSA(index).simData.Data(:,pos_Cx_met);
            Treg_met=simDataPSA(index).simData.Data(:,pos_Treg_met);
            MacM1_met=simDataPSA(index).simData.Data(:,pos_MacM1_met);
            MacM2_met=simDataPSA(index).simData.Data(:,pos_MacM2_met);
            mAPC_met=simDataPSA(index).simData.Data(:,pos_mAPC_met);
            APC_met=simDataPSA(index).simData.Data(:,pos_APC_met);
        
            total_teff_met=Teff_met+Teff_exh_met;
            total_th_met=Th_met+Th_exh_met;
            total_tcells_met=Teff_met+Th_met+Treg_met+Teff_exh_met+Th_exh_met;
            total_met=MDSC_met+Teff_met+Th_met+Teff_exh_met+Th_exh_met+Ctotal_met+Cx_met...
            +Treg_met+MacM1_met+MacM2_met+mAPC_met+APC_met;
            total_cd4_met=total_th_met+Treg_met;
        
            MDSC_ratio_temp=(MDSC_met./total_met)./(MDSC_primary./total_primary);
            Teff_ratio_temp=(Teff_met./total_met)./(Teff_primary./total_primary);
            Treg_ratio_temp=(Treg_met./total_met)./(Treg_primary./total_primary);
            Th_ratio_temp=(Th_met./total_met)./(Th_primary./total_primary);
            Mac_ratio_temp=((MacM1_met+MacM2_met)./total_met)./((MacM1_primary+MacM2_primary)./total_primary);
            MacM1_ratio_temp=(MacM1_met./total_met)./(MacM1_primary./total_primary);
            MacM2_ratio_temp=(MacM2_met./total_met)./(MacM2_primary./total_primary);
            mAPC_ratio_temp=(mAPC_met./total_met)./(mAPC_primary./total_primary);
            APC_ratio_temp=(APC_met./total_met)./(APC_primary./total_primary);
            allAPC_ratio_temp=((APC_met+mAPC_met)./total_met)./((APC_primary+mAPC_primary)./total_primary);
            total_tcells_ratio_temp=(total_tcells_met./total_met)./(total_tcells_primary./total_primary);
            total_teff_ratio_temp=(total_teff_met./total_met)./(total_teff_primary./total_primary);
            total_th_ratio_temp=(total_th_met./total_met)./(total_th_primary./total_primary);
            total_cd4_ratio_temp=(total_cd4_met./total_met)./(total_cd4_primary./total_primary);


            MDSC_ratio=[MDSC_ratio,MDSC_ratio_temp];
            Teff_ratio=[Teff_ratio,Teff_ratio_temp];
            Treg_ratio=[Treg_ratio,Treg_ratio_temp];
            Th_ratio=[Th_ratio,Th_ratio_temp];
            Mac_ratio=[Mac_ratio,Mac_ratio_temp];
            MacM1_ratio=[MacM1_ratio,MacM1_ratio_temp];
            MacM2_ratio=[MacM2_ratio,MacM2_ratio_temp];
            mAPC_ratio=[mAPC_ratio,mAPC_ratio_temp];
            APC_ratio=[APC_ratio,APC_ratio_temp];
            allAPC_ratio=[allAPC_ratio,allAPC_ratio_temp];
            total_tcells_ratio=[total_tcells_ratio,total_tcells_ratio_temp];
            total_teff_ratio=[total_teff_ratio,total_teff_ratio_temp];
            total_th_ratio=[total_th_ratio,total_th_ratio_temp];
            total_cd4_ratio=[total_cd4_ratio,total_cd4_ratio_temp];

        end
    
    end
    
    writematrix(MDSC_ratio,append(organs(j),'_prim_MDSC.csv'));
    writematrix(Teff_ratio,append(organs(j),'_prim_Teff.csv'));
    writematrix(Treg_ratio,append(organs(j),'_prim_Treg.csv'));
    writematrix(Th_ratio,append(organs(j),'_prim_Th.csv'));
    writematrix(Mac_ratio,append(organs(j),'_prim_Mac.csv'));
    writematrix(MacM1_ratio,append(organs(j),'_prim_MacM1.csv'));
    writematrix(MacM2_ratio,append(organs(j),'_prim_MacM2.csv'));
    writematrix(mAPC_ratio,append(organs(j),'_prim_mAPC.csv'));
    writematrix(APC_ratio,append(organs(j),'_prim_APC.csv'));
    writematrix(allAPC_ratio,append(organs(j),'_prim_allAPC.csv'));
    writematrix(total_tcells_ratio,append(organs(j),'_prim_total_tcells.csv'));
    writematrix(total_teff_ratio,append(organs(j),'_prim_total_teff.csv'));
    writematrix(total_th_ratio,append(organs(j),'_prim_total_th.csv'));
    writematrix(total_cd4_ratio,append(organs(j),'_prim_total_cd4.csv'));

end