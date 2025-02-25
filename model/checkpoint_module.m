% Checkpoint Module
%
% Models PD1 and CTLA4 Interactions
%
% Inputs: model       -- simbio model object with four compartments
%         params      -- object containing the default parameters
%         Tname       -- name of the T cell forming the checkpoint synapse
%         Cname       -- name of the cancer or APC cell forming the
%                        checkpoint synapse
%      TumorArray -- Array of suffixes for tumor comp. names
%      LNArray    -- Array of suffixes for LN comp. names
% Outputs: model -- simbio model object with new PD1 module
%
% Note: This only works for Ti (i >= 1) interaction with APC and Cj. We
% would need to add a check for Tregs if we would like to generalize it to checkpoint
% model of Treg.

function model = checkpoint_module(model,params,Tname,Cname,TumorArray,LNArray)

% number of tumors and LNs
ntumors=length(TumorArray);
uniqLNArray = unique(LNArray);
nLNs = length(uniqLNArray);

% select the right compartment based on cancer or APC
if Cname(1)=='C'
    compDrug = sbioselect(model, 'Name', 'V_T'); 
    gamma = 'gamma_T';
    Param_name = 'C';
    TumorLNArray=TumorArray;
    ntimes=ntumors;
elseif Cname(1)=='A'
    compDrug = sbioselect(model, 'Name', 'V_LN'); 
    gamma = 'gamma_LN';
    Param_name = 'APC';
    TumorLNArray=uniqLNArray;
    ntimes=nLNs;
end

for array=1:ntimes
    % Add the synapse compartment
    comp = addcompartment(model,['syn_',Tname,'_',Cname,TumorLNArray{array}],params.A_syn.Value,'CapacityUnits',params.A_syn.Units);
        set(comp,'Notes',['synapse comparment between ',Tname,' and ',Cname,' ', params.A_syn.Notes]);
end

% Determine if checkpoint sizes have been defined before in antigen module
first_call = true;

try % see if synapse exist
    p = addparameter(model,'A_syn' ,params.A_syn.Value ,'ValueUnits',params.A_syn.Units);
    set(p,'Notes',['Surface area of the synapse ' params.A_syn.Notes]);
catch
    first_call = false;
end

if first_call
    % Add surface areas
    p = addparameter(model,'A_Tcell' ,params.A_Tcell.Value ,'ValueUnits',params.A_Tcell.Units);
        set(p,'Notes',['Surface area of the T cell ' params.A_Tcell.Notes]);
    p = addparameter(model,'A_cell' ,params.A_cell.Value ,'ValueUnits',params.A_cell.Units);
        set(p,'Notes',['Surface area of the Cancer cell ' params.A_cell.Notes]);
    p = addparameter(model,'A_APC' ,params.A_APC.Value ,'ValueUnits',params.A_APC.Units);
        set(p,'Notes',['Surface area of the APC ' params.A_APC.Notes]);
    p = addparameter(model,'d_syn' ,params.d_syn.Value ,'ValueUnits',params.d_syn.Units);
        set(p,'Notes',['The synapse gap distance (kd2D = kd3D*d_syn) ' params.d_syn.Notes]);
end

% Determine if first call
first_call = true;
try % see if synapse exist
    kon = addparameter(model,'kon_PD1_PDL1',params.kon_PD1_PDL1.Value,'ValueUnits',params.kon_PD1_PDL1.Units);
        set(kon,'Notes',['kon of PD1-PDL1 binding ' params.kon_PD1_PDL1.Notes]);
catch
    first_call = false;
end

if first_call
    % Add Pharmacokinetics
    params_aPD1    = pk_parameters('pembrolizumab');
    params_aPDL1   = pk_parameters('atezolizumab');
    params_aCTLA4  = pk_parameters('tremelimumab');
    model = pk_module(model,'aPD1',TumorArray,LNArray,params_aPD1);
    model = pk_module(model,'aPDL1',TumorArray,LNArray,params_aPDL1);
    model = pk_module(model,'aCTLA4',TumorArray,LNArray,params_aCTLA4);
    
    k_out_PDL1 = addparameter(model,'k_out_PDL1',params.k_out_PDL1.Value,'ValueUnits',params.k_out_PDL1.Units); % estimated by 5e4 PMID 29034543
        set(k_out_PDL1,'Notes',['Expression rate of PDL1 on tumor cells ' params.k_out_PDL1.Notes]);
    k_in_PDL1 = addparameter(model,'k_in_PDL1',params.k_in_PDL1.Value,'ValueUnits',params.k_in_PDL1.Units); % 1 PMID: 30442814
        set(k_in_PDL1,'Notes',['Degradation rate of PDL1 on tumor cells ' params.k_in_PDL1.Notes]);
    r_PDL1_IFNg = addparameter(model,'r_PDL1_IFNg',params.r_PDL1_IFNg.Value,'ValueUnits',params.r_PDL1_IFNg.Units); % 4-6
        set(r_PDL1_IFNg,'Notes',['Number of folds increase of PDL1 expression by IFNg ' params.r_PDL1_IFNg.Notes]);
    
    % Add kon Values
    kon = addparameter(model,'kon_PD1_PDL2',params.kon_PD1_PDL2.Value,'ValueUnits',params.kon_PD1_PDL2.Units);
        set(kon,'Notes',['kon of PD1-PDL2 binding ' params.kon_PD1_PDL2.Notes]);
    kon = addparameter(model,'kon_PD1_aPD1',params.kon_PD1_aPD1.Value,'ValueUnits',params.kon_PD1_aPD1.Units);
        set(kon,'Notes',['kon of PD1-aPD1 binding ' params.kon_PD1_aPD1.Notes]);
    kon = addparameter(model,'kon_PDL1_aPDL1',params.kon_PDL1_aPDL1.Value,'ValueUnits',params.kon_PDL1_aPDL1.Units);
        set(kon,'Notes',['kon of PDL1-aPDL1 binding ' params.kon_PDL1_aPDL1.Notes]);
    kon = addparameter(model,'kon_CD28_CD80',params.kon_CD28_CD80.Value,'ValueUnits',params.kon_CD28_CD80.Units);
        set(kon,'Notes',['kon of CD28-CD80 binding ' params.kon_CD28_CD80.Notes]);
    kon = addparameter(model,'kon_CD28_CD86',params.kon_CD28_CD86.Value,'ValueUnits',params.kon_CD28_CD86.Units);
        set(kon,'Notes',['kon of CD28-CD86 binding ' params.kon_CD28_CD86.Notes]);
    kon = addparameter(model,'kon_CTLA4_CD80',params.kon_CTLA4_CD80.Value,'ValueUnits',params.kon_CTLA4_CD80.Units);
        set(kon,'Notes',['kon of CTLA4-CD80 binding ' params.kon_CTLA4_CD80.Notes]);
    kon = addparameter(model,'kon_CTLA4_CD86',params.kon_CTLA4_CD86.Value,'ValueUnits',params.kon_CTLA4_CD86.Units);
        set(kon,'Notes',['kon of CTLA4-CD86 binding ' params.kon_CTLA4_CD86.Notes]);
    kon = addparameter(model,'kon_CD80_PDL1',params.kon_CD80_PDL1.Value,'ValueUnits',params.kon_CD80_PDL1.Units);
        set(kon,'Notes',['kon of CD80-PDL1 binding ' params.kon_CD80_PDL1.Notes]);
    kon = addparameter(model,'kon_CTLA4_aCTLA4',params.kon_CTLA4_aCTLA4.Value,'ValueUnits',params.kon_CTLA4_aCTLA4.Units);
        set(kon,'Notes',['kon of CTLA4-aCTLA4 binding ' params.kon_CTLA4_aCTLA4.Notes]);
    kon = addparameter(model,'kon_CD80_CD80',params.kon_CD80_CD80.Value,'ValueUnits',params.kon_CD80_CD80.Units);
        set(kon,'Notes',['kon of CD80 self-association ' params.kon_CD80_CD80.Notes]);
    
    % Add koff Values
    koff = addparameter(model,'koff_PD1_PDL1' ,params.koff_PD1_PDL1.Value ,'ValueUnits',params.koff_PD1_PDL1.Units);
        set(koff,'Notes',['koff of PD1-PDL1 binding ' params.koff_PD1_PDL1.Notes]);
    koff = addparameter(model,'koff_PD1_PDL2' ,params.koff_PD1_PDL2.Value ,'ValueUnits',params.koff_PD1_PDL2.Units);
        set(koff,'Notes',['koff of PD1-PDL2 binding ' params.koff_PD1_PDL2.Notes]);
    koff = addparameter(model,'koff_PD1_aPD1' ,params.koff_PD1_aPD1.Value ,'ValueUnits',params.koff_PD1_aPD1.Units);
        set(koff,'Notes',['koff of PD1-aPD1 binding ' params.koff_PD1_aPD1.Notes]);
    koff = addparameter(model,'koff_PDL1_aPDL1',params.koff_PDL1_aPDL1.Value,'ValueUnits',params.koff_PDL1_aPDL1.Units);
        set(koff,'Notes',['koff of PDL1-aPDL1 binding ' params.koff_PDL1_aPDL1.Notes]);
    koff = addparameter(model,'koff_CD28_CD80' ,params.koff_CD28_CD80.Value ,'ValueUnits',params.koff_CD28_CD80.Units);
        set(koff,'Notes',['koff of CD28-CD80 binding ' params.koff_CD28_CD80.Notes]);
    koff = addparameter(model,'koff_CD28_CD86' ,params.koff_CD28_CD86.Value ,'ValueUnits',params.koff_CD28_CD86.Units);
        set(koff,'Notes',['koff of CD28-CD86 binding ' params.koff_CD28_CD86.Notes]);
    koff = addparameter(model,'koff_CTLA4_CD80' ,params.koff_CTLA4_CD80.Value ,'ValueUnits',params.koff_CTLA4_CD80.Units);
        set(koff,'Notes',['koff of CTLA4-CD80 binding ' params.koff_CTLA4_CD80.Notes]);
    koff = addparameter(model,'koff_CTLA4_CD86' ,params.koff_CTLA4_CD86.Value ,'ValueUnits',params.koff_CTLA4_CD86.Units);
        set(koff,'Notes',['koff of CTLA4-CD86 binding ' params.koff_CTLA4_CD86.Notes]);
    koff = addparameter(model,'koff_CD80_PDL1' ,params.koff_CD80_PDL1.Value ,'ValueUnits',params.koff_CD80_PDL1.Units);
        set(koff,'Notes',['koff of CD80-PDL1 binding ' params.koff_CD80_PDL1.Notes]);
    koff = addparameter(model,'koff_CTLA4_aCTLA4',params.koff_CTLA4_aCTLA4.Value,'ValueUnits',params.koff_CTLA4_aCTLA4.Units);
        set(koff,'Notes',['koff of CTLA4-aCTLA4 binding ' params.koff_CTLA4_aCTLA4.Notes]);
    koff = addparameter(model,'koff_CD80_CD80',params.koff_CD80_CD80.Value,'ValueUnits',params.koff_CD80_CD80.Units);
        set(koff,'Notes',['koff of CD80 self-association ' params.koff_CD80_CD80.Notes]); 
    
    % Bivalent anibody parameters
    p = addparameter(model,'Chi_PD1_aPD1' ,params.Chi_PD1_aPD1.Value,'ValueUnits',params.Chi_PD1_aPD1.Units);
        set(p,'Notes',['Antibody cross-arm binding efficiency that also includes the conversion of kon from 3D to 2D ' params.Chi_PD1_aPD1.Notes]);
    p = addparameter(model,'Chi_PDL1_aPDL1' ,params.Chi_PDL1_aPDL1.Value,'ValueUnits',params.Chi_PDL1_aPDL1.Units);
        set(p,'Notes',['Antibody cross-arm binding efficiency that also includes the conversion of kon from 3D to 2D ' params.Chi_PDL1_aPDL1.Notes]);
    p = addparameter(model,'Chi_CTLA4_aCTLA4' ,params.Chi_CTLA4_aCTLA4.Value,'ValueUnits',params.Chi_CTLA4_aCTLA4.Units);
        set(p,'Notes',['Antibody cross-arm binding efficiency that also includes the conversion of kon from 3D to 2D ' params.Chi_CTLA4_aCTLA4.Notes]);
    
    % PD1-related Hill parameters
    p = addparameter(model,'PD1_50',params.PD1_50.Value,'ValueUnits',params.PD1_50.Units);
        set(p,'Notes',['PD1/PDL1 concentration for half-maximal T cell inactivation ' params.PD1_50.Notes]);
    p = addparameter(model,'n_PD1',params.n_PD1.Value,'ValueUnits',params.n_PD1.Units);
        set(p,'Notes',['Hill coefficient for PD1/PDL1 half-maximal T cell inactivation ' params.n_PD1.Notes]);
    p = addparameter(model,'CD28_CD8X_50',params.CD28_CD8X_50.Value,'ValueUnits',params.CD28_CD8X_50.Units);
        set(p,'Notes',['CD28-CD80/CD28-CD86 concentration for half-maximal T cell co-estimulation ' params.CD28_CD8X_50.Notes]);
    p = addparameter(model,'n_CD28_CD8X',params.n_CD28_CD8X.Value,'ValueUnits',params.n_CD28_CD8X.Units);
        set(p,'Notes',['Hill coefficient for CD28-CD80/CD28-CD86 half-maximal T cell co-estimulation ' params.n_CD28_CD8X.Notes]);

end

% Checkpoint Expressions
% Check if T cell was defined before
first_Tcell_call = true;
try
    p = addparameter(model,[Tname,'_PD1_total'],params.T8_PD1.Value,'ValueUnits',params.T8_PD1.Units,'ConstantValue',false);
catch
    first_Tcell_call = false;
end
if first_Tcell_call
        set(p,'Notes',['concentration of PD1 on ',Tname,' cells ' params.T8_PD1.Notes]);
    p = addparameter(model,[Tname,'_CD28_total'],params.T8_CD28.Value,'ValueUnits',params.T8_CD28.Units,'ConstantValue',false);
        set(p,'Notes',['concentration of CD28 on ',Tname,' cells ' params.T8_CD28.Notes]);
    p = addparameter(model,[Tname,'_CTLA4_syn'],params.T8_CTLA4.Value,'ValueUnits',params.T8_CTLA4.Units,'ConstantValue',false);
        set(p,'Notes',['concentration of CTLA4 on ',Tname,' cells ' params.T8_CTLA4.Notes]);
    p = addparameter(model,[Tname,'_PDL1_total'],params.T8_PDL1.Value,'ValueUnits',params.T8_PDL1.Units,'ConstantValue',false);
        set(p,'Notes',['concentration of PDL1 on ',Tname,' cells ' params.T8_PDL1.Notes]);
end

% Check if Cancer or APC was defined before
first_Ccell_call = true;
try
    p = addparameter(model,[Cname,'_PDL1_base'],params.([Param_name,'_PDL1']).Value,'ValueUnits',params.([Param_name,'_PDL1']).Units,'ConstantValue',false);
        set(p,'Notes',['baseline number of PDL1 molecules per ',Cname,' cell ' params.([Param_name,'_PDL1']).Notes]);
catch
    first_Ccell_call = false;
end

if first_Ccell_call
    p = addparameter(model,['r_PDL2',Cname],params.(['r_PDL2',Param_name]).Value,'ValueUnits',params.(['r_PDL2',Param_name]).Units,'ConstantValue',false);
        set(p,'Notes',['PDL2/PDL1 molecule ratio on ', Cname,'cell ', params.(['r_PDL2',Param_name]).Notes]);
    p = addparameter(model,[Cname,'_CD80_total'],params.([Param_name,'_CD80']).Value,'ValueUnits',params.([Param_name,'_CD80']).Units,'ConstantValue',false);
        set(p,'Notes',['number of CD80 molecules per ',Cname,' cell ' params.([Param_name,'_CD80']).Notes]);
    p = addparameter(model,[Cname,'_CD86_total'],params.([Param_name,'_CD86']).Value,'ValueUnits',params.([Param_name,'_CD86']).Units,'ConstantValue',false);
        set(p,'Notes',['number of CD86 molecules per ',Cname,' cell ' params.([Param_name,'_CD86']).Notes]);
end

for array=1:ntimes

    % select the compartment
    comp = sbioselect(model,'Name',['syn_',Tname,'_',Cname,TumorLNArray{array}]);
    
    % Add Species
    x = addspecies(comp,'PDL1_total',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of bound and unbound PDL1 molecules ');
    x = addspecies(comp,'PDL2_total',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of bound and unbound PDL2 molecules ');
    x = addspecies(comp,'PD1_PDL1',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of PD1-PDL1 complex');
    x = addspecies(comp,'PD1_PDL2',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of PD1-PDL2 complex');
    x = addspecies(comp,'PD1',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of PD1 in synapse');
    x = addspecies(comp,'PDL1',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of PDL1 in synapse');
    x = addspecies(comp,'PDL2',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of PDL2 in synapse');
    x = addspecies(comp,'PD1_aPD1',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of PD1-aPD1 complex');
    x = addspecies(comp,'PD1_aPD1_PD1',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of PD1-aPD1-PD1 complex');
    x = addspecies(comp,'PDL1_aPDL1',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of PDL1-aPDL1 complex');
    x = addspecies(comp,'PDL1_aPDL1_PDL1',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of PDL1-aPDL1-PDL1 complex');
    x = addspecies(comp,'TPDL1',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of PDL1 in synapse of T cell');
    x = addspecies(comp,'TPDL1_aPDL1',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of TPDL1-aPDL1 complex');
    x = addspecies(comp,'TPDL1_aPDL1_TPDL1',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of TPDL1-aPDL1-TPDL1 complex');
    
    x = addspecies(comp,'CD28_CD80',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CD28-CD80 complex');
    x = addspecies(comp,'CD28_CD80_CD28',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CD28-CD80-CD28 complex');
    x = addspecies(comp,'CD28_CD86',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CD28-CD86 complex');
    x = addspecies(comp,'CD80_CTLA4',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CD80-CTLA4 complex');
    x = addspecies(comp,'CD80_CTLA4_CD80',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CD80-CTLA4-CD80 complex');
    x = addspecies(comp,'CTLA4_CD80_CTLA4',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CTLA4-CD80-CTLA4 complex');
    x = addspecies(comp,'CD80_CTLA4_CD80_CTLA4',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CD80-CTLA4-CD80-CTLA4 complex');
    x = addspecies(comp,'CD86_CTLA4',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CD86-CTLA4 complex');
    x = addspecies(comp,'CD86_CTLA4_CD86',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CD86-CTLA4-CD86 complex');
    x = addspecies(comp,'PDL1_CD80',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of PDL1-CD80 complex');
    x = addspecies(comp,'PDL1_CD80_CD28',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of PDL1-CD80-CD28 complex');
    x = addspecies(comp,'PDL1_CD80_CTLA4',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of PDL1-CD80-CTLA4 complex');
    x = addspecies(comp,'CD28',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CD28 in synapse');
    x = addspecies(comp,'CTLA4',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CTLA4 in synapse');
    x = addspecies(comp,'CD80',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CD80 dimer in synapse');
    x = addspecies(comp,'CD80m',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CD80 monomer in synapse');
    x = addspecies(comp,'CD86',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CD86 in synapse');
    x = addspecies(comp,'CTLA4_aCTLA4',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CTLA4-aCTLA4 complex');
    x = addspecies(comp,'CTLA4_aCTLA4_CTLA4',0,'InitialAmountUnits','molecule/micrometer^2');
        set(x,'Notes','concentration of CTLA4-aCTLA4-CTLA4 complex');
    
    % Update Input Parameters
    addrule(model,[comp.Name,'.PD1',' = '  ,Tname,'_PD1_total /A_Tcell' ] ,'initialAssignment');
    addrule(model,[comp.Name,'.CD28',' = ' ,Tname,'_CD28_total /A_Tcell'] ,'initialAssignment');
    addrule(model,[comp.Name,'.CTLA4',' = ',Tname,'_CTLA4_syn /A_Tcell'   ] ,'initialAssignment');
    addrule(model,[comp.Name,'.TPDL1',' = ',Tname,'_PDL1_total /A_Tcell'] ,'initialAssignment');
    
    if Cname(1)=='C'
        addrule(model,[comp.Name,'.PDL1 = ',Cname,'_PDL1_base /A_cell'] ,'initialAssignment');
        addrule(model,[comp.Name,'.PDL2 = ',Cname,'_PDL1_base*r_PDL2',Cname,' /A_cell'] ,'initialAssignment');
        addrule(model,[comp.Name,'.PDL1_total = ',comp.Name,'.PDL1+',comp.Name,'.PD1_PDL1+',comp.Name,'.PDL1_aPDL1+2*',comp.Name,'.PDL1_aPDL1_PDL1+',...
                                              comp.Name,'.PDL1_CD80+',comp.Name,'.PDL1_CD80_CD28+',comp.Name,'.PDL1_CD80_CTLA4'] ,'repeatedAssignment');
        addrule(model,[comp.Name,'.PDL2_total = ',comp.Name,'.PD1_PDL2+',comp.Name,'.PDL2'] ,'repeatedAssignment');
        addrule(model,[comp.Name,'.CD80 = ',Cname,'_CD80_total /A_cell'] ,'initialAssignment');
        addrule(model,[comp.Name,'.CD86 = ',Cname,'_CD86_total /A_cell'] ,'initialAssignment');
    elseif Cname(1)=='A'
        addrule(model,[comp.Name,'.PDL1 = ',Cname,'_PDL1_base /A_APC'] ,'initialAssignment');
        addrule(model,[comp.Name,'.PDL2 = ',Cname,'_PDL1_base*r_PDL2',Cname,' /A_APC'] ,'initialAssignment');
        addrule(model,[comp.Name,'.PDL1_total = ',comp.Name,'.PDL1+',comp.Name,'.PD1_PDL1+',comp.Name,'.PDL1_aPDL1+2*',comp.Name,'.PDL1_aPDL1_PDL1+',...
                                              comp.Name,'.PDL1_CD80+',comp.Name,'.PDL1_CD80_CD28+',comp.Name,'.PDL1_CD80_CTLA4'] ,'repeatedAssignment');
        addrule(model,[comp.Name,'.PDL2_total = ',comp.Name,'.PD1_PDL2+',comp.Name,'.PDL2'] ,'repeatedAssignment');
        addrule(model,[comp.Name,'.CD80 = ',Cname,'_CD80_total /A_APC'] ,'initialAssignment');
        addrule(model,[comp.Name,'.CD86 = ',Cname,'_CD86_total /A_APC'] ,'initialAssignment');
    end
    
    % PDL1 Secretion
    R = addreaction(model,['null -> ',comp.Name,'.PDL1']);
        set (R, 'ReactionRate', ['k_out_PDL1*V_T' TumorArray{array} '.IFNg/(V_T' TumorArray{array} '.IFNg+IFNg_50_ind)*(1-',...
        comp.Name,'.PDL1_total/(',Cname,'_PDL1_base*r_PDL1_IFNg/A_cell))']);
        set (R, 'Notes'       , 'Translocation of PDL1 between cell surface and cytoplasm');
    R = addreaction(model,['null -> ',comp.Name,'.PDL2']);
        set (R, 'ReactionRate', ['k_out_PDL1*r_PDL2',Cname,'*V_T' TumorArray{array} '.IFNg/(V_T' TumorArray{array} '.IFNg+IFNg_50_ind)*(1-',...
        comp.Name,'.PDL2_total/(',Cname,'_PDL1_base*r_PDL1_IFNg/A_cell*r_PDL2',Cname,'))']);
        set (R, 'Notes'       , 'Translocation of PDL2 between cell surface and cytoplasm');
    R = addreaction(model,['null -> ' comp.Name,'.PDL1']);
        set (R, 'ReactionRate', ['k_in_PDL1*(',Cname,'_PDL1_base/A_cell-',comp.Name,'.PDL1_total)']);
        set (R, 'Notes'       , 'Translocation of PDL1 between cell surface and cytoplasm');
    R = addreaction(model,['null -> ' comp.Name,'.PDL2']);
        set (R, 'ReactionRate', ['k_in_PDL1*(',Cname,'_PDL1_base/A_cell*r_PDL2',Cname,'-',comp.Name,'.PDL2_total)']);
        set (R, 'Notes'       , 'Translocation of PDL2 between cell surface and cytoplasm');
    
    % Dynamics of PD1/PDL1/PDL2/aPD1/aPDL1
     R = addreaction(model,[comp.Name,'.PD1 + ',comp.Name,'.PDL1 <-> ',comp.Name,'.PD1_PDL1']);
        set (R, 'ReactionRate', ['kon_PD1_PDL1*(',comp.Name,'.PD1)*(',comp.Name,'.PDL1)  -  koff_PD1_PDL1*',comp.Name,'.PD1_PDL1']);
        set (R, 'Notes'       , 'binding and unbinding of PD1 PDL1 in synapse');
     R = addreaction(model,[comp.Name,'.PD1 + ',comp.Name,'.PDL2 <-> ',comp.Name,'.PD1_PDL2']);
        set (R, 'ReactionRate', ['kon_PD1_PDL2*(',comp.Name,'.PD1)*(',comp.Name,'.PDL2)  -  koff_PD1_PDL2*',comp.Name,'.PD1_PDL2']);
        set (R, 'Notes'       , 'binding and unbinding of PD1 PDL2 in synapse');
     R = addreaction(model,[comp.Name,'.PD1 <-> ',comp.Name,'.PD1_aPD1']);
        set (R, 'ReactionRate', ['2*kon_PD1_aPD1*(',comp.Name,'.PD1 * ',compDrug.Name,TumorLNArray{array},'.aPD1/',gamma,'_aPD1) -  koff_PD1_aPD1*',comp.Name,'.PD1_aPD1']);
        set (R, 'Notes'       , ['binding and unbinding of PD1 to aPD1 on ',Tname,' surface in synapse']);
     R = addreaction(model,[comp.Name,'.PD1_aPD1 + ',comp.Name,'.PD1 <-> ',comp.Name,'.PD1_aPD1_PD1']);
        set (R, 'ReactionRate', ['Chi_PD1_aPD1*kon_PD1_aPD1*(',comp.Name,'.PD1 * ',comp.Name,'.PD1_aPD1) -  2*koff_PD1_aPD1*',comp.Name,'.PD1_aPD1_PD1']);
        set (R, 'Notes'       , ['binding and unbinding of PD1 to PD1-aPD1 on ',Tname,' surface in synapse']);
     R = addreaction(model,[comp.Name,'.PDL1 <-> ',comp.Name,'.PDL1_aPDL1']);
        set (R, 'ReactionRate', ['2*kon_PDL1_aPDL1*(',comp.Name,'.PDL1 * ',compDrug.Name,TumorLNArray{array},'.aPDL1/',gamma,'_aPDL1) -  koff_PDL1_aPDL1*',comp.Name,'.PDL1_aPDL1']);
        set (R, 'Notes'       , ['binding and unbinding of PDL1 to aPDL1 on ',Cname,' surface in synapse']);
     R = addreaction(model,[comp.Name,'.PDL1_aPDL1 + ',comp.Name,'.PDL1 <-> ',comp.Name,'.PDL1_aPDL1_PDL1']);
        set (R, 'ReactionRate', ['Chi_PDL1_aPDL1*kon_PDL1_aPDL1*(',comp.Name,'.PDL1 * ',comp.Name,'.PDL1_aPDL1) -  2*koff_PDL1_aPDL1*',comp.Name,'.PDL1_aPDL1_PDL1']);
        set (R, 'Notes'       , ['binding and unbinding of PDL1 to PDL1-aPDL1 on ',Cname,' surface in synapse']);
    
    % Dynamics of CD28/CTLA4/CD80/CD86/aCTLA4/aPDL1
    % CD28-CD80
     R = addreaction(model,[comp.Name,'.CD28 + ',comp.Name,'.CD80 <-> ',comp.Name,'.CD28_CD80']);
        set (R, 'ReactionRate', ['2*kon_CD28_CD80*(',comp.Name,'.CD28)*(',comp.Name,'.CD80)  -  koff_CD28_CD80*',comp.Name,'.CD28_CD80']);
        set (R, 'Notes'       , 'binding and unbinding of CD28 and CD80 in synapse');
     R = addreaction(model,[comp.Name,'.CD28_CD80 + ',comp.Name,'.CD28 <-> ',comp.Name,'.CD28_CD80_CD28']);
        set (R, 'ReactionRate', ['kon_CD28_CD80*(',comp.Name,'.CD28)*(',comp.Name,'.CD28_CD80)  -  2*koff_CD28_CD80*',comp.Name,'.CD28_CD80_CD28']);
        set (R, 'Notes'       , 'binding and unbinding of CD28-CD80 and CD28 in synapse');
    % CD28-CD86
     R = addreaction(model,[comp.Name,'.CD28 + ',comp.Name,'.CD86 <-> ',comp.Name,'.CD28_CD86']);
        set (R, 'ReactionRate', ['kon_CD28_CD86*(',comp.Name,'.CD28)*(',comp.Name,'.CD86)  -  koff_CD28_CD86*',comp.Name,'.CD28_CD86']);
        set (R, 'Notes'       , 'binding and unbinding of CD28 and CD86 in synapse');
    % CTLA4-CD80
     R = addreaction(model,[comp.Name,'.CTLA4 + ',comp.Name,'.CD80 <-> ',comp.Name,'.CD80_CTLA4']);
        set (R, 'ReactionRate', ['4*kon_CTLA4_CD80*(',comp.Name,'.CTLA4)*(',comp.Name,'.CD80)  -  koff_CTLA4_CD80*',comp.Name,'.CD80_CTLA4']);
        set (R, 'Notes'       , 'binding and unbinding of CTLA4 and CD80 in synapse');
     R = addreaction(model,[comp.Name,'.CTLA4 + ',comp.Name,'.CD80_CTLA4 <-> ',comp.Name,'.CTLA4_CD80_CTLA4']);
         set (R, 'ReactionRate', ['2*kon_CTLA4_CD80*(',comp.Name,'.CTLA4)*(',comp.Name,'.CD80_CTLA4)  -  2*koff_CTLA4_CD80*',comp.Name,'.CTLA4_CD80_CTLA4']);
        set (R, 'Notes'       , 'binding and unbinding of CTLA4 and CD80-CTLA4 in synapse');
     R = addreaction(model,[comp.Name,'.CD80 + ',comp.Name,'.CTLA4_CD80_CTLA4 <-> ',comp.Name,'.CD80_CTLA4_CD80_CTLA4']);
        set (R, 'ReactionRate', ['4*kon_CTLA4_CD80*(',comp.Name,'.CD80)*(',comp.Name,'.CTLA4_CD80_CTLA4)  -  koff_CTLA4_CD80*',comp.Name,'.CD80_CTLA4_CD80_CTLA4']);
        set (R, 'Notes'       , 'binding and unbinding of CD80 and CTLA4-CD80-CTLA4 in synapse');
     R = addreaction(model,[comp.Name,'.CD80_CTLA4 + ',comp.Name,'.CD80 <-> ',comp.Name,'.CD80_CTLA4_CD80']);
        set (R, 'ReactionRate', ['2*kon_CTLA4_CD80*(',comp.Name,'.CD80_CTLA4)*(',comp.Name,'.CD80)  -  2*koff_CTLA4_CD80*',comp.Name,'.CD80_CTLA4_CD80']);
        set (R, 'Notes'       , 'binding and unbinding of CD80-CTLA4 and CD80 in synapse');
     R = addreaction(model,[comp.Name,'.CTLA4 + ',comp.Name,'.CD80_CTLA4_CD80 <-> ',comp.Name,'.CD80_CTLA4_CD80_CTLA4']);
        set (R, 'ReactionRate', ['4*kon_CTLA4_CD80*(',comp.Name,'.CTLA4)*(',comp.Name,'.CD80_CTLA4_CD80)  -  koff_CTLA4_CD80*',comp.Name,'.CD80_CTLA4_CD80_CTLA4']);
        set (R, 'Notes'       , 'binding and unbinding of CTLA4 and CD80-CTLA4-CD80 in synapse');
    % CTLA4-CD86
     R = addreaction(model,[comp.Name,'.CTLA4 + ',comp.Name,'.CD86 <-> ',comp.Name,'.CD86_CTLA4']);
        set (R, 'ReactionRate', ['2*kon_CTLA4_CD86*(',comp.Name,'.CTLA4)*(',comp.Name,'.CD86)  -  koff_CTLA4_CD86*',comp.Name,'.CD86_CTLA4']);
        set (R, 'Notes'       , 'binding and unbinding of CTLA4 and CD86 in synapse');
     R = addreaction(model,[comp.Name,'.CD86_CTLA4 + ',comp.Name,'.CD86 <-> ',comp.Name,'.CD86_CTLA4_CD86']);
        set (R, 'ReactionRate', ['kon_CTLA4_CD86*(',comp.Name,'.CD86_CTLA4)*(',comp.Name,'.CD86)  -  2*koff_CTLA4_CD86*',comp.Name,'.CD86_CTLA4_CD86']);
        set (R, 'Notes'       , 'binding and unbinding of CD86-CTLA4 and CD86 in synapse');
    % CTLA4-aCTLA4
     R = addreaction(model,[comp.Name,'.CTLA4 <-> ',comp.Name,'.CTLA4_aCTLA4']);
        set (R, 'ReactionRate', ['4*kon_CTLA4_aCTLA4*(',comp.Name,'.CTLA4 * ',compDrug.Name,TumorLNArray{array},'.aCTLA4/',gamma,'_aCTLA4) -  koff_CTLA4_aCTLA4*',comp.Name,'.CTLA4_aCTLA4']);
        set (R, 'Notes'       , ['binding and unbinding of CTLA4 to aCTLA4 on ',Tname,' surface in synapse']);
     R = addreaction(model,[comp.Name,'.CTLA4_aCTLA4 + ',comp.Name,'.CTLA4 <-> ',comp.Name,'.CTLA4_aCTLA4_CTLA4']);
        set (R, 'ReactionRate', ['2*Chi_CTLA4_aCTLA4*kon_CTLA4_aCTLA4*(',comp.Name,'.CTLA4 * ',comp.Name,'.CTLA4_aCTLA4) -  2*koff_CTLA4_aCTLA4*',comp.Name,'.CTLA4_aCTLA4_CTLA4']);
        set (R, 'Notes'       , ['binding and unbinding of CTLA4 to aCTLA4 on ',Tname,' surface in synapse']);
    % CD80 dimer dissociation
     R = addreaction(model,[comp.Name,'.CD80m + ',comp.Name,'.CD80m <-> ',comp.Name,'.CD80']);
        set (R, 'ReactionRate', ['kon_CD80_CD80*(',comp.Name,'.CD80m)*(',comp.Name,'.CD80m) - koff_CD80_CD80*',comp.Name,'.CD80']);
        set (R, 'Notes'       , 'self-association and dissociation of CD80 monomers in synapse');
    % cis PDL1-CD80-CD28
     R = addreaction(model,[comp.Name,'.CD80m + ',comp.Name,'.PDL1 <-> ',comp.Name,'.PDL1_CD80']);
        set (R, 'ReactionRate', ['kon_CD80_PDL1*(',comp.Name,'.CD80m)*(',comp.Name,'.PDL1)  -  koff_CD80_PDL1*',comp.Name,'.PDL1_CD80']);
        set (R, 'Notes'       , 'binding and unbinding of CD80 and PDL1 in synapse');
     R = addreaction(model,[comp.Name,'.PDL1_CD80 + ',comp.Name,'.CD28 <-> ',comp.Name,'.PDL1_CD80_CD28']);
        set (R, 'ReactionRate', ['kon_CD28_CD80*(',comp.Name,'.PDL1_CD80)*(',comp.Name,'.CD28)  - koff_CD28_CD80*',comp.Name,'.PDL1_CD80_CD28']);
        set (R, 'Notes'       , 'binding and unbinding of PDL1-CD80 and CD28 in synapse');
     R = addreaction(model,[comp.Name,'.PDL1_CD80 + ',comp.Name,'.CTLA4 <-> ',comp.Name,'.PDL1_CD80_CTLA4']);
        set (R, 'ReactionRate', ['2*kon_CTLA4_CD80*(',comp.Name,'.PDL1_CD80)*(',comp.Name,'.CTLA4)  - koff_CTLA4_CD80*',comp.Name,'.PDL1_CD80_CTLA4']);
        set (R, 'Notes'       , 'binding and unbinding of PDL1-CD80 and CTLA4 in synapse'); 
    % TPDL1-aPDL1
     R = addreaction(model,[comp.Name,'.TPDL1 <-> ',comp.Name,'.TPDL1_aPDL1']);
        set (R, 'ReactionRate', ['2*kon_PDL1_aPDL1*(',comp.Name,'.TPDL1 * ',compDrug.Name,TumorLNArray{array},'.aPDL1/',gamma,'_aPDL1) -  koff_PDL1_aPDL1*',comp.Name,'.TPDL1_aPDL1']);
        set (R, 'Notes'       , ['binding and unbinding of PDL1 to aPDL1 on ',Cname,' surface in synapse']);
     R = addreaction(model,[comp.Name,'.TPDL1_aPDL1 + ',comp.Name,'.TPDL1 <-> ',comp.Name,'.TPDL1_aPDL1_TPDL1']);
        set (R, 'ReactionRate', ['Chi_PDL1_aPDL1*kon_PDL1_aPDL1*(',comp.Name,'.TPDL1 * ',comp.Name,'.TPDL1_aPDL1) -  2*koff_PDL1_aPDL1*',comp.Name,'.TPDL1_aPDL1_TPDL1']);
        set (R, 'Notes'       , ['binding and unbinding of PDL1 to PDL1-aPDL1 on ',Cname,' surface in synapse']);
    
    % Update PD1 Hill Function
    addrule(model,['H_PD1_',Cname,TumorLNArray{array},' = ((',comp.Name,'.PD1_PDL1+',comp.Name,'.PD1_PDL2)/PD1_50)^n_PD1/(((',comp.Name,'.PD1_PDL1+',comp.Name,'.PD1_PDL2)/PD1_50)^n_PD1 + 1)'],'repeatedAssignment');
    addrule(model,['H_CD28_',Cname,TumorLNArray{array},' = ((',comp.Name,'.CD28_CD80 + ',comp.Name,'.CD28_CD86 + 2*',comp.Name,'.CD28_CD80_CD28 + ',comp.Name,'.PDL1_CD80_CD28)/CD28_CD8X_50)^n_CD28_CD8X/(((',comp.Name,'.CD28_CD80 + ',comp.Name,'.CD28_CD86 + 2*',comp.Name,'.CD28_CD80_CD28 + ',comp.Name,'.PDL1_CD80_CD28)/CD28_CD8X_50)^n_CD28_CD8X + 1)'],'repeatedAssignment');

end