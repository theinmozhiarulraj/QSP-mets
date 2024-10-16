% Antigen Module
%
% Models antigen presentation on APCs
%
% Requirements: cancer_module, Tcell_module and APC_module
%
% Inputs: model   -- SimBiology model object with four compartments
%         ID      -- T cell-antigen ID number [must be unique]
%         params  -- object containing model parameter Values, Units, and Notes:
%                    - N_MHC--number of types of MHC molecules
%                    - MHC_T--total amount of MHC per area
%                    - kin--rate of MHC internalization
%                    - kout--rate of MHC exernalization
%                    - V_e--endosomal volume
%                    - A_e--endosomal surface area
%                    - A_s--APC surface area
%                    - k_up--rate of antigen uptake by APCs
%                    - k_xP_deg--rate of extracellular antigen degradation
%                    - k_P_deg--rate of endosomal antigen degradation
%                    - k_p_deg--rate of endosomal epitope degradation
%                    - k_on--rate of antigen_MHC binding
%                    - p_50--epitope concentration for half-maximal T cell activation
%         antigen -- antigen structure
%         TumorArray -- Array of suffixes for tumor comp. names
%         LNArray    -- Array of suffixes for LN comp. names
%
% Outputs: model -- SimBiology model object with new antigen module


function model = antigen_module(model,ID,params,antigen,TumorArray,LNArray)

% number of tumors
ntumors=length(TumorArray);

% Names
antigen_name = ['P' ID];
epitope_name = ['p' ID];
Tcell_name = ['T' ID];

% Get Number of T Cell Clones
nTcells = howManyClones(model); %thein: this is not used

% Add MHCs on First Call
N_MHC = params.N_MHC.Value;
if (length(model.Compartment)<(2+length(unique(LNArray))+ntumors+1)) 
    first_call = true;
else
    first_call = false;
end

if (first_call)
    for array=1:ntumors
        % Add Endosomal And Surface Compartments 
        comp = addcompartment(model,['V_e' TumorArray{array}],params.V_e.Value,'CapacityUnits',params.V_e.Units);
            set(comp,'Notes',['APC endosomal comparment ' params.V_e.Notes]);
        comp = addcompartment(model,['A_e' TumorArray{array}],params.A_e.Value,'CapacityUnits',params.A_e.Units);
            set(comp,'Notes',['APC endosomal surface comparment ' params.A_e.Notes]);
        comp = addcompartment(model,['A_s' TumorArray{array}],params.A_s.Value,'CapacityUnits',params.A_s.Units);
            set(comp,'Notes',['APC surface comparment ' params.A_s.Notes]);
    end

    % Add Translocation Rates
    kin = addparameter(model,'kin',params.kin.Value,'ValueUnits',params.kin.Units);
        set(kin,'Notes',['Rate of MHC internalization ' params.kin.Notes]);
    kout = addparameter(model,'kout',params.kout.Value,'ValueUnits',params.kout.Units);
        set(kout,'Notes',['Rate of MHC externalization ' params.kout.Notes]);

    % Add MHCs
    for i = 1:N_MHC
        model = MHC_module(model,i,params.MHC_T,TumorArray,LNArray);
    end
end

% Add Parameters
k_up = addparameter(model,'k_up',params.k_up.Value,'ValueUnits',params.k_up.Units);
    set(k_up,'Notes',['Rate of antigen uptake by APCs ' params.k_up.Notes]);
k_xP_deg = addparameter(model,'k_xP_deg', params.k_xP_deg.Value,'ValueUnits',params.k_xP_deg.Units);
    set(k_xP_deg,'Notes',['Rate of extracellular antigen degradation ' params.k_xP_deg.Notes]);
k_P_deg = addparameter(model,'k_P_deg',params.k_P_deg.Value,'ValueUnits',params.k_P_deg.Units);
    set(k_P_deg,'Notes',['Rate of endosomal antigen degradation ' params.k_P_deg.Notes]);
k_p_deg = addparameter(model,'k_p_deg',params.k_p_deg.Value,'ValueUnits',params.k_p_deg.Units);
    set(k_p_deg,'Notes',['Rate of endosomal epitope degradation ' params.k_p_deg.Notes]);
k_on = addparameter(model,'k_on',params.k_on.Value,'ValueUnits',params.k_on.Units);
    set(k_on,'Notes',['Rate of antigen-MHC binding ' params.k_on.Notes]);
for i = 1:N_MHC
    k_d = addparameter(model,['k_' antigen_name '_d' num2str(i)],antigen.kd(i).Value,'ValueUnits',antigen.kd(i).Units);
        set(k_d,'Notes',['Antigen-MHC kd ' antigen.kd(i).Notes]);
end
p_50 = addparameter(model,[epitope_name '_50'],params.p_50.Value,'ValueUnits',params.p_50.Units);
    set(p_50,'Notes',params.p_50.Notes);

% Determine if checkpoint sizes have been defined
first_call = true;
try % see if synapse exist
    parameter = addparameter(model,'A_syn' ,params.A_syn.Value ,'ValueUnits',params.A_syn.Units);
    set(parameter,'Notes',['Surface area of the synapse ' params.A_syn.Notes]);
catch
    first_call = false;
end

if first_call
    % Add surface areas
    parameter = addparameter(model,'A_Tcell' ,params.A_Tcell.Value ,'ValueUnits',params.A_Tcell.Units);
        set(parameter,'Notes',['Surface area of the T cell ' params.A_Tcell.Notes]);
    parameter = addparameter(model,'A_cell' ,params.A_cell.Value ,'ValueUnits',params.A_cell.Units);
        set(parameter,'Notes',['Surface area of the Cancer cell ' params.A_cell.Notes]);
    parameter = addparameter(model,'A_APC' ,params.A_APC.Value ,'ValueUnits',params.A_APC.Units);
        set(parameter,'Notes',['Surface area of the APC ' params.A_APC.Notes]);
end

for array=1:ntumors

    k_dep = antigen_rate(model,ID,antigen.cancers,antigen.concentration,nTcells,TumorArray{array});

    % thein: changed initial values of V_T.P, V_e.P & V_e.p from 1e-18 to 0
    % Add Antigen and Epitope
    Pf = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'P',0,'InitialAmountUnits','molarity'); 
        set(Pf,'Notes',['Concentration of  free antigen (' antigen_name ') in the LN compartment']);
    Pe = addspecies(sbioselect(model, 'Name', ['V_e' TumorArray{array}]),'P',0,'InitialAmountUnits','molarity'); 
        set(Pe,'Notes',['Concentration of antigen ' antigen_name ' in the APC endosomes']);
    p = addspecies(sbioselect(model, 'Name', ['V_e' TumorArray{array}]),'p',0,'InitialAmountUnits','molarity'); 
        set(p,'Notes',['Concentration of epitope ' antigen_name ' in the APC endosomes']);
    
    % Add Reactions
    % Antigen
    % Antigen Source
    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.P']);
        set(reaction,'ReactionRate',k_dep);
        set(reaction,'Notes','Antigen deposition from dying cancer cells');
    
    % Free Antigen Degradation
    reaction = addreaction(model,['V_T' TumorArray{array} '.P -> null']);
        set(reaction,'ReactionRate',['k_xP_deg*V_T' TumorArray{array} '.P*V_T' TumorArray{array}]);
        set(reaction,'Notes','Free antigen degradation');
    
    % Antigen Uptake by mAPCs
    reaction = addreaction(model,['V_T' TumorArray{array} '.P -> null']);
        set(reaction,'ReactionRate',['k_up*V_T' TumorArray{array} '.APC*V_T' TumorArray{array} '.P*V_T' TumorArray{array}]);
        set(reaction,'Notes','Antigen uptake by mature antigen presenting cells');
    
    reaction = addreaction(model,['null -> V_e' TumorArray{array} '.P']);
        set(reaction,'ReactionRate',['k_up*cell*V_T' TumorArray{array} '.P*V_e' TumorArray{array}]);
        set(reaction,'Notes','Antigen uptake by mature antigen presenting cells');
    
    % Antigen Degradation in APC Endosomes
    reaction = addreaction(model,['V_e' TumorArray{array} '.P -> V_e' TumorArray{array} '.p']);
        set(reaction,'ReactionRate',['k_P_deg*V_e' TumorArray{array} '.P*V_e' TumorArray{array}]);
        set(reaction,'Notes','Antigen degradation in APC endosomes');

    % Epitope
    reaction = addreaction(model,['V_e' TumorArray{array} '.p -> null']);
        set(reaction,'ReactionRate',['k_p_deg*V_e' TumorArray{array} '.p*V_e' TumorArray{array}]);
        set(reaction,'Notes','Epitope degradation in APC endosomes');

    % MHC-Epitope Complexes
    for i = 1:N_MHC
        %thein: changed the initial values to 0 from 1e-6
        Mp_e = addspecies(sbioselect(model, 'Name', ['A_e' TumorArray{array}]),'Mp',0,'InitialAmountUnits',params.MHC_T.Units); 
            set(Mp_e,'Notes','Antigen-MHC complex');
        Mp_s = addspecies(sbioselect(model, 'Name', ['A_s' TumorArray{array}]),'Mp',0,'InitialAmountUnits',params.MHC_T.Units); 
            set(Mp_s,'Notes','Antigen-MHC complex');
        % Binding in Endosome
        reaction = addreaction(model,['V_e' TumorArray{array} '.p + A_e' TumorArray{array} '.M' num2str(i) ' -> A_e' TumorArray{array} '.Mp']);
            set(reaction,'ReactionRate',['k_on*V_e' TumorArray{array} '.p*A_e' TumorArray{array} '.M' num2str(i) '*A_e' TumorArray{array}]);
            set(reaction,'Notes','Antigen-MHC binding in endosome');
        % Unbinding in Endosome
        reaction = addreaction(model,['A_e' TumorArray{array} '.Mp -> V_e' TumorArray{array} '.p + A_e' TumorArray{array} '.M' num2str(i)]);
            set(reaction,'ReactionRate',['k_' antigen_name '_d' num2str(i) '*k_on*A_e' TumorArray{array} '.Mp*A_e' TumorArray{array}]);
            set(reaction,'Notes','Antigen-MHC unbinding in endosome');
        % Unbinding on Surface
        reaction = addreaction(model,['A_s' TumorArray{array} '.Mp -> A_s' TumorArray{array} '.M' num2str(i)]);
            set(reaction,'ReactionRate',['k_' antigen_name '_d' num2str(i) '*k_on*A_s' TumorArray{array} '.Mp*A_s' TumorArray{array}]);
            set(reaction,'Notes','Antigen-MHC unbinding on APC surface');
        % Antigen-MHC Translocation
        reaction = addreaction(model,['A_e' TumorArray{array} '.Mp -> A_s' TumorArray{array} '.Mp']);
            set(reaction,'ReactionRate',['kout*A_e' TumorArray{array} '.Mp*A_e' TumorArray{array}]);
            set(reaction,'Notes','Antigen-MHC translocation');

        comp_e = sbioselect(model,'Name',['A_e' TumorArray{array}]);
        comp_s = sbioselect(model,'Name',['A_s' TumorArray{array}]);
        rename(sbioselect(comp_e,'Name','Mp'),['M' num2str(i) epitope_name]);
        rename(sbioselect(comp_s,'Name','Mp'),['M' num2str(i) epitope_name]);
    end

end

% Update Hill Function for Antigens
if ID == '0'
    % Adds TCR kinetic proofreading or receptor occupancy to the model for Treg
    %model = TCR_RO_module(model,ID,params,TumorArray);
    model = TCR_KPR_module(model,ID,params,TumorArray);
else
    % Adds TCR kinetic proofreading to the model for Teff
    %(it only works with 1 MHC at the moment)
    model = TCR_KPR_module(model,ID,params,TumorArray);
    % model = TCR_RO_module(model,ID,params);
end


for array=1:ntumors
    % Rename Objects with 'species_name'
    comp_T=sbioselect(model,'Name',['V_T' TumorArray{array}]);
    comp_e=sbioselect(model,'Name',['V_e' TumorArray{array}]);
    rename(sbioselect(comp_T,'Name','P'),antigen_name);
    rename(sbioselect(comp_e,'Name','P'),antigen_name);
    rename(sbioselect(comp_e,'Name','p'),epitope_name);  
end
rename(k_up,['k_' antigen_name '_up']);
rename(k_xP_deg,['k_x' antigen_name '_deg']);
rename(k_P_deg,['k_' antigen_name '_deg']);
rename(k_p_deg,['k_' epitope_name '_deg']);
rename(k_on,['k_' antigen_name '_on']);

