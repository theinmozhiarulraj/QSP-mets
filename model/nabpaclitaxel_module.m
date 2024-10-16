% Nab-Paclitaxel Module
%
% Models nab-paclitaxel PK-PD
%
% Inputs: model        -- SimBiology model object
%         params       -- object containing the default parameters
%         TumorArray   -- Array of suffixes for tumor comp. names
% Outputs: model -- SimBiology model object with new nab-paclitaxel module

% thein: not yet adapted for multiple cancer clones

function model = nabpaclitaxel_module(model, params,TumorArray,LNArray)

ntumors=length(TumorArray);

% Setup Compartments
comp_V1 = addcompartment(model,'V_1',15.8,'CapacityUnits','liter','ConstantCapacity',false);
    set(comp_V1,'Notes',['Central compartment (V1) for nab-paclitaxel PK']);
comp_V2 = addcompartment(model,'V_2',1650,'CapacityUnits','liter','ConstantCapacity',false);
    set (comp_V2,'Notes',['Peripheral compartment (V2) for nab-paclitaxel PK']);
comp_V3 = addcompartment(model,'V_3',75.4,'CapacityUnits','liter','ConstantCapacity',false);
    set (comp_V3,'Notes',['Peripheral compartment (V3) for nab-paclitaxel PK']);

% Model Parameters
Vmcl = addparameter(model,'Vmcl',params.Vmcl.Value,'ValueUnits',params.Vmcl.Units);
    set(Vmcl,'Notes','Maximum elimination rate of nab-paclitaxel from the central compartment V1');
Kcl = addparameter(model,'Kcl',params.Kcl.Value,'ValueUnits',params.Kcl.Units);
    set(Kcl,'Notes','Nab-paclitaxel concentration in the central compartment V1 at 50% Vmcl');
Q2 = addparameter(model,'Q2',params.Q2.Value,'ValueUnits',params.Q2.Units);
    set(Q2,'Notes','Intercompartmental clearance between the central compartment V1 and the second peripheral compartment V3');
Vmt = addparameter(model,'Vmt',params.Vmt.Value,'ValueUnits',params.Vmt.Units);
    set(Vmt,'Notes','Maximum intercompartmental distribution rate between the central compartment V1 and the first peripheral compartment V2');
Kt = addparameter(model,'Kt',params.Kt.Value,'ValueUnits',params.Kt.Units);
    set(Kt,'Notes','Nab-paclitaxel concentration in the central compartment V1 at 50% Vmt');
BSA = addparameter(model,'BSA',params.BSA.Value,'ValueUnits',params.BSA.Units);
    set(BSA,'Notes','Body surface area');

% Model Species
s = addspecies(comp_V1,'NabP',0,'InitialAmountUnits','nanogram/milliliter');
    set(s,'Notes','Nab-paclitaxel in the central compartment V1');
s = addspecies(comp_V2,'NabP',0,'InitialAmountUnits','nanogram/milliliter');
    set(s,'Notes','Nab-paclitaxel in the peripheral compartment V2');
s = addspecies(comp_V3,'NabP',0,'InitialAmountUnits','nanogram/milliliter');
    set(s,'Notes','Nab-paclitaxel in the peripheral compartment V3');

% Model Reaction
r = addreaction(model,'V_1.NabP <-> V_2.NabP');
    set(r,'ReactionRate','Vmt/(V_1.NabP+Kt)*V_1.NabP - Vmt/(V_2.NabP+Kt)*V_2.NabP');
    set(r,'Notes','Intercompartmental distribution of nab-paclitaxel between V1 and V2 compartment');
r = addreaction(model,'V_1.NabP <-> V_3.NabP');
    set(r,'ReactionRate','Q2*V_1.NabP - Q2*V_3.NabP');
    set(r,'Notes','Intercompartmental clearance of nab-paclitaxel between V1 and V3 compartment');
r = addreaction(model,'V_1.NabP -> null');
    set(r,'ReactionRate','Vmcl/(V_1.NabP+Kcl)*V_1.NabP');
    set(r,'Notes','Clearance of nab-paclitaxel from V1 compartment');

% PD Cytotoxicity
r_nabp = addparameter(model,'r_nabp',params.r_nabp.Value,'ValueUnits',params.r_nabp.Units,'ConstantValue',false);
    set(r_nabp,'Notes',['Tumour to plasma concentration ratio of nab-paclitaxel ' params.r_nabp.Notes]);
k_C_nabp = addparameter(model,'k_C_nabp',params.k_C_nabp.Value,'ValueUnits',params.k_C_nabp.Units,'ConstantValue',false);
    set(k_C_nabp,'Notes',['Cancer cell killing rate by nab-paclitaxel ' params.k_C_nabp.Notes]);
IC50_nabp = addparameter(model,'IC50_nabp',params.IC50_nabp.Value,'ValueUnits',params.IC50_nabp.Units,'ConstantValue',false);
    set(IC50_nabp,'Notes',['Half-maximal nab-paclitaxel concentration for cancer cell killing ' params.IC50_nabp.Notes]);
MW_nabp = addparameter(model,'MW_nabp',params.MW_nabp.Value,'ValueUnits',params.MW_nabp.Units);
    set(MW_nabp,'Notes',['Molecular weight of nab-paclitaxel ' params.MW_nabp.Notes]);
Kc_nabp = addparameter(model,'Kc_nabp',params.Kc_nabp.Value,'ValueUnits',params.Kc_nabp.Units);
    set(Kc_nabp,'Notes',['Half-Maximal cancer cell number for cytotoxic drug diffusion ' params.Kc_nabp.Notes]);

for array=1:ntumors

    NabP = addspecies(sbioselect(model, 'Name', ['V_T' TumorArray{array}]),'NabP',0,'InitialAmountUnits','nanomolarity'); 
        set(NabP,'Notes','Nab-paclitaxel concentration in tumour ');
    addrule(model,['V_T' TumorArray{array} '.NabP = r_nabp*V_1.NabP/MW_nabp'],'repeatedAssignment');
    
    reaction = addreaction(model,['V_T' TumorArray{array} '.C1 -> V_T' TumorArray{array} '.C_x']);
        set(reaction,'ReactionRate',['k_C_nabp*V_T' TumorArray{array} '.C1*(V_T' TumorArray{array} '.NabP/(V_T' TumorArray{array} '.NabP+IC50_nabp))*min(C_total' TumorArray{array} ',Kc_nabp)/C_total' TumorArray{array}]); % *(1-C_total/(C_total+Kc_nabp))
        set(reaction,'Notes','Cancer cell death by nab-paclitaxel ');  
    addrule(model,['k_C1_therapy' TumorArray{array} ' = k_C_nabp*(V_T' TumorArray{array} '.NabP/(V_T' TumorArray{array} '.NabP+IC50_nabp))*min(C_total' TumorArray{array} ',Kc_nabp)/C_total' TumorArray{array}],'repeatedAssignment');

end

%% Resistance
k_C_resist = addparameter(model,'k_C_resist',params.k_C_resist.Value,'ValueUnits',params.k_C_resist.Units,'ConstantValue',false);
    set(k_C_resist,'Notes',['Cancer resistance to nab-paclitaxel ' params.k_C_resist.Notes]);
r_resist = addparameter(model,'r_resist',params.r_resist.Value,'ValueUnits',params.r_resist.Units);
    set(r_resist,'Notes',['Number of folds increase of nab-paclitaxel EC50 in resistant cancer clones ' params.r_resist.Notes]);

for array=1:ntumors

    reaction = addreaction(model,['V_T' TumorArray{array} '.C1 -> V_T' TumorArray{array} '.C2']);
        set(reaction,'ReactionRate',['k_C_resist*V_T' TumorArray{array} '.C1*H_TGFb']); % *H_TGF_CTL
        set(reaction,'Notes','Cancer cell resistance to nab-paclitaxel ');
    reaction = addreaction(model,['V_T' TumorArray{array} '.C2 -> V_T' TumorArray{array} '.C_x']);
        set(reaction,'ReactionRate',['k_C_nabp*V_T' TumorArray{array} '.C2*(V_T' TumorArray{array} '.NabP/(V_T' TumorArray{array} '.NabP+IC50_nabp*r_resist))*min(C_total' TumorArray{array} ',Kc_nabp)/C_total' TumorArray{array}]);
        set(reaction,'Notes','Resistant cancer cell death by nab-paclitaxel ');
    
    addrule(model,['k_C2_therapy' TumorArray{array} ' = k_C_nabp*(V_T' TumorArray{array} '.NabP/(V_T' TumorArray{array} '.NabP+IC50_nabp*r_resist))*min(C_total' TumorArray{array} ',Kc_nabp)/C_total' TumorArray{array}],'repeatedAssignment');

end

% set tumour growth rate of the resistant clone same as the sensitive clone
% thein: commented these lines uncomment when nabp module is used
addrule(model,'k_C2_growth = k_C1_growth','repeatedAssignment');
addrule(model,'k_C2_death = k_C1_death','repeatedAssignment');

%% Vascularization
p = addparameter(model,'k_vas_nabp',params.k_vas_nabp.Value,'ValueUnits',params.k_vas_nabp.Units,'ConstantValue',false);
    set(p,'Notes',['Secretion rate of angiogenic factors induced by nab-paclitaxel ' params.k_vas_nabp.Notes]);
p = addparameter(model,'IC50_nabp_vas',params.IC50_nabp_vas.Value,'ValueUnits',params.IC50_nabp_vas.Units,'ConstantValue',false);
    set(p,'Notes',['Half-maximal conc. of nab-paclitaxel on angiogenic factor induction ' params.IC50_nabp_vas.Notes]);
p = addparameter(model,'k_K_nabp',params.k_K_nabp.Value,'ValueUnits',params.k_K_nabp.Units,'ConstantValue',false);
    set(p,'Notes',['Inhibition rate of maximal tumor capacity by nab-paclitaxel ' params.k_K_nabp.Notes]);

for array=1:ntumors

    reaction = addreaction(model,['null -> V_T' TumorArray{array} '.c_vas']);
        set(reaction,'ReactionRate',['k_vas_nabp*C_total' TumorArray{array} '*V_T' TumorArray{array} '.NabP/(V_T' TumorArray{array} '.NabP+IC50_nabp_vas)']);
        set(reaction,'Notes','Angiogenic factor release in response to nab-paclitaxel ');
    reaction = addreaction(model,['V_T' TumorArray{array} '.K -> null']);
        set(reaction,'ReactionRate',['k_K_nabp*V_T' TumorArray{array} '.K*V_T' TumorArray{array} '.NabP']);
        set(reaction,'Notes','Inhibition of tumor vasculature due to endothelial cell death by nab-paclitaxel ');

end