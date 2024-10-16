% MHC Module
%
% Sub-module used by antigen module
%
% Inputs: model -- SimBiology model object with four compartments
%         ID    -- ID number for MHC [must be unique]
%         ICs   -- array containing fraction of MHC in
%                  - edosomal compartment
%                  - surface compartment
%         MHC_T -- total amount of MHC [parameter]
%    TumorArray -- Array of suffixes for tumor comp. names
% Outputs: model -- SimBiology model object with new MHC module

function model = MHC_module(model,ID,MHC_T,TumorArray,LNArray)

% number of tumor compartments
ntumors=length(TumorArray);

for array=1:ntumors

    % Add MHC
    M_e = addspecies(sbioselect(model, 'Name', ['A_e' TumorArray{array}]),'M',MHC_T.Value,'InitialAmountUnits',MHC_T.Units); 
        set(M_e,'Notes','Amount of MHC per area on the endosomal surface');
    M_s = addspecies(sbioselect(model, 'Name', ['A_s' TumorArray{array}]),'M',1e-6,'InitialAmountUnits',MHC_T.Units); 
        set(M_s,'Notes','Amount of MHC per area on the cell surface');
    
    % Add Reactions
    reaction = addreaction(model,['A_e' TumorArray{array} '.M -> A_s' TumorArray{array} '.M']);
        set(reaction,'ReactionRate',['kout*A_e' TumorArray{array} '.M*A_e' TumorArray{array} '-kin*A_s' TumorArray{array} '.M*A_s' TumorArray{array}]);
        set(reaction,'Notes','MHC translocation');
    
    % Rename MHC Molecules
    rename(M_e,['M' num2str(ID)]);
    rename(M_s,['M' num2str(ID)]);

end