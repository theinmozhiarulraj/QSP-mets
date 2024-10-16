function pos = find_index(simDataPSA,species_name,compartment_name,patient_id)

    indices=find(ismember(simDataPSA(patient_id).simData.DataNames,species_name));
    
    pos=-1;
    
    for i=1:length(indices)
    
      index=indices(i);

      %disp(index)
      %disp(compartment_name)

      if(simDataPSA(patient_id).simData.DataInfo{index,1}.Compartment==compartment_name)
          pos=index;
      end
    
    end


end

