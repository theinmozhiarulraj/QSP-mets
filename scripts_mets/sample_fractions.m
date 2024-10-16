function fractions_shuffled = sample_fractions(nAgs)

    ag_fraction=zeros(1,nAgs);
    range_max=1;
    
    for i=1:nAgs
        % r = a + (b-a).*rand(N,1).
        rno=range_max*rand(1,1);

        if(i==nAgs)
            rno=range_max;
        end
        range_max=range_max-rno;
        
        ag_fraction(i)=rno;
    
    end
    
    % ag_fraction
    indices = randperm(length(ag_fraction));  
    fractions_shuffled = ag_fraction(indices);
    
end
