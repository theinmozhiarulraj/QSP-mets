% function to sample the fraction of different antigens in cancer clones
% Input - ncancer_clones -- Number of cancer clones
%         nAgs           -- Number of neoAgs
% Output - sampled fractions of different antigens in each cancer clone

function ag_fraction_clones = sample_ag_proportions(nclones_cancer,nAgs)

ag_fraction_clones=[];

for j=1:nclones_cancer

    ag_fraction=zeros(1,nAgs);
    range_max=1;
    
    for i=1:nAgs
        % r = a + (b-a).*rand(N,1).
        rno=range_max*rand(1,1);

        % if nAgs=1, fraction of antigen = 1
        if(i==nAgs)
            rno=range_max;
        end
        range_max=range_max-rno;    
        ag_fraction(i)=rno;
    
    end

    % randomly shuffle the indices
    indices = randperm(length(ag_fraction));  
    ag_fraction_shuffled = ag_fraction(indices);
    
    %ag_fraction_shuffled
    ag_fraction_clones=[ag_fraction_clones;ag_fraction_shuffled];

end

