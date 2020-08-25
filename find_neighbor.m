function z = find_neighbor(x)
% x is the new sample; z is the sample point with the extremum
% seed that is the closest to the new extremum seed.
global ObSet

num = 1;
    for j=2:length(ObSet)-1
        if(ObSet(j).level==1)
            distt(num)=norm(x.seed-ObSet(j).seed,2);  % compute the disttance of the new extremum seed to the other extremum in the first level
            sample_ind(num)=j;
            num = num+1;
        end
    end 
   
    if(exist('distt','var'))
        ind = find(distt==min(distt));                    % find the index of the extremum seed with the minimum disttance to the new extremum seed
        z = ObSet(sample_ind(ind(1)));                  % obtain the corresponding sample point
    else
        z = x;
    end
end