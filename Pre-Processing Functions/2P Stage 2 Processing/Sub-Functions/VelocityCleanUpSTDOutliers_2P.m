
function [fixed_v]=VelocityCleanUpSTDOutliers_2P(v, thresh)
%finds points outside a threshold velocity and interpolated between last
% and next good points, unless at the begining or end, then puts in mean
% velocities
%pd 2-21-08
upper_bound=thresh*std(v)+mean(v);
lower_bound=-thresh*std(v)+mean(v);
badpoints=union(find(v>upper_bound),find(v<lower_bound));
noutliers=length(badpoints);
fixed_v=v;
while noutliers>0
    %badpoints=find(abs(v)>thresh);
    goodpoints=union(setdiff(badpoints+1,badpoints),setdiff(badpoints-1,badpoints));
    for p=1:length(badpoints)
        lastgood=max(find(goodpoints<badpoints(p)));
        nextgood=min(find(goodpoints>badpoints(p)));
        if ((goodpoints(1)>0)&&(goodpoints(end)<=length(fixed_v)))
            if((numel(lastgood)>0)&&(isempty(nextgood)==0))
                fixed_v(badpoints(p))=.5*fixed_v(goodpoints(lastgood))+.5*fixed_v(goodpoints(nextgood));
            elseif((numel(lastgood)>0))
                fixed_v(badpoints(p))=fixed_v(goodpoints(nextgood));
            else
                fixed_v(badpoints(p))=fixed_v(goodpoints(lastgood));
            end
        elseif (goodpoints(1)==0)
            fixed_v(badpoints(p))=mean(fixed_v);%(goodpoints(nextgood));
        elseif (goodpoints(end)==length(fixed_v)+1)
            fixed_v(badpoints(p))=mean(fixed_v);%fixed_v(goodpoints(lastgood));
        elseif numel(goodpoints)==0
            fixed_v(badpoints(p))=-fixed_v;
        end      
    end
    upper_bound=thresh*std(fixed_v)+mean(fixed_v);
    lower_bound=-thresh*std(fixed_v)+mean(fixed_v);
    badpoints=union(find(fixed_v>upper_bound),find(fixed_v<lower_bound));
    noutliers=length(badpoints);   
end
end