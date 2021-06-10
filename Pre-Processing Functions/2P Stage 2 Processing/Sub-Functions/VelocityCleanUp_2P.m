
function [fixed_v]=VelocityCleanUp_2P(v, thresh)
%finds points faster than a threshold velocity and interpolated between last
% and next good points
%pd 11-30-07

badpoints=find(abs(v)>thresh);
goodpoints=union(setdiff(badpoints+1,badpoints),setdiff(badpoints-1,badpoints));
fixed_v=v;
for p=1:length(badpoints)
    lastgood=max(find(goodpoints<badpoints(p)));
    nextgood=min(find(goodpoints>badpoints(p)));
   if ((goodpoints(1)>0)&&(goodpoints(end)<=length(v)))        
        if((numel(lastgood)>0)&&(isempty(nextgood)==0))
            fixed_v(badpoints(p))=.5*v(goodpoints(lastgood))+.5*v(goodpoints(nextgood));
        elseif((numel(lastgood)>0))
            fixed_v(badpoints(p))=v(goodpoints(nextgood));
        else
            fixed_v(badpoints(p))=v(goodpoints(lastgood));
        end
    elseif (goodpoints(1)==0)
        fixed_v(badpoints(p))=v(goodpoints(nextgood));
    elseif (goodpoints(end)==length(v)+1)
        fixed_v(badpoints(p))=v(goodpoints(lastgood));
    elseif numel(goodpoints)==0
        fixed_v(badpoints(p))=-v;
    end
    
end
end