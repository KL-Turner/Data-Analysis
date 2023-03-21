function [comparison] = FindSolenoidComparison(hemisphere,solenoid)

if strcmp(solenoid,'AudSol') == true
    comparison = 'aud';
elseif strcmp(solenoid,'LPadSol')
    if any(strcmp(hemisphere,{'LH','frontalLH'}))
        comparison = 'ipsi';
    else
        comparison = 'contra';
    end
elseif strcmp(solenoid,'RPadSol')
    if any(strcmp(hemisphere,{'RH','frontalRH'}))
        comparison = 'ipsi';
    else
        comparison = 'contra';
    end
end

end