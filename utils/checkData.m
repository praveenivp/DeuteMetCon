function checkData(twix_2H,twix_B0)

% check receive coils
sig=squeeze(sos(twix_2H.image{''},[1 3:16]));

if(any(sig<0.15*mean(sig)))
fprintf('SNR low Sire! on these %d channels(<15%)\n', find(sig<0.15*mean(sig)));
end
if(any(sig<0.05*mean(sig)))
warning('SNR too low Sire! on these %d channels(<5%)\n', find(sig<0.05*mean(sig)));
end


if(nargin>1)
    % check shim currents
    if(any(getShimCurrents(twix_2H)~=getShimCurrents(twix_B0)))
        warning('Shim currents changed\n');
        fprintf('Shim currents twix_2H : %d \n',getShimCurrents(twix_2H));
        fprintf('Shim currents twix_B0 : %d \n',getShimCurrents(twix_B0));
    else
        fprintf('Shim looks alright\n')
    end

    % check adjust frequency
    if(strcmp(twix_2H.hdr.Phoenix.sTXSPEC.asNucleusInfo{1}.tNucleus,...
            twix_B0.hdr.Phoenix.sTXSPEC.asNucleusInfo{1}.tNucleus))

        freq_equal=twix_2H.hdr.Phoenix.sTXSPEC.asNucleusInfo{1}.lFrequency==...
            twix_B0.hdr.Phoenix.sTXSPEC.asNucleusInfo{1}.lFrequency;

        if(~freq_equal)
            warning('adj frequency changed\n');
        end
    end

end

end

function sim_currents=getShimCurrents(twix)
shim_1st=[twix.hdr.Phoenix.sGRADSPEC.asGPAData{1}.lOffsetX,...
    twix.hdr.Phoenix.sGRADSPEC.asGPAData{1}.lOffsetY,...
    twix.hdr.Phoenix.sGRADSPEC.asGPAData{1}.lOffsetZ];

sim_currents=[shim_1st,...
    cell2mat(twix.hdr.Phoenix.sGRADSPEC.alShimCurrent)];

end