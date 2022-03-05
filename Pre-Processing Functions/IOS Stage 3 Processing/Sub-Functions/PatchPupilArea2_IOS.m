function [] = PatchPupilArea2_IOS(procDataFileID)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: The pupil camera occasionally drops packets of frames. We can calculate the difference in the number
%            of expected frames as well as the indeces that LabVIEW found the packets lost.
%________________________________________________________________________________________________________________________

load(procDataFileID)
if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
    if isfield(ProcData.data.Pupil,'pupilPatch2') == false %#ok<NODEF>
        [animalID,~,fileID] = GetFileInfo_IOS(procDataFileID);
        % pupil area and derivative
        pupilArea = ProcData.data.Pupil.pupilArea;
        diffArea = abs(diff(pupilArea));
        % threshold for interpolation
        threshold = 250;
        diffIndex = diffArea > threshold;
        % derivative figure
        summaryFigure = figure;
        subplot(2,1,1)
        p1 = plot(pupilArea);
        hold on
        p2 = plot(diffArea);
        p3 = yline(threshold);
        axis tight
        legend([p1,p2,p3],'original','abs(1st deriv)','manual threshold')
        % derivative index figure
        %         figure
        %         p1 = plot(pupilArea);
        %         hold on
        %         p2 = plot(diffIndex*max(pupilArea));
        %         axis tight
        %         legend([p1,p2],'original','diff index')
        % link adjacent indeces
        [linkedDiffIndex] = LinkBinaryEvents_IOS(gt(diffIndex,0),[30,0]);
        %         figure
        %         p1 = plot(pupilArea);
        %         hold on
        %         p2 = plot(linkedDiffIndex*max(pupilArea));
        %         axis tight
        %         legend([p1,p2],'original','linked diff index')
        % obtain edge pairs for interpolation
        edgeFound = false;
        xx = 1;
        for aa = 1:length(linkedDiffIndex)
            if edgeFound == false
                if linkedDiffIndex(1,aa) == 1
                    startEdge(xx,1) = aa;
                    edgeFound = true;
                end
            elseif edgeFound == true
                if linkedDiffIndex(1,aa) == 0
                    endEdge(xx,1) = aa;
                    xx = xx + 1;
                    edgeFound = false;
                end
            end
        end
        % edge figure
        %         figure
        %         p1 = plot(pupilArea,'k');
        %         hold on
        %         for aa = 1:length(startEdge)
        %             x1 = xline(startEdge(aa,1),'r');
        %         end
        %         for aa = 1:length(endEdge)
        %             x2 = xline(endEdge(aa,1),'b');
        %         end
        %         axis tight
        %         legend([p1,x1,x2],'original','start edge','end edge')
        % set time between edges as NaN
        nanPupilArea = pupilArea;
        for aa = 1:length(startEdge)
            startTime = startEdge(aa,1);
            endTime = endEdge(aa,1);
            nanPupilArea(startTime:endTime) = NaN;
            patchedPupilArea = fillmissing(nanPupilArea,'spline');
        end
        % comparison figure
        subplot(2,1,2)
        p1 = plot(pupilArea);
        hold on
        p2 = plot(patchedPupilArea);
        axis tight
        legend([p1,p2],'original','patched')
        % save the file to directory.
        [pathstr,~,~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Pupil Data Patching 2/'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(summaryFigure,[dirpath animalID '_' fileID '_PupilPatch2'])
        close(summaryFigure)
        % check any nans in final array
        %         disp(['Number of NaNs in final array: ' num2str(sum(isnan(patchedPupilArea)))]); disp(' ')
        % save
        ProcData.data.Pupil.patchedPupilArea = patchedPupilArea;
        ProcData.data.Pupil.pupilPatch2 = 'y';
        save(procDataFileID,'ProcData')
    end
end

end
