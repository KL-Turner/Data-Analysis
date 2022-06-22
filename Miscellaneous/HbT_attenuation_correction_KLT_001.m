%% READ ME
% 
% DO NOT USE THIS SCRIPT IT IS GARBAGE- KWG 10/04/2021
%
%% Hemodynamic response correction code using Mle
    if strcmpi(correctionFlag,'y') 
        binEdges = (0:0.001:1);
        figure;
        dataHist = histogram2(rescaledData(:,3),rescaledData(:,2),binEdges,binEdges);
        binCounts = dataHist.BinCounts;
        allCounts = sum(binCounts,'all');
        allHist = (binCounts/allCounts)*100;
        xCounts = sum(binCounts,2);
        mleFit.RH.allHist = allHist;
        for binNum = 2:size(binCounts,1)
            if binNum == 2
                percData(1) = (xCounts(1)/sum(xCounts))*100;
            end
            percData(binNum) = (sum(xCounts(1:binNum))/sum(xCounts))*100;
        end
        startInd = find(percData <= 2,1,'last') + 1; % exclude lower 2% of data points
        endInd = find(percData >= 98,1,'first') - 1; % exclude upper 2% of data points
        mleFit.RH.startInd = startInd;
        mleFit.RH.endInd = endInd;
        for binNum = startInd:endInd
            index = (binNum - startInd) + 1;
            leftEdge = binEdges(binNum);
            rightEdge = binEdges(binNum + 1);
            dataInds = rescaledData(:,3) >= leftEdge & rescaledData(:,3) <= rightEdge;
            theData = rescaledData(dataInds,1);
            [phat,pci] = mle(theData,'distrbution','normal');
            figure(101);
            normHist = histogram(theData,binEdges);
            normCounts = normHist.BinCounts./sum(normHist.BinCounts);
            theFit = pdf('normal',normHist.BinEdges(1:(end - 1)),phat(1),phat(2));
            theFit = theFit./sum(theFit);
            rsqr = 1 - (sum((normCounts - theFit).^2)/sum((normCounts - mean(normCounts)).^2));
            mleFit.RH.fitnotes.phat(index,:) = phat;
            mleFit.RH.fitnotes.pci(:,:,index) = pci;
            mleFit.RH.goodness.rsqr(index) = rsqr;
            mleFit.RH.fitData.eGFPData{index} = theData;
            mleFit.RH.fitData.BinEdges(index,:) = [leftEdge,rightEdge];
            mleFit.RH.fitData.normHist(index,:) = normCounts;
            mleFit.RH.fitData.HistFit(index,:) = theFit;
            mleFit.RH.fitData.fitHistEdges(index,:) = normHist.BinEdges;
        end
        [RH_linFit] = fitlm(binEdges(startInd:endInd)',mleFit.RH.fitnotes.phat(:,1),'RobustOpts','on');
        mleFit.RH.linFit = RH_linFit;
        % Plot the fit
        RH_yInt = table2array(mleFit.RH.linFit.Coefficients(1,1));
        RH_slope = table2array(mleFit.RH.linFit.Coefficients(2,1));
        RH_finalGCaMP = rescaledData(:,2) - (RH_slope*rescaledData(:,3) + RH_yInt);
        Corrected465.RH = RH_finalGCaMP;
        % LH
        binEdges = 0:0.001:1;
        figure;
        dataHist = histogram2(rescaledData(:,6),rescaledData(:,5),binEdges,binEdges);
        binCounts = dataHist.BinCounts;
        allCounts = sum(binCounts,'all');
        allHist = (binCounts/allCounts)*100;
        xCounts = sum(binCounts,2);
        mleFit.LH.allHist = allHist;
        for binNum = 2:size(binCounts,1)
            if binNum == 2
                percData(1) = (xCounts(1)/sum(xCounts))*100;
            end
            percData(binNum) = (sum(xCounts(1:binNum))/sum(xCounts))*100;
        end
        startInd = find(percData <= 2,1,'last') + 1; % exclude lower 2% of data points
        endInd = find(percData >= 98,1,'first') - 1; % exclude upper 2% of data points
        mleFit.LH.startInd = startInd;
        mleFit.LH.endInd = endInd;
        for binNum = startInd:endInd
            index = (binNum - startInd) + 1;
            leftEdge = binEdges(binNum);
            rightEdge = binEdges(binNum + 1);
            dataInds = rescaledData(:,6) >= leftEdge & rescaledData(:,6) <= rightEdge;
            theData = rescaledData(dataInds,1);
            [phat,pci] = mle(theData,'distrbution','normal');
            figure(102);
            normHist = histogram(theData,binEdges);
            normCounts = normHist.BinCounts./sum(normHist.BinCounts);
            theFit = pdf('normal',normHist.BinEdges(1:(end - 1)),phat(1),phat(2));
            theFit = theFit./sum(theFit);
            rsqr = 1 - (sum((normCounts - theFit).^2)/sum((normCounts - mean(normCounts)).^2));
            mleFit.LH.fitnotes.phat(index,:) = phat;
            mleFit.LH.fitnotes.pci(:,:,index) = pci;
            mleFit.LH.goodness.rsqr(index) = rsqr;
            mleFit.LH.fitData.eGFPData{index} = theData;
            mleFit.LH.fitData.BinEdges(index,:) = [leftEdge,rightEdge];
            mleFit.LH.fitData.normHist(index,:) = normCounts;
            mleFit.LH.fitData.HistFit(index,:) = theFit;
            mleFit.LH.fitData.fitHistEdges(index,:) = normHist.BinEdges;
        end
        [LH_linFit] = fitlm(binEdges(startInd:endInd)',mleFit.LH.fitnotes.phat(:,1),'RobustOpts','on');
        mleFit.LH.linFit = LH_linFit;
        % Plot the fit
        LH_yInt = table2array(mleFit.LH.linFit.Coefficients(1,1));
        LH_slope = table2array(mleFit.LH.linFit.Coefficients(2,1));
        LH_finalGCaMP = rescaledData(:,5) - (LH_slope*rescaledData(:,6) + LH_yInt);
        Corrected465.LH = LH_finalGCaMP;
    else
        Corrected465.RH = rescaledData(:,2);
        Corrected465.LH = rescaledData(:,5);
    end