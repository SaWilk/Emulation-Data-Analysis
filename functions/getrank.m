function tmprank2 = getrank(tmpdata)
    
    tmprank = rank(tmpdata);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Here: alternate computation of the rank by Sven Hoffmann
    %tmprank = rank(tmpdata(:,1:min(3000, size(tmpdata,2)))); old code
    covarianceMatrix = cov(tmpdata', 1);
    [E, D] = eig (covarianceMatrix);
    rankTolerance = 1e-7;
    tmprank2=sum (diag (D) > rankTolerance);
    if tmprank ~= tmprank2
        fprintf('Warning: fixing rank computation inconsistency (%d vs %d) most likely because running under Linux 64-bit Matlab\n', tmprank, tmprank2);
        tmprank2 = min(tmprank, tmprank2);
    end

