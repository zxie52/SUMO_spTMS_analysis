function [c, mapl]=clustsize(h)

%counts number of pixels in the each cluster of ones in the binary matrix h
%returns numbered clusters of ones in mapl

    if nansum(h(:)) == 0
        c=0;
    else
    mapl = bwlabel(h);
    nblobs = max(mapl(:));
    clustcount = zeros(1,nblobs);
        for i=1:nblobs
            clustcount(i) = sum(mapl(:)==i);
        end
    c = clustcount;
    end
