function [threshSize] = clustthresh1D(cond1, cond2, nperm, tail, cfalpha, csalpha)


%returns cluster threshold based on cluster size for group-level paired t-test betweeen cond 1 and cond 2 (must
%be 2D and same size with subjects in second dim), or between cond1 and 0
%if second input is zero.

switch nargin
    case 1
        error ('needs at least 2 inputs')
    case 2
        nperm = 5000;
        tail = 'both';
        cfalpha = 0.05;
        csalpha = 0.05;
    case 3 
        tail = 'both';
        cfalpha = 0.05;
        csalpha = 0.05;
    case 4
        cfalpha = 0.05;
        csalpha = 0.05;
    case 5
        csalpha = 0.05;
end
                

    if isscalar(cond2)
        
        for prm = 1:nperm
            flipind = logical(round(rand(1,size(cond1,2))));
            permcond = cond1;
            permcond(:,flipind) = permcond(:,flipind).*-1;
            h = ttest(permcond,cond2,'dim',2,'tail',tail,'alpha',cfalpha);
            cs(prm) = max(clustsize_thisfunc(h));           
        end
        
        threshSize = prctile(cs,100*(1-csalpha));
    
    
    else
        allcond = cat(2, cond1, cond2);
        
        for prm = 1:nperm 
            permind = randsample(size(allcond,2),size(allcond,2));
            [h,p,ci,stat] = ttest(allcond(:,permind(1:end/2)), allcond(:,end/2+1:end),'dim',2,'tail',tail,'alpha',cfalpha);
            cs(prm) = max(clustsize_thisfunc(h)); 
        end
        
        threshSize = prctile(cs,100*(1-csalpha));
            
        
    end
    
    fprintf(['nperm ' num2str(nperm) '\ntail ' tail '\ncluster-forming alpha ' num2str(cfalpha), '\ncluster-size alpha ' num2str(csalpha) '\n'])

    
    function [c, mapl]=clustsize_thisfunc(h)
        
        %counts number of pixels in each cluster of ones in the binary matrix h
        %returns numbered clusters of ones in mapl

        [mapl,nblobs] = bwlabeln(h);
        if nblobs >0
            clustcount = zeros(1,nblobs);
                for i=1:nblobs
                    clustcount(i) = sum(mapl(:)==i);
                end
            c = clustcount;
        else 
            c = 0;

        end
    end

end
