function [CDK2RATIO_FULL, NUCAREA_FULL] = LoadTraceDataEllipTrack(path,rows,sites,numcols)
%Load in all tracedata
CDK2RATIO_FULL = [];
NUCAREA_FULL = [];
for n = 1:numcols
    filepath = path;
    cd(filepath);
    folders = dir(filepath);
    CDK2RATIO = [];
    NUCAREA = [];
    DAUGHTERS = [];
    adder = []; z=1;
    for row=rows
        for col=n
            for site=sites
                try
                fprintf(['\n' num2str(row) '_' num2str(col) '\n'])
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                load([filepath,'\', shot, '\signal.mat']);
                cells = all_signals{row,col,site};
                adder = [adder;length(cells)];
                CDK2RATIO_singles = [];
                NUCAREA_singles = [];
                DAUGHTERS_singles = [];
                    for L = 1:length(cells)
                        CDK2RATIO_singles(L,:) = ((cells{L,1}.CDK2_cytoring_mean)./(cells{L,1}.CDK2_nuc_mean))';
                        NUCAREA_singles(L,:) = ((cells{L,1}.nuc_area))';
                        if ~isempty(find(~cellfun(@isempty,cells{L,1}.daughters)))
                            DAUGHTERS_singles{L,1} = find(~cellfun(@isempty,cells{L,1}.daughters));
                            DAUGHTERS_singles{L,2} = cells{L,1}.daughters{DAUGHTERS_singles{L,1}};
                        else
                            DAUGHTERS_singles{L,1} = nan;
                            DAUGHTERS_singles{L,2} = nan;
                        end
                        if z > 1 & ~isempty(find(~cellfun(@isempty,cells{L,1}.daughters)))
                            DAUGHTERS_singles{L,1} = DAUGHTERS_singles{L,1} + 0;
                            DAUGHTERS_singles{L,2} = DAUGHTERS_singles{L,2} + adder(z-1);
                        end
                    end
                CDK2RATIO = [CDK2RATIO;CDK2RATIO_singles];
                NUCAREA = [NUCAREA;NUCAREA_singles];
                DAUGHTERS = [DAUGHTERS;DAUGHTERS_singles];
                z=z+1;
                catch
                end
            end
        end
    end
    
    %Arrays for mitotic timing and genealogy linkage
    MITOSIS = cell2mat(DAUGHTERS(:,1));
    GENEALOGY = [];
    for L = 1:length(MITOSIS)
        GENEALOGY(L,1) = DAUGHTERS{L,2}(1,1);
        try
        GENEALOGY(L,2) = DAUGHTERS{L,2}(1,2);
        catch
        end
    end
    GENEALOGY(GENEALOGY == 0) = NaN;
    
    %Start from end of lineage and find family tree IDs
    mitotic_estimate = 40; %Estimate for max # of mitosis during movie
    loops = mitotic_estimate-2;
    BackTrack = find(~isnan(CDK2RATIO(:,end)));
    BackTrack = num2cell(BackTrack);
    for L = 1:length(BackTrack)
       BackTrack{L,2} = find(GENEALOGY == BackTrack{L,1});
       z = 1;
           for Z = 1:loops
                if ~isempty(BackTrack{L,z+1}) 
                    BackTrack{L,z+2} = find(GENEALOGY == BackTrack{L,z+1});
                else
                    BackTrack{L,z+2} = [];
                end
                z=z+1;
           end
    end
    Lineage = [];
    for L = 1:length(BackTrack)
        TEMP = (BackTrack(L,:));
        TEMP = TEMP(~cellfun('isempty',TEMP));
        for N = 1:length(TEMP)
            if TEMP{1,N} > size(CDK2RATIO,1)
               TEMP{1,N} = TEMP{1,N} - size(CDK2RATIO,1);
            end
        end
        Lineage{L,1} = flip(TEMP);
    end
    
    %Concatenate traces based on lineage IDs
    CDK2Full = [];
    AreaFull = [];
    for W = 1:length(BackTrack)
        t = []; a = [];
        for L = 1:length(Lineage{W,1})
            t = [t,CDK2RATIO(Lineage{W,1}{1,L},:)];
            a = [a,NUCAREA(Lineage{W,1}{1,L},:)];
        end
        
        t = t(~isnan(t));
        a = a(~isnan(a));
            if length(t) < size(CDK2RATIO,2)
                ending = nan(1,size(CDK2RATIO,2) - length(t));
                t = [t,ending];
            elseif length(t) > size(CDK2RATIO,2)
                t = nan(1:size(CDK2RATIO,2));
            end
            
            if length(a) < size(CDK2RATIO,2)
                ending = nan(1,size(CDK2RATIO,2) - length(a));
                a = [a,ending];
            elseif length(a) > size(CDK2RATIO,2)
                a = nan(1:size(CDK2RATIO,2));
            end  
            
        CDK2Full = [CDK2Full;t];
        AreaFull = [AreaFull;a];
    end
    
    %Array of all mitoses that occur along the lineage for each cell
    MitosisFull = [];
    for W = 1:length(BackTrack)
        l = [];
        for L = 1:length(Lineage{W,1})-1
            out = cell2mat(DAUGHTERS(Lineage{W,1}{1,L},1));
            l = [l;out];
        end
        MitosisFull{W,1} = l;
    end
    
    %Trim to full traces only
    AreaFull(any(isnan(CDK2Full), 2), :) = [];
    MitosisFull(any(isnan(CDK2Full), 2), :) = [];
    CDK2Full(any(isnan(CDK2Full), 2), :) = [];
    
    CDK2RATIO_FULL = [CDK2RATIO_FULL;CDK2Full];
    NUCAREA_FULL = [NUCAREA_FULL;AreaFull];
end
