function [R] = ATLSearchlight(condition, subjects, varargin)
% NOTE: If running as a batch of jobs, the set of windows should be held
% constant over the batch, unless you define sub-directories where chuncks
% of windows will be written.
% DEPENDENCIES
% 1. glmnet
% 2. sqrt_truncate_r.m (from WholeBrain_RSA/src)
    p = inputParser();
    if isdeployed
        addRequired(p, 'condition');
        addOptional(p, 'subjects', '1 2 3 5 7 8 9 10', @ischar);
        addParameter(p, 'radius', '6', @ischar);
        addParameter(p, 'parallel', '0', @ischar);
        addParameter(p, 'perms', '1 2 3 4 5 6 7 8 9 10', @ischar);
        addParameter(p, 'permutationindex', 'permutation_indexes.mat', @ischar);
        addParameter(p, 'WRITE_PERM_INDEXES', '0', @ischar);
        addParameter(p, 'dataroot', '/mnt/sw01-home01/mbmhscc4/scratch/data/Naming_ECoG/avg');
        parse(p, condition, subjects, varargin{:});
        fprintf('condition: %s\n', p.Results.condition);
        fprintf('subjects: %s\n', p.Results.subjects);
        fprintf('parallel: %s\n', p.Results.parallel);
        fprintf('permutations: %s\n\n', p.Results.permutationindex);

    else
        addRequired(p, 'condition');
        addOptional(p, 'subjects', [1:3,5,7:10]);
        addParameter(p, 'radius', 6);
        addParameter(p, 'parallel', 0);
        addParameter(p, 'perms', 1:10, @isnumeric);
        addParameter(p, 'permutationindex', 'permutation_indexes.mat', @ischar);
        addParameter(p, 'dataroot', fullfile('D:','MRI','SoundPicture','data','MAT','avg','bystudy'));
        addParameter(p, 'WRITE_PERM_INDEXES', false, @ischar);
        parse(p, condition, subjects, varargin{:});
        fprintf('condition: %s\n', p.Results.condition);
        fprintf('subjects: %d\n', p.Results.subjects);
        fprintf('parallel: %d\n', p.Results.parallel);
        fprintf('permutations: %s\n\n', p.Results.permutationindex);
    end

    if isdeployed
        PARALLEL = logical(str2double(p.Results.parallel));
        subjects_set = str2double(strsplit(p.Results.subjects));
        dataroot = p.Results.dataroot;
        perms = str2double(strsplit(p.Results.perms));
        permfile = p.Results.permutationindex;
        radius = str2double(p.Results.radius);
        WRITE_PERM_INDEXES = logical(p.Results.WRITE_PERM_INDEXES);
    else
        addpath('C:\Users\mbmhscc4\MATLAB\src\WholeBrain_RSA\src\');
        addpath('C:\Users\mbmhscc4\MATLAB\Toolboxes\glmnet');
        PARALLEL = p.Results.parallel;
        subjects_set = p.Results.subjects;
        dataroot = p.Results.dataroot;
        perms = p.Results.perms;
        permfile = p.Results.permutationindex;
        radius = p.Results.radius;
        WRITE_PERM_INDEXES = logical(p.Results.WRITE_PERM_INDEXES);
    end

    load(fullfile(dataroot, 'metadata_avg.mat'));
    NSUBJ = numel(subjects_set);
    R = struct('subject',[],'modality',[],'roi',[],'X',[]);
    
    if PARALLEL && isempty(gcp('nocreate'))
        ppp = parpool('local');
    end

    nperm = numel(perms);
    if isempty(permfile)
        PERMS = zeros(37,nperm);
        for i = 1:nperm
            PERMS(:,i) = randperm(37);
        end
    else
        load(permfile, 'PERMS');
    end
    if WRITE_PERM_INDEXES
        save('permutation_indexes.mat', 'PERMS');
    end
    
    fortran_error_log = fopen('fortran_error.log', 'w');
    for i = 1:NSUBJ
        s = subjects_set(i);
        M = selectbyfield(metadata, 'subject', s);
        f = fullfile(dataroot, sprintf('s%02d_avg.mat', s));
        S = selectbyfield(M.targets, 'label', 'semantic', 'type', 'similarity', 'sim_source', 'featurenorms', 'sim_metric', 'cosine');
        C = sqrt_truncate_r(S.target, 0.3);
        COORDS = selectbyfield(M.coords, 'orientation', 'orig');
        colfilter = [selectbyfield(M.filters, 'label', 'colfilter_aud'); selectbyfield(M.filters, 'label', 'ROI_semantic')];
        rowfilter = selectbyfield(M.filters, 'label', 'rowfilter_aud');
        cf = all(cat(1,colfilter.filter));
        rf = rowfilter.filter;

        D = squareform(pdist(COORDS.xyz(cf,:)));
        SL = cell(size(D,1),1);
        for j = 1:size(D,1)
            SL{j} = find(D(j,:) <= radius);
        end

        CV = M.cvind(rf,1);

        switch condition
            case 'audio'
                %%%%%%%
                % Audio
                load(f, 'audio');
                R(i,1).subject = s;
                R(i,1).modality = 'audio';
                R(i,1).roi = 'semantic';
                R(i,1).X = audio(rf, cf);
                R(i,1).SL = SL;
                R(i,1).radius = radius;
                R(i,1).err1 = nan(9, numel(SL));
                R(i,1).err2 = nan(9, numel(SL));
                R(i,1).perm_err1 = nan(9, numel(SL), 100);
                R(i,1).perm_err2 = nan(9, numel(SL), 100);
                R(i,1).perms = perms;
                R(i,1).permutation_indexes = PERMS(:,perms);
                nvox = numel(SL);
                disp('Audio')
                for j = 1:9
                    fprintf('cv: %d\n', j);
                    nchar = 0;
                    for k = 1:nvox
                        if mod(k,10) == 1;
                            fprintf(repmat('\b',1,nchar));
                            nchar = fprintf('%0.0f%%', (k/numel(SL)) * 100);
                        end
                        test = CV == j;
                        train = ~test;
                        glmnet_cv = CV(train);
                        glmnet_cv(glmnet_cv>j) = (glmnet_cv(glmnet_cv>j)) - 1;
                        Ct = C(train,:);
                        X = R(i,1).X(:,SL{k});
                        Xt = X(train,:);
                        try
                            opts_cv = glmnetSet(struct('mtype','grouped','alpha',1));
                            fitobj_cv = cvglmnet(Xt, Ct, 'mgaussian', opts_cv, 'mse', 8, glmnet_cv, PARALLEL);
                            opts = glmnetSet(struct('mtype','grouped','lambda',fitobj_cv.lambda_min,'alpha',1));
                            fitobj = glmnet(Xt, Ct, 'mgaussian', opts);
                            Cz = glmnetPredict(fitobj,X);
                            R(i,1).err1(j,k) = norm(C(test,:)-Cz(test,:), 'fro') ./ norm(C(test,:), 'fro');
                            R(i,1).err2(j,k) = norm(C(train,:)-Cz(train,:), 'fro') ./ norm(C(train,:), 'fro');
                        catch ME
                            fprintf(fortran_error_log, '%d,%d,%d,%d\n', s,j,k,-1);
                        end
                    end
                    fprintf(repmat('\b',1,nchar));
                    fprintf('%0.0f%%\n\n', 100);
                end

                %%%%%%%%%%%%%%%%%%%%
                % AUDIO PERMUTATIONS
                disp('Audio (permutations)')
                for j = 1:9
                    fprintf('cv: %d\n', j);
                    nchar = 0;
                    for k = 1:nvox
                        if mod(k,10) == 1;
                            fprintf(repmat('\b',1,nchar));
                            nchar = fprintf('%0.0f%%', (k/numel(SL)) * 100);
                        end
                        for q = 1:nperm
                            p = perms(q);
                            test = CV == j;
                            train = ~test;
                            glmnet_cv = CV(train);
                            glmnet_cv(glmnet_cv>j) = (glmnet_cv(glmnet_cv>j)) - 1;
                            Ct = C(train,:);
                            % HERE IS THE PERMUTATION
                            pix = PERMS(rf,p);
                            [~,ix] = sort(pix);
                            [~,pix] = sort(ix); % this is a hack to compensate of the rowfilter
                            X = R(i,1).X(pix,SL{k});
                            Xt = X(train,:);                         
                            try
                                opts_cv = glmnetSet(struct('mtype','grouped','alpha',1));
                                fitobj_cv = cvglmnet(Xt, Ct, 'mgaussian', opts_cv, 'mse', 8, glmnet_cv, PARALLEL);
                                opts = glmnetSet(struct('mtype','grouped','lambda',fitobj_cv.lambda_min,'alpha',1));
                                fitobj = glmnet(Xt, Ct, 'mgaussian', opts);
                                Cz = glmnetPredict(fitobj,X);
                                R(i,1).err1(j,k) = norm(C(test,:)-Cz(test,:), 'fro') ./ norm(C(test,:), 'fro');
                                R(i,1).err2(j,k) = norm(C(train,:)-Cz(train,:), 'fro') ./ norm(C(train,:), 'fro');
                            catch ME
                                fprintf(fortran_error_log, '%d,%d,%d,%d\n',s,j,k,p);
                            end
                        end
                    end
                    fprintf(repmat('\b',1,nchar));
                    fprintf('%0.0f%%\n\n', 100);
                end
                
            case 'visual'
                %%%%%%%%
                % Visual
                load(f, 'visual');
                colfilter = [selectbyfield(M.filters, 'label', 'colfilter_vis'); selectbyfield(M.filters, 'label', 'ROI_semantic')];
                rowfilter = selectbyfield(M.filters, 'label', 'rowfilter_vis');
                cf = all(cat(1,colfilter.filter));
                rf = rowfilter.filter;
                R(i,2).subject = s;
                R(i,2).modality = 'visual';
                R(i,2).roi = 'semantic';
                R(i,2).X = visual(rf, cf);
                R(i,2).SL = SL;
                R(i,2).radius = radius;
                R(i,2).err1 = nan(9, numel(SL));
                R(i,2).err2 = nan(9, numel(SL));
                R(i,2).perm_err1 = nan(9, numel(SL), 100);
                R(i,2).perm_err2 = nan(9, numel(SL), 100);
                R(i,2).perms = perms;
                R(i,2).permutation_indexes = PERMS(:,perms);
                nvox = numel(SL);
                CV = M.cvind(rf,1);
                disp('Visual');
                for j = 1:9
                    fprintf('cv: %d\n', j);
                    nchar = 0;
                    for k = 1:nvox
                        if mod(k,10) == 1;
                            fprintf(repmat('\b',1,nchar));
                            nchar = fprintf('%0.0f%%', (k/numel(SL)) * 100);
                        end
                        test = CV == j;
                        train = ~test;
                        glmnet_cv = CV(train);
                        glmnet_cv(glmnet_cv>j) = (glmnet_cv(glmnet_cv>j)) - 1;
                        Ct = C(train,:);
                        X = R(i,2).X(:,SL{k});
                        Xt = X(train,:);
                        try
                            opts_cv = glmnetSet(struct('mtype','grouped','alpha',1));
                            fitobj_cv = cvglmnet(Xt, Ct, 'mgaussian', opts_cv, 'mse', 8, glmnet_cv, PARALLEL);
                            opts = glmnetSet(struct('mtype','grouped','lambda',fitobj_cv.lambda_min,'alpha',1));
                            fitobj = glmnet(Xt, Ct, 'mgaussian', opts);
                            Cz = glmnetPredict(fitobj,X);
                            R(i,2).err1(j,k) = norm(C(test,:)-Cz(test,:), 'fro') ./ norm(C(test,:), 'fro');
                            R(i,2).err2(j,k) = norm(C(train,:)-Cz(train,:), 'fro') ./ norm(C(train,:), 'fro');
                        catch ME
                            fprintf(fortran_error_log, '%d,%d,%d,%d\n', s,j,k,-1);
                        end
                    end
                    fprintf(repmat('\b',1,nchar));
                    fprintf('%0.0f%%\n\n', 100);
                end

                %%%%%%%%%%%%%%%%%%%%%
                % VISAUL PERMUTATIONS
                disp('Visual');
                for j = 1:9
                    fprintf('cv: %d\n', j);
                    nchar = 0;
                    for k = 1:nvox
                        if mod(k,10) == 1;
                            fprintf(repmat('\b',1,nchar));
                            nchar = fprintf('%0.0f%%', (k/numel(SL)) * 100);
                        end
                        for q = 1:nperm
                            p = perms(q);
                            test = CV == j;
                            train = ~test;
                            glmnet_cv = CV(train);
                            glmnet_cv(glmnet_cv>j) = (glmnet_cv(glmnet_cv>j)) - 1;
                            Ct = C(train,:);
                            % HERE IS THE PERMUTATION
                            pix = PERMS(rf,p);
                            [~,ix] = sort(pix);
                            [~,pix] = sort(ix); % this is a hack to compensate of the rowfilter
                            X = R(i,2).X(pix,SL{k});
                            Xt = X(train,:);
                            try
                                opts_cv = glmnetSet(struct('mtype','grouped','alpha',1));
                                fitobj_cv = cvglmnet(Xt, Ct, 'mgaussian', opts_cv, 'mse', 8, glmnet_cv, PARALLEL);
                                opts = glmnetSet(struct('mtype','grouped','lambda',fitobj_cv.lambda_min,'alpha',1));
                                fitobj = glmnet(Xt, Ct, 'mgaussian', opts);
                                Cz = glmnetPredict(fitobj,X);
                                R(i,2).err1(j,k) = norm(C(test,:)-Cz(test,:), 'fro') ./ norm(C(test,:), 'fro');
                                R(i,2).err2(j,k) = norm(C(train,:)-Cz(train,:), 'fro') ./ norm(C(train,:), 'fro');
                            catch ME
                                fprintf(fortran_error_log, '%d,%d,%d,%d\n', s,j,k,p);
                            end
                        end
                    end
                    fprintf(repmat('\b',1,nchar));
                    fprintf('%0.0f%%\n\n', 100);
                end
        end
    end
    save('results.mat', 'R');
    if PARALLEL
        delete(ppp);
    end
    fclose(fortran_error_log);
end