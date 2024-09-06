% :Usage:
% ::
%
%     [y, pI, S, I] = hpfilter(y, TR, HP, spersess, varargin)
%
% :Inputs:
%
%   **y:**
%        data to be filtered
%
%   **TR:**
%        TR in seconds
%
%   **HP:**
%        Highpass filter cutoff in seconds
%
%   **spersess:**
%        scans per sessions; vector with length = # of sessions
%
% :Optional inputs: in fixed order -
%
%   **pI:**
%        pseudoinverse of I matrix, for fast application to new y 
%
%   **dummy:**
%        vector of fixed time points to model with dummy regressors in each run
%        e.g., 1:2 models first two time points in each run with
%        separate dummy regressor (appended to I)
%
%   **dooutliers:**
%        flag (1/0).  If 1, imputes session mean to time points
%        with data > 4 median absolute deviations from the session median (after
%        filtering). Not part of scnlab standard preprocessing. Use with caution.
%
% :Outputs:
%
%   **y:**
%        filtered data, session intercepts removed
%
%   **pI:**
%        intercept model X*pinv, such that y - pI * y removes intercept
%
%   **I:**
%        intercept and dummy scan design matrix
%
%   **S:**
%        smoothing model, such that S * y does HP filtering
%           if sess lengths are unequal, make as long as the longest
%           session; may do funny things to shorter ones?
%
% :Examples:
% ::
%
%    % slowest, creates intercept and smoothing matrix
%    [y, I, S] = hpfilter(data, 2, 120, [169 169 169 173 173]);
%
%    % for subsequent voxels with the same session info,
%    y = hpfilter(data, [], S, [169 169 169 173 173], I);
%
%    y = clusters(1).indiv_timeseries(:, 1);
%    [y, I, S] = hpfilter(y, 2, 120, [169 169 169 173 173]);
%    y = clusters(1).indiv_timeseries(:, 1);
%    [y] = hpfilter(y, [], S, [169 169 169 173 173], I);
%
%    % Regress out average activity in first 2 scans of each session (artifacts)
%    y = hpfilter(raw, 2, 100, EXPT.FIR.nruns, [], 1:2);
%
%    % Another example set up first, then run on multi-voxel matrix:
%    [y, pI, S, I] = hpfilter(data{1}(:,1), HPDESIGN.TR, HPDESIGN.HP, HPDESIGN.spersess, [], 1);
%    tic, y = hpfilter(data{1}, [], S, HPDESIGN.spersess, pI); toc
%
%    % But if you just have one matrix, no need to set up, so this is just as fast:
%    tic , [y, pI, S, I] = hpfilter(data{1}, HPDESIGN.TR, HPDESIGN.HP,HPDESIGN.spersess, [], 1);, toc
%
% ..
%    tor wager
%    Modified 5/12/06 to include dummy images for each session
%    Modified April 07 to add sep. dummy covs for each session and add outlier option
%    Also: y can now be a matrix; hpfilter operates on columns (faster than looping)
% ..

function [y, pI, HP, incpt] = hpfilter(y, TR, HP, spersess, varargin)
    incpt = [];
    pI = [];
    dummyscans = [];
    dooutliers = 0;
    
    spersess = spersess(:)'; % enforce as row vector

    if ~isempty(varargin), pI = varargin{1}; end
    if length(varargin) > 1, dummyscans = varargin{2}; end
    if length(varargin) > 2, dooutliers = varargin{3}; end

    if isempty(pI) || ~ismatrix(pI)
        % if we haven't entered this, then compute it from spersess
        % spersess = [number of images in each run] vector

        incpt = intercept_model(spersess, dummyscans);
        pI = incpt * pinv(incpt);
    end

    n = size(y, 1);
    if size(pI, 1) > n
        disp('Warning: Intercept has more observations than data. spersess is wrong?')
        if size(incpt, 1) > n
            incpt = incpt(1:n, :);
        end
        
        pI = pI(1:n, 1:n);

        tmp = cumsum(spersess);
        
        if tmp(end) > length(y)
            disp('Warning: hpfilter: Images in spersess does not match length of data to be filtered.');
        end
        
        while tmp(end) > length(y)
            disp('Trying removal of a session.');
            spersess = spersess(1:end-1);
            tmp = cumsum(spersess);
        end
    end

    % remove intercept
    y = y - pI * y;

    if ~ismatrix(HP) || size(HP, 1) ~= size(y, 1)
        % it's not already an SPM smoothing matrix
        len = max(spersess); % max number of images in a session (run)
        HP = use_spm_filter(TR, len, 'none', 'specify', HP);
    end

    % starting and ending images for each session
    st = cumsum([1 spersess(1:end-1)]);
    en = cumsum(spersess);

    % high-pass filter here
    for i = 1:length(st)
        y(st(i):en(i), :) = HP(1:spersess(i), 1:spersess(i)) * y(st(i):en(i), :);
    end
    
    % because intercept removal is applied before HP filtering, some residual
    % effects of intercepts may remain, so remove these:
    y = y - pI * y;

if dooutliers

    n = size(y,1);

    % identify outliers and replace with timeseries mean
    % (mean has no influence on model betas)
    mady = repmat(mad(y), n, 1);
    wh = (abs(y) > 4 * mady);

    subjectmeans = repmat(mean(y, 1), n, 1);
    y(wh) = subjectmeans(wh);
end

    
end % HPFILTER


% ISMATRIX: Returns 1 if the input matrix is 2+ dimensional, 0 if it is a scalar 
%           or vector.
%
%     Usage ismat = ismatrix(X)
%
% RE Strauss, 5/19/00

function ismat = ismatrix(X)
  [r,c] = size(X);
  if (r>1 && c>1)
    ismat = 1;
  else
    ismat = 0;
  end

end % ISMATRIX





function [S,KL,KH] = use_spm_filter(TR,dims,LChoice,HChoice,HParam,varargin)
% :Usage:
% ::
%
%     function [S,KL,KH] = use_spm_filter(TR,dim of filter,LChoice,HChoice,HP filter in s,[LP Gauss len in s])
%
% :Inputs:
%
%   **K{s}.LChoice:**
%        Low-pass  filtering {'hrf' 'Gaussian' 'none'}
%   **K{s}.LParam:**
%        Gaussian parameter in seconds
%   **K{s}.HChoice:**
%        High-pass filtering {'specify' 'none'}
% ..
%    05/22/01 Tor Wager
% ..

K{1}.RT = TR; 
K{1}.LChoice = LChoice;
K{1}.HChoice = HChoice;
K{1}.HParam = HParam; 
K{1}.row = ones(dims, 1);

if length(varargin) > 0
        K{1}.LParam = varargin{1};
end

    
KL = []; KH = [];

spmS = spm_filter('set', K);

S = eye(length(K{1}.row));

if ~strcmp(HChoice,'none')
    
    KH = full(spmS{1}.KH);
    S = S - KH * pinv(KH);
    
end

if ~strcmp(LChoice,'none')
    KL = full(spmS{1}.KL);
    S = KL * S;             % lowpass * highpass; hp = I - res forming mtx of S.KH
end

return
end % USE_SPM_FILTER




function [vargout] = spm_filter(Action,K,Y)
% filter routine
% FORMAT [K] = spm_filter('set',K)
% FORMAT [Y] = spm_filter('apply',K,Y)
%
% Action    - 'set'   fills in filter structure K
% Action    - 'apply' applies K to Y = K*Y
% K         - filter convolution matrix or:
% K{s}      - cell of structs containing session-specific specifications
%
% K{s}.RT       - repeat time in seconds
% K{s}.row      - row of Y constituting session s
% K{s}.LChoice  - Low-pass  filtering {'hrf' 'Gaussian' 'none'}
% K{s}.LParam   - Gaussian parameter in seconds
% K{s}.HChoice  - High-pass filtering {'specify' 'none'}
% K{s}.HParam   - cut-off period in seconds
%
% K{s}.HP       - low frequencies to be removed
% K{s}.LP       - sparse toepltz low-pass convolution matrix
% 
% Y         - data matrix
%
% K         - filter structure
% Y         - filtered data K.K*Y
%___________________________________________________________________________
%
% spm_filter implements band pass filtering in an efficient way by
% using explicitly the projector matrix form of the High pass
% component.  spm_filter also configures the filter structure in
% accord with the specification fields if required
%___________________________________________________________________________
% @(#)spm_filter.m  2.4 Karl Friston 99/08/31


% set or apply
%---------------------------------------------------------------------------
switch Action

    case 'set'
    %-------------------------------------------------------------------
    for s = 1:length(K)

        % matrix order
        %-----------------------------------------------------------
        k     = length(K{s}.row);

        % make low pass filter
        %-----------------------------------------------------------
        switch K{s}.LChoice

            case 'none'
            %---------------------------------------------------
            h       = 1;
            d       = 0;

            case 'hrf'
            %---------------------------------------------------
            h       = spm_hrf(K{s}.RT);
            h       = [h; zeros(size(h))];
            g       = abs(fft(h));
            h       = real(ifft(g));
            h       = fftshift(h)';
            n       = length(h);
            d       = [1:n] - n/2 - 1;

            case 'Gaussian'
            %---------------------------------------------------
            sigma   = K{s}.LParam/K{s}.RT;
            h       = round(4*sigma);
            h       = exp(-[-h:h].^2/(2*sigma^2));
            n       = length(h);
            d       = [1:n] - (n + 1)/2;
            if      n == 1, h = 1; end

            otherwise
            %---------------------------------------------------
            error('Low pass Filter option unknown');
            return

        end

        % create and normalize low pass filter
        %-----------------------------------------------------------
        K{s}.KL = spdiags(ones(k,1)*h,d,k,k);
        K{s}.KL = spdiags(1./sum(K{s}.KL')',0,k,k)*K{s}.KL;


        % make high pass filter
        %-----------------------------------------------------------
        switch K{s}.HChoice

            case 'none'
            %---------------------------------------------------
            K{s}.KH = [];

            case 'specify'
            %---------------------------------------------------
            n       = fix(2*(k*K{s}.RT)/K{s}.HParam + 1);
            X       = spm_dctmtx(k,n);
            K{s}.KH = sparse(X(:,[2:n]));

            otherwise
            %---------------------------------------------------
            error('High pass Filter option unknown');
            return

        end

    end

    % return structure
    %-------------------------------------------------------------------
    vargout = K;


    case 'apply'
    %-------------------------------------------------------------------
    if iscell(K)


        % ensure requisite feild are present
        %-----------------------------------------------------------
        if ~isfield(K{1},'KL')
            K = spm_filter('set',K);
        end

        for s = 1:length(K)

            % select data
            %---------------------------------------------------
            y = Y(K{s}.row,:);

            % apply low pass filter
            %---------------------------------------------------
            y = K{s}.KL*y;

            % apply high pass filter
            %---------------------------------------------------
            if ~isempty(K{s}.KH)
                y = y - K{s}.KH*(K{s}.KH'*y);
            end

            % reset filtered data in Y
            %---------------------------------------------------
            Y(K{s}.row,:) = y;

        end

    % K is simply a convolution matrix
    %-------------------------------------------------------------------
    else
        Y = K*Y;
    end

    % return filtered data
    %-------------------------------------------------------------------
    vargout   = Y;


    otherwise
    %-------------------------------------------------------------------
    warning('Filter option unknown');


end

end % SPM_FILTER



function C = spm_dctmtx(N,K,n,f)
% Creates basis functions for Discrete Cosine Transform.
% FORMAT C = spm_dctmtx(N,K,n)
%     OR C = spm_dctmtx(N,K)
%     OR D = spm_dctmtx(N,K,n,'diff')
%     OR D = spm_dctmtx(N,K,'diff')
% N - dimension
% K - order
% n - optional points to sample
%____________________________________________________________________________
% spm_dctmtx creates a matrix for the first few basis functions of a one
% dimensional discrete cosine transform.
% With the 'diff' argument, spm_dctmtx produces the derivatives of the
% DCT.
%
% See:    Fundamentals of Digital Image Processing (p 150-154).
%         Anil K. Jain 1989.
%____________________________________________________________________________
% @(#)spm_dctmtx.m  1.3 John Ashburner MRCCU/FIL 96/08/14

d = 0;

if (nargin == 2)
    n = (0:(N-1))';
    if (nargin == 3)
        d = 1;
    end
elseif (nargin == 3)
    if (strcmp(n,'diff'))
        d = 1;
        n = (0:(N-1))';
    else
        n = n(:);
    end
elseif (nargin == 4)
    n = n(:);
    if (strcmp(f,'diff'))
        d = 1;
    else
        error('Incorrect Usage');
    end
else
    error('Incorrect Usage');
end

C = zeros(size(n,1),K);

if (d == 0)
    C(:,1)=ones(size(n,1),1)/sqrt(N);
    for k=2:K
        C(:,k) = sqrt(2/N)*cos(pi*(2*n+1)*(k-1)/(2*N));
    end
else
    for k=2:K
        C(:,k) = -2^(1/2)*(1/N)^(1/2)*sin(1/2*pi*(2*n*k-2*n+k-1)/N)*pi*(k-1)/N;
    end
end
end % SPM_DCTMTX


function x = intercept_model(nvols_per_run, varargin)
% Build design matrix X for intercepts
% given vector of session lengths [s1 s2 s3] in images
%
% :Usage:
% ::
%
%     x = intercept_model(nvols_per_run, [indx of dummy scans in each session])
%
% :Examples:
% ::
%
%    nvols_per_run = [166 166 144 137];
%    x = intercept_model(nvols_per_run);
%
%    x = intercept_model(repmat(166, 1, 5));
%
%    Xi = intercept_model(EXPT.FIR.nruns, 1:2);
%
% ..
%    tor modified april 07: separate column for each run
% ..

    if ~isvector(nvols_per_run), error('nvols_per_run must be a vector of the number of volumes for each run'); end
    nvols_per_run = nvols_per_run(:)'; %ensure row vector
    
    npoints = sum(nvols_per_run); % number of time points total
    nruns = length(nvols_per_run); % number of runs

    x = zeros(npoints, nruns);

    st = cumsum([1 nvols_per_run]);
    en = st(2:end) - 1; % ending values
    st = st(1:end-1); % starting values

    for i = 1:nruns
        x(st(i):en(i), i) = 1;
    end

    if length(varargin) > 0 && ~isempty(varargin{1})
        % Model dummy regressors
        x = [x model_dummy(npoints, st', varargin{1})];
    end
end % INTERCEPT_MODEL



function x = model_dummy(npoints, st, wh)
    % dummy regressors for first scan or two (at least, that's the intended
    % use)
    % separate column for each run!

    k = length(wh) .* length(st);
    x = zeros(npoints, k);

    dosamecol = 0;

    if dosamecol
        for i = 1:k
            ind = st(r) + wh(i) - 1; % which to filter
            ind(ind < 1) = []; % if negative numbers (end of runs), this makes it ok
            x(ind, i) = 1;
        end

    else

        colindx = 1;
        
        for r = 1:length(st) % for each run
            for i = 1:length(wh)
        
                ind = st(r) + wh(i) - 1; % which to filter
                ind(ind < 1) = []; % if negative numbers (end of runs), this makes it ok
                
                x(ind, colindx) = 1;
                colindx = colindx + 1;
            end
        end

    end

end % MODEL_DUMMY

