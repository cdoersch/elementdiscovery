function Z = linkage(Y, method, pdistArg, varargin)
%LINKAGE Create hierarchical cluster tree.
%   Z = LINKAGE(X), where X is a matrix with two or more rows, creates a
%   matrix Z defining a tree of hierarchical clusters of the rows of X.
%   Clusters are based on the single linkage algorithm using Euclidean
%   distances between the rows of X. Rows of X correspond to observations
%   and columns to variables.
%
%   Z = LINKAGE(X,METHOD) creates a hierarchical cluster tree using the
%   the specified algorithm. The available methods are:
%
%      'single'    --- nearest distance (default)
%      'complete'  --- furthest distance
%      'average'   --- unweighted average distance (UPGMA) (also known as
%                      group average)
%      'weighted'  --- weighted average distance (WPGMA)
%      'centroid'  --- unweighted center of mass distance (UPGMC)
%      'median'    --- weighted center of mass distance (WPGMC)
%      'ward'      --- inner squared distance (min variance algorithm)
%
%   Z = LINKAGE(X,METHOD,METRIC) performs clustering based on the distance
%   metric METRIC between the rows of X. METRIC can be any of the distance
%   measures accepted by the PDIST function. The default is 'euclidean'.
%   For more information on PDIST and available distances, type HELP PDIST.
%   The 'centroid', 'median', and 'ward' methods are intended only for the
%   Euclidean distance metric.
%
%   Z = LINKAGE(X,METHOD,PDIST_INPUTS) enables you to pass extra input
%   arguments to PDIST. PDIST_INPUTS should be a cell array containing
%   input arguments to be passed to PDIST.
%
%   Z = LINKAGE(X,METHOD,METRIC,'savememory',VALUE) specifies the value for
%   the optional argument 'savememory'. VALUE is a string 'on' or 'off'.
%   When VALUE is 'on', LINKAGE performs clustering without computing the
%   distance matrix internally, so it usually requires less memory for data
%   with large number of observations. VALUE 'on' requires both of the
%   following two conditions to be satisfied:
%      * METHOD is 'ward', 'centroid', or 'median'.
%      * METRIC is 'euclidean'.
%   When the above two conditions are satisfied, the default is 'on' when 
%   the number of columns in X is not greater than 20, or the machine
%   doesn't have enough memory to save the distance matrix.
%
%   Z = LINKAGE(Y) and Z = LINKAGE(Y,METHOD) are alternative syntaxes that
%   accept a vector representation Y of a distance matrix. Y may be a
%   distance matrix as computed by PDIST, or a more general dissimilarity
%   matrix conforming to the output format of PDIST.
%
%   The output matrix Z contains cluster information. Z has size m-1 by 3,
%   where m is the number of observations in the original data. Column 1
%   and 2 of Z contain cluster indices linked in pairs to form a binary
%   tree. The leaf nodes are numbered from 1 to m. They are the singleton
%   clusters from which all higher clusters are built. Each newly-formed
%   cluster, corresponding to Z(i,:), is assigned the index m+i, where m is
%   the total number of initial leaves. Z(i,1:2) contains the indices of
%   the two component clusters which form cluster m+i. There are m-1 higher
%   clusters which correspond to the interior nodes of the output
%   clustering tree. Z(i,3) contains the corresponding linkage distances
%   between the two clusters which are merged in Z(i,:), e.g. if there are
%   total of 30 initial nodes, and at step 12, cluster 5 and cluster 7 are
%   combined and their distance at this time is 1.5, then row 12 of Z will
%   be (5,7,1.5). The newly formed cluster will have an index 12+30=42. If
%   cluster 42 shows up in a latter row, that means this newly formed
%   cluster is being combined again into some bigger cluster.
%
%   The 'centroid' and 'median' methods can produce a cluster tree that is
%   not monotonic. This occurs when the distance from the union of two
%   clusters, r and s, to a third cluster is less than the distance between
%   r and s. In such a case, in a dendrogram drawn with the default
%   orientation, the path from a leaf to the root node takes some downward
%   steps. You may want to use another method when that happens.
%
%   You can provide the output Z to other functions including DENDROGRAM to
%   display the tree, CLUSTER to assign points to clusters, INCONSISTENT to
%   compute inconsistent measures, and COPHENET to compute the cophenetic
%   correlation coefficient.
%
%   Examples:
%      % Compute four clusters of the Fisher iris data using Ward linkage
%      % and ignoring species information, and see how the cluster
%      % assignments correspond to the three species.
%      load fisheriris
%      Z = linkage(meas,'ward','euclidean');
%      c = cluster(Z,'maxclust',4);
%      crosstab(c,species)
%      dendrogram(Z)
%
%     % Create a hierarchical cluster tree for a data with 20000
%     % observations using Ward linkage. If you set the value of
%     % 'savememory' to 'off', you may get an out-of-memory error if your
%     % machine doesn't have enough memory to hold the distance matrix.
%     X = rand(20000,3);
%     Z = linkage(X,'ward','euclidean','savememory', 'on')
%
%   See also PDIST, INCONSISTENT, COPHENET, DENDROGRAM, CLUSTER,
%   CLUSTERDATA, KMEANS, SILHOUETTE.

%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.10 $

% Check for size and type of input
[k, n] = size(Y);
m = ceil(sqrt(2*n)); % m = (1+sqrt(1+8*n))/2, but works for large n
if k>1  % data matrix input
    if nargin<2
        method = 'single';
    end
    if nargin<3
        pdistArg = 'euclidean';
    end
    nargs = 3;
else % distance matrix input or bad input
    nargs = nargin;
end

if nargs>=3 % should be data input
    if k == 1 && m*(m-1)/2 == n
         warning(message('stats:linkage:CallingPDIST'));
    end
    if k < 2
        error(message('stats:linkage:TooFewDistances'));
    end
    passX = true; % the input is the data matrix.
    
    if ischar(pdistArg)
        pdistArg = {pdistArg};
    elseif ~iscell(pdistArg)
        error(message('stats:linkage:BadPdistArgs'));
    end
        
else % should be distance input
    passX = false;
    if n < 1
        error(message('stats:linkage:TooFewDistances'));
    end
    if k ~= 1 || m*(m-1)/2 ~= n
        error(message('stats:linkage:BadSize'));
    end
    
end

% Selects appropriate method
if nargs == 1 % set default switch to be 'si'
    methodStr = 'single';
else
    %           Preferred   Synonym
    methods = {'single',   'nearest'; ...
               'complete', 'farthest'; ...
               'average',  'upgma'; ...
               'weighted', 'wpgma'; ...
               'centroid', 'upgmc'; ...
               'median',   'wpgmc'; ...
               'ward''s',  'incremental'};
    [methodStr,s] = internal.stats.getParamVal(method,methods(:),'METHOD');

    if s>size(methods,1)
        methodStr = methods{s-size(methods,1), 1};
    end
end
method = methodStr(1:2);

% The recursive distance updates for the three methods 'cetroid','median'
% and 'ward' only make sense when the distance matrix contains Euclidean
% distances (which will be squared) or the distance metric is Euclidean
nonEuc = false;
isCeMeWa = false;
if any(strcmp(method,{'ce' 'me' 'wa'}))
    isCeMeWa = true;
    if ~passX % The distance matrix is passed
        if (any(~isfinite(Y)) || ~iseuclidean(Y))
            warning(message('stats:linkage:NotEuclideanMatrix', methodStr));
        end
    else % the data matrix is passed
      
        if (~isempty(pdistArg))
            if (~ischar (pdistArg{1}))
                nonEuc = true;
            else
                %  method may be a known name or the name of a user function
                distMethods = {'euclidean'; 'minkowski';'mahalanobis'};
                i = find(strncmpi(pdistArg{1}, distMethods, length(pdistArg{1})));
                if length(i) > 1
                    error(message('stats:linkage:BadDistance', pdistArg{ 1 }));
                elseif (isempty(i) || i == 3 || ...
                  (i == 2 && length(pdistArg) ~= 1 && isscalar(pdistArg{2}) && pdistArg{2} ~= 2) )
                    nonEuc = true;
                end
            end
            
        end
        if (nonEuc)
            warning(message('stats:linkage:NotEuclideanMethod', methodStr));
        end
    end
end

%Parse the memory efficient option 'Savememory'
pnames ={'savememory'};
dflts = {[]};

[memEff] = internal.stats.parseArgs(pnames, dflts, varargin{:});

if ~isempty(memEff)
    memEff = internal.stats.getParamVal(memEff,{'on'; 'off'},'SaveMemory');
    memEff = strcmp(memEff,'on');
end

% The memeory efficient option can be valid only when the following three
% conditions are satisfied:
% * The first input argument is data matrix
% * The linkage is either 'centroid, 'median' or 'ward'.
% * The distance metric is euclidean 
if passX && isCeMeWa && ~nonEuc
    if isempty(memEff) %Set the default choice for memEff
        if  n<=20
            memEff = true;
        else
            memEff = false;
        end
        
        try
            tempDistMat=zeros(1,k*(k-1)/2,class(Y));
        catch ME
            if (strcmp(ME.identifier,'MATLAB:nomem') || ...
                strcmp(ME.identifier,'MATLAB:pmaxsize'));
             memEff = true;          
            else % Otherwise, just let the error propagate.
                throw(ME);
            end
        end
        clear tempDistMat;
    end
else
    if ~isempty(memEff) && memEff
       warning(message('stats:linkage:IgnoringMemEff'));
    end
    memEff = false; 
end

if exist('linkagemex','file')==3
    % call mex file
    if passX
        if (memEff)
            %call the memory efficient algorithm
            Z = linkagemex(Y',method,{'euc'},memEff); %note that Y is transposed here
        else
           Z = linkagemex(Y,method,pdistArg, memEff);
        end
    else
        Z = linkagemex(Y,method);
    end
    
    if any(strcmp(method,{'av' 'wa' 'co' 'we'}))
    % if 'ave','ward', 'com' or 'weighted average' is used, we need to 
    % re-arrange the rows in Z matrix by sorting the third rows of Z matrix.
      Z = rearrange(Z);
    end
else
    warning(message('stats:linkage:NoMexFilePresent'));
    if passX
        Y = pdist(Y,pdistArg{:});
    end
    % optional old linkage function (use if mex file is not present)
    Z = linkageold(Y,method);
end

% Check if the tree is monotonic and warn if not.  Z is built so that the rows
% are in non-decreasing height order, thus we can look at the heights in order
% to determine monotonicity, rather than having to explicitly compare each parent
% with its children.
zdiff = diff(Z(:,3));
if any(zdiff<0)
    % With distances that are computed recursively (average, weighted, median,
    % centroid, ward's), errors can accumulate.  Two nodes that are really
    % at the same height in the tree may have had their heights calculated in
    % different ways, making them differ by +/- small amounts.  Make sure that
    % doesn't produce false non-monotonicity warnings.
    negLocs = find(zdiff<0);
    if any(abs(zdiff(negLocs)) > eps(Z(negLocs,3))) % eps(the larger of the two values)
        warning(message('stats:linkage:NonMonotonicTree', methodStr));
    end
end

% %re-arrange the Z matrix by sorting the third rows of Z matrix.
% function Z2 = rearrange(Z)
%   %Get indices to sort Z
%     [~,idx] = sort(Z(:,3));
%     
%     % Get indices in the reverse direction
%     revidx(idx) = 1:length(idx);
%     
%     % Get vector of desired cluster numbers for Z2
%     nrows=size(Z,1);
%     v2 = [1:nrows+1,nrows+1+revidx];
%     
%     % Put Z2 into sorted order, without renumbering the clusters
%     Z0 = Z(idx,:);
%     Z2 = Z0;
%     % Renumber the clusters
%     Z2(:,1:2) = v2(Z2(:,1:2));
%     
%     % Make sure the lower-numbered cluster is in column 1
%     t = Z2(:,1)>Z2(:,2);
%     Z2(t,1:2) = Z2(t,[2 1]);


%%%%%%%%%%%%%%%%%%%%%%%%% OLD LINKAGE FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Z = linkageold(Y, method)
%LINKAGEOLD Create hierarchical cluster tree using only MATLAB code.

n = size(Y,2);
m = ceil(sqrt(2*n)); % (1+sqrt(1+8*n))/2, but works for large n
if isa(Y,'single')
   Z = zeros(m-1,3,'single'); % allocate the output matrix.
else
   Z = zeros(m-1,3); % allocate the output matrix.
end

% during updating clusters, cluster index is constantly changing, R is
% a index vector mapping the original index to the current (row, column)
% index in Y.  N denotes how many points are contained in each cluster.
N = zeros(1,2*m-1);
N(1:m) = 1;
n = m; % since m is changing, we need to save m in n.
R = 1:n;

% Square the distances so updates are easier.  The cluster heights will be
% square-rooted back to the original scale after everything is done.
if any(strcmp(method,{'ce' 'me' 'wa'}))
   Y = Y .* Y;
end

for s = 1:(n-1)
   if strcmp(method,'av')
      p = (m-1):-1:2;
      I = zeros(m*(m-1)/2,1);
      I(cumsum([1 p])) = 1;
      I = cumsum(I);
      J = ones(m*(m-1)/2,1);
      J(cumsum(p)+1) = 2-p;
      J(1)=2;
      J = cumsum(J);
      W = N(R(I)).*N(R(J));
      [v, k] = min(Y./W);
   else
      [v, k] = min(Y);
   end

   i = floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)));
   j = k - (i-1)*(m-i/2)+i;

   Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A

   % Update Y. In order to vectorize the computation, we need to compute
   % all the indices corresponding to cluster i and j in Y, denoted by I
   % and J.
   I1 = 1:(i-1); I2 = (i+1):(j-1); I3 = (j+1):m; % these are temp variables
   U = [I1 I2 I3];
   I = [I1.*(m-(I1+1)/2)-m+i i*(m-(i+1)/2)-m+I2 i*(m-(i+1)/2)-m+I3];
   J = [I1.*(m-(I1+1)/2)-m+j I2.*(m-(I2+1)/2)-m+j j*(m-(j+1)/2)-m+I3];

   switch method
   case 'si' % single linkage
      Y(I) = min(Y(I),Y(J));
   case 'co' % complete linkage
      Y(I) = max(Y(I),Y(J));
   case 'av' % average linkage
      Y(I) = Y(I) + Y(J);
   case 'we' % weighted average linkage
      Y(I) = (Y(I) + Y(J))/2;
   case 'ce' % centroid linkage
      K = N(R(i))+N(R(j));
      Y(I) = (N(R(i)).*Y(I)+N(R(j)).*Y(J)-(N(R(i)).*N(R(j))*v)./K)./K;
   case 'me' % median linkage
      Y(I) = (Y(I) + Y(J))/2 - v /4;
   case 'wa' % Ward's linkage
      Y(I) = ((N(R(U))+N(R(i))).*Y(I) + (N(R(U))+N(R(j))).*Y(J) - ...
	  N(R(U))*v)./(N(R(i))+N(R(j))+N(R(U)));
   end
   J = [J i*(m-(i+1)/2)-m+j];
   Y(J) = []; % no need for the cluster information about j.

   % update m, N, R
   m = m-1;
   N(n+s) = N(R(i)) + N(R(j));
   R(i) = n+s;
   R(j:(n-1))=R((j+1):n);
end

if any(strcmp(method,{'ce' 'me' 'wa'}))
   Z(:,3) = sqrt(Z(:,3));
end

Z(:,[1 2])=sort(Z(:,[1 2]),2);



function Z = rearrange(Z)
  %Get indices to sort Z
    [~,idx] = sort(Z(:,3));
    
    % Get indices in the reverse direction
    revidx(idx) = 1:length(idx);
    
    % Get vector of desired cluster numbers for Z2
    nrows=size(Z,1);
    v2 = [1:nrows+1,nrows+1+revidx];
    
    % Put Z2 into sorted order, without renumbering the clusters
    Z = Z(idx,:);
    
    % Renumber the clusters
    Z(:,1:2) = v2(Z(:,1:2));
    
    % Make sure the lower-numbered cluster is in column 1
    t = Z(:,1)>Z(:,2);
    Z(t,1:2) = Z(t,[2 1]);
    
    % is Z in chronological order?
    if any((Z(:,2) >= (1:nrows)'+nrows+1))
        Z = fixnonchronologicalZ(Z);
        t = Z(:,1)>Z(:,2);
        Z(t,1:2) = Z(t,[2 1]);
    end

function Z = fixnonchronologicalZ(Z)
% Fixes a binary tree that has branches defined in a non-chronological order  
 
nl = size(Z,1)+1; % number of leaves
nb = size(Z,1); % number of branches
last = nl; % last defined node, we start only with leaves
for i = 1:nb
    tn = nl+i; %this node
    if any(Z(i,[1,2])>last)
        % this node (tn) uses nodes not defined yet, find a node (h) that
        % does use nodes already defined so we can interchange them:
        h = find(all(Z(i+1:end,[1 2])<=last,2),1)+i;
        % change nodes:
        Z([i h],:) = Z([h i],:);
        nn = nl+h; % new node
        % change references to such nodes
        %Z([find(Z(1:2*nb) == nn,1) find(Z(1:2*nb) == tn,1)]) = [tn nn];          
        %try
        tonn=find(Z(1:2*nb) == tn,1);
        Z([find(Z(1:2*nb) == nn,1)]) = [tn];          
        %catch,keyboard;end
        if(~isempty(tonn))
          Z(tonn) = [nn];          
        end
    end
    last = tn;
end



