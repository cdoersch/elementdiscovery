% The core optimizer.  It performs gradient ascent on 
% sum_{i where Y(i)>0} max(w(:)'*X(:,i)*weights(i),0) + nrmlambda*norm(w(1:end-1))  ['the objective']
% subject to sum_{i where Y(i)<=0} max(w(:)'*X(:,i)*weights(i),0)=beta. ['the constraint']
% Note that this is not a convex optimization, and this function does not
% do anything special to correct for bad initial points.  If you give it
% a random initialization, don't expect it to produce sane results!
function [ wOut outLabels miscstate ] = doGradDescentproj( X, Y, wIn ,weights,roundid)
global ds;
disp([num2str(size(X,2)) ' points, ' num2str(sum(Y>0)) ' pos'])
try
Y = Y(:);

% read the parameters from dswork.  
dsload('ds.round.lambda');
if(dsfield(ds,'round','lambda'))
    nrmlambda=ds.round.lambda;
elseif(dsfield(ds,'conf','params','lambda'))
    nrmlambda=ds.conf.params.lambda;
else
    nrmlambda=.002;
end
if(dsfield(ds,'conf','params','stepsize'))
    stepsize=ds.conf.params.stepsize;
else
    stepsize=.1;
end
if(dsfield(ds,'conf','params','optimizationComputeLimit'))
    optimizationComputeLimit=ds.conf.params.optimizationComputeLimit;
end
dsload('ds.round.beta');
if(dsfield(ds,'round','beta'))
  beta=ds.round.beta;
elseif(dsfield(ds,'conf','params','beta'))
  beta=ds.conf.params.beta;
else
  beta=1;
end

% In the current implementation, the feature vectors don't have the
% row of -1's appended
if (~all(X(end,:) == 1))
    fprintf('Appending row of -ones to feature vector\n');
    X = [X; -ones(1, size(X,2))];
end

if(sum(Y<=0)==0)
  error('no negative points!');
end

% If this is one of the initialization rounds, we only update the bias term b.
% TODO: don't hardcode the 3.
if(roundid>-1 && roundid<=3)
  wOut=wIn(:);
  if(isempty(X))
    scores=[];
    return;
  end
  scores=(wOut(:)')*X;
  negs=scores(Y<=0);
  % fminsearch on the squared amount that the constraint is violated.  However, 
  % fminsearch can miss the zero point if x gets too large and the gradient vanishes.
  % Hence give a huge penalty for x in the region with zero gradient.
  toadd=fminsearch(@(x) (sum(x+negs((x+negs)>0))-beta)^2 +1e10*(sum((x+negs)>0)==0),0,optimset('tolX',eps))
  wOut(end)=wOut(end)-toadd;
  scores=scores+toadd;
  toadd=.0001-min(max(scores(Y<=0)),max(scores(Y>0)));
  toadd=max(0,toadd);
  wOut(end)=wOut(end)-toadd;
  outLabels=scores+toadd;
  return;
end

% If X has more rows than columns (more features per patch than patches)
% use a qr decomposition to rotate the space such that the last size(X,1)-size(X,2)
% dimensions are all zero.  Then we can discard those dimensions.  This can provide
% a big speed up for only marginal loss in numeric precision.
if(size(X,1)>size(X,2))
    X=X(1:end-1,:);
    [B,~]=qr(X);
    B=B(:,1:size(X,2));
    X=(B')*X;
    X(end+1,:)=-1;
    b=wIn(end);
    wIn=((wIn(1:end-1)')*B)';
    wIn(end+1)=b;
    ranqr=1;
else
    ranqr=0;
end

% If the weights (alphas) aren't set, then set them to all ones.
% Note that even in the version of the algorithm without cross-talk, we still
% use the weights to handle redundant detections in an image.
if(~exist('weights','var'))
  weights=ones(size(Y));
  if(dsbool(ds.conf.params,'posweight'))
    weights(Y==1)=ds.conf.params.posweight;
  end
end
weights=weights(:)';

w = wIn;

w_orig=w;

% Keep track of the number of vector-matrix multiplies, approximately
% tracking the compute time.
nprods=0;
i=0;
% This should never be needed in the current implementation, but this counts
% the number of times we've 'fixed' a bad bias (see below).
nfails=0;
% The norm of the update to w.
nrmstep=Inf;
% The number of times we've taken a step without crossing any hinges (i.e.
% a point where w'*X(:,i) changes sign).
nquadmins=0;
while(nprods<optimizationComputeLimit && nrmstep > 1e-8 && nquadmins<10)
    i=i+1;
    wX = w'*X;
    nprods=nprods+1;

    % If we received a bad initialization, it might be the case that none
    % of the negative points have a score above zero, and so the gradient
    % of the constraint will be zero and everything explodes.  Hence we 'fix'
    % the bias by increasing it until at least something is above zero.
    while(sum(weights((wX)'>0 & Y<=0).*(wX((wX)'>0 & Y<=0)))==0)
      disp('too few members, fixing bias...');
      nfails=nfails+1;
      if(nfails>100000)
          save([dsload('.ds.conf.dispoutpath') '/biasfail' num2str(floor(rand*1000))]);
          error('too many bias fails');
      end
      toadd=.01;
      wX=wX+toadd;
      w(end)=w(end)-toadd;
    end
    
    % On the first round only, we use fminsearch to set the bias directly
    % in a way that satisfies the constraint.  
    if(i==1)
      negs=wX(Y<=0);
      negweights=weights(Y<=0);
      toadd=fminsearch(@(x) ((x+negs((x+negs)>0))*negweights((x+negs)>0)'-beta)^2 +1e10*(sum((x+negs)>0)==0),0,optimset('tolX',eps))
      w(end)=w(end)-toadd;
      wX=wX+toadd;
    end

    % compute the gradient of the objective sum_{i where Y(i)>0} max(w(:)'*X(:,i)*weights(i),0)
    % which is just the weighted sum of the positive datapoints whose score is above zero.
    lfull = wX.*weights;
    lfull(wX<0)=0;
    ldashfull = (wX>=0).*weights;
    onSet  = lfull' > 0;
    contribIdx = Y>0;
    lminus = sum(lfull((~contribIdx) & onSet));
    ldashplus  = X(:,contribIdx& onSet) * ldashfull(contribIdx & onSet)'; 

    % compute the gradient of norm(w)
    ldashnrm=2*w*nrmlambda;
    ldashnrm(end)=0;

    % project the gradient onto the tangent space of the constraint at the point w.
    projspace=X(:,(~contribIdx) & onSet)*weights((~contribIdx) & onSet)';
    nprods=nprods+1;
    projspace=projspace/norm(projspace);
    nabla=ldashplus-ldashnrm;
    nabla = nabla-projspace*(nabla'*projspace);

    % Take the affine approximation of the current constraint set--i.e. assume we can
    % take as large a step as we want in the curent direction nabla without violating
    % the constraint.  It's a simple quadratic equation. Compute the best step we can take.
    normquad=-[sum(nabla(1:end-1).^2),sum(2.*w(1:end-1).*nabla(1:end-1)),sum(w(1:end-1).^2)]*nrmlambda;
    ldashplusquad=ldashplus'*nabla;
    quadmin=(normquad(2)+ldashplusquad)/(-2*normquad(1));
    % numerical inaccuracy may cause quadmin to be slightly negative, in which case we
    % should simply correct it to zero.  If it's too negative (<-1e-8), something is broken
    % in the optimizer.
    if(quadmin<-1e-8)
      error('quadmin was negative!');
    end
    if(quadmin<0)
      quadmin=0;
    end
    if(isnan(quadmin))
      error('quadmin was nan');
    end
    % TODO: get rid of this section
    % next, compute how big a step we can take without causing the sign of the score any patch
    % to change.
    % dall contains the gradient of the score of each patch with respect to the step size
    %dall=nabla'*X;
    % now dall contains the step which causes a sign change for a particular point
    %dall=-wX./dall;
    % any negative points are actually getting farther from the hinge as we take a larger
    % step; ignore them.
    %dall(dall<0)=Inf;
    % now get the smallest step.
    %[dall,dallidx]=sort(dall,'ascend');

    % The alternative is a fixed step size.  In practice the optimization is 
    % slightly faster if the step size is smaller when the weights are 
    % larger, but altering the fixed step size using this heuristic
    % probably destroys all convergence guarantees.  Note the optimization 
    % still works without this heuristic.
    step=stepsize/sum(weights(contribIdx&onSet));
    wold=w;

    if(step<quadmin)
      nquadmins=0;
      % take the step.
      w=w+nabla*step;
      nrmstep=norm(nabla*step);
      negWx=w'*X(:,(~contribIdx));
      nprods=nprods+1;

      % In the (mathematically impossible) event that we took a step which made all 
      % negative datapoints less than zero, we need to fix the bias again.
      while(sum(weights(~contribIdx).*negWx)==0)
        disp('too few members after step, fixing bias...');
        nfails=nfails+1;
        if(nfails>100000)
            save([dsload('.ds.conf.dispoutpath') '/biasfail' num2str(floor(rand*1000))]);
            error('too many bias fails');
        end
        toadd=.01;
        negwX=negwX+toadd;
        w(end)=w(end)-toadd;
      end

      % Now that we have taken a step with respect to the objective function,
      % the constraint is probably violated.  Hence, we perform gradient descent
      % on f(w)=|sum_{i where Y(i)<=0} max(w(:)'*X(:,i)*weights(i),0)-beta|
      negweight=weights(~contribIdx);
      totneg=negWx(negWx>0)*negweight(negWx>0)';
      neggt0=~contribIdx;

      % We need to do careful bookkeeping on the set of patches whose score
      % is greater than zero--i.e. those patches that contribute to the gradient
      % of f.  The gradient is constant except when we hit a hinge, at which
      % point a patch will either begin or stop contributing to the gradient
      % of f.  Hence, whenever we take a step, we step all the way to one
      % of these hinges.  However, once we've taken the step, we can't just
      % test whether w'*X(:,i)>0 for all i, because one patch will
      % score close to 0, but not exactly zero due to precision issues.
      % We hence need to *force* that particular patch to either start or stop
      % contributing to the gradient of f.  neggt0 does that bookkeeping.
      %
      % I am actually not entirely sure that this does the bookkeeping correctly--
      % i.e. I think that if two patches have their scores hit 0 at the same time,
      % one of them might not have their corresponding entry in neggt0 flipped
      % properly.  In practice, however, I have not noticed this becoming a problem,
      % and any errors will get corrected on the next iteration anyway when neggt0
      % gets recomputed from scratch.
      neggt0(~contribIdx)=negWx>0;
      nprojs=0;
      while(1)
        % negWx contains the score of each negative point, and projdir contains
        % the gradient of f.  Note, however, that both of these quantities
        % can be computed incrementally, since we know the gradient of the
        % score of each point at every step of the optimization, and we know that
        % this gradient stays constant, and we know the contribution of each point
        % to the gradient.  Error accumulates, however, so we
        % re-compute negWx every 20 iterations.  The constant 20 has not been tuned.
        if(mod(nprojs,20)==0)
          projdir=X(:,neggt0)*weights(neggt0)';
          nprods=nprods+1;
          if(nprojs>0)
            negWx=w'*X(:,(~contribIdx));
            nprods=nprods+1;
          end
        end
        % compute the gradient of the score of each of the positives with respect to the
        % current projection direction.
        % Note that this could be computed incrementally if we pre-computed the kernel
        % matrix X'*X, and this would probably lead to a big speedup.  Haven't tried it
        % though.
        dnegs=projdir'*X(:,(~contribIdx));
        nprods=nprods+1;
        % compute how far we have to step before the constraint is satisfied, assuming
        % we don't hit any hinges.
        zer=(beta-totneg)./(dnegs(neggt0(~contribIdx))*weights(neggt0)');
        % If it turns out we're stepping in the wrong direction, step in the right direction.
        % Necessary because of the absolute value.
        if(zer<0)
          projdir=-projdir;
          zer=-zer;
          dnegs=-dnegs;
        end

        % Figure out how far we can step before each of the datapoints hits their hinge.
        signchg=-negWx./dnegs;
        % If the sign change is negative, it means the current gradient direction is going
        % away from that point's hinge.
        signchg(signchg<=1e-8)=Inf;
        % If we hit a hinge before we satisfy the constraint, step to that hinge.  Otherwise,
        % step to the point which satisfies the constraint.
        if(min(signchg)<zer && zer>1e-10)
          nprojs=nprojs+1;
          [msc,idx]=min(signchg);
          fContribIdx=find(~contribIdx);
          % wXidx is the index (i.e. column of X) of the point that hit its hinge
          wXidx=fContribIdx(idx);
          w=w+projdir*msc;
          negWx=negWx+dnegs*msc;
          newwx=zeros(size(wX));
          newwx(~contribIdx)=negWx;
          if(neggt0(wXidx))
            projdir=projdir-X(:,wXidx)*weights(wXidx);
          else
            projdir=projdir+X(:,wXidx)*weights(wXidx);
          end
          neggt0(wXidx)=~neggt0(wXidx);
          totneg=newwx(neggt0)*weights(neggt0)';
        else
          w=w+projdir*zer;
          break;
        end
      end
    else
      w=w+nabla*quadmin;
      nrmstep=norm(nabla*quadmin);
      disp('quadminned');
    end
    nrmw=(norm(w(1:end-1)));
    
    % In general, we're running this optimization on a tiny subset of the total
    % patch space.  Hence, if w changes too much, we're computing ratios in a 
    % part of patch space where we have few or no samples, and so our estimate
    % of the gradient is going to be bad.  This is an ugly hack to terminate
    % if it seems like we've drifted too far, which we ran out of space
    % in the paper to explain.  Let p be a point on the intersection of the decision
    % boundary defined by w_orig and the unit hypersphere.  Let theta be the angle
    % between w_orig(1:end-1) and p.  Let Q be the set of all points q on the unit 
    % hypershpere such that the angle between w_orig(1:end-1) and q is theta/2.
    % Make sure that Q lies entirely on the positive side of the hyperplane
    % defined by w.  Perform the same test with the roles of w and w_orig reversed.
    % If both tests pass, we keep going; otherwise we return for more patches.
    anglechange=getangle(w(end),w(1:end-1))*.5>getangle(w_orig(end),w_orig(1:end-1))-acos(dot(w_orig(1:end-1),w(1:end-1))/(norm(w(1:end-1))*norm(w_orig(1:end-1))));
    anglechange=anglechange || getangle(w_orig(end),w_orig(1:end-1))*.5>getangle(w(end),w(1:end-1))-acos(dot(w_orig(1:end-1),w(1:end-1))/(norm(w(1:end-1))*norm(w_orig(1:end-1))));
    angorig=getangle(w_orig(end),w_orig(1:end-1));
    angnew=getangle(w(end),w(1:end-1));
    angdiffval=acos(dot(w_orig(1:end-1),w(1:end-1))/(norm(w(1:end-1))*norm(w_orig(1:end-1))));
    if(anglechange&&~dsbool(ds.conf.params,'nocheck'))
      w=wold;
      if(anglechange) 
        disp('optimization wants to leave safe zone.  returning to get more patches.');
        break;
      end
    end
    % print out a bunch of debug info.
    if(mod(i,25)==0||nprods>=optimizationComputeLimit)
      disp(['done: ' num2str(i) ' ' num2str(toc) ' seconds']);
      disp(['lminus:' num2str(lminus)]);
      disp(['w_end:' num2str(w(end))]);
      disp(['ang_orig:' num2str(angorig)]);
      disp(['ang_new:' num2str(angnew)]);
      disp(['ang_dot:' num2str(angdiffval)]);
      disp(['nabla_end:' num2str(nabla(end))]);
      disp(['nrmstep:' num2str(nrmstep)]);
      disp(['step:' num2str(step)]);
      disp(['quadmin:' num2str(quadmin)]);
      disp(['norm w:' num2str(nrmw)]);
      disp(['active wX (pos neg): (' num2str(sum(wX(:)>0 & Y(:)>0)) ' ' num2str(sum(wX(:)>0 & Y(:)<=0)) ')']);
    end
end
disp(['made it through ' num2str(i) ' rounds']);
wX = w'*X;
outLabels=wX;
if(sum(outLabels>0)<2)
  error('too few patches in cluster...');
end

wOut = w;

if(ranqr);
    b=wOut(end);
    wOut(end)=[];
    wOut=B*wOut;
    wOut(end+1)=b;
end
wOut(end)=(wOut(end));
catch ex;dsprinterr;end
end

function res=getangle(bias,vec)
  res=acos(-bias/norm(vec));
end
