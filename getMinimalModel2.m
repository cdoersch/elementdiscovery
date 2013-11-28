function minModel = getMinimalModel(model,suppVec)
% precomputes the weight vectors.
if(~exist('suppVec','var'))
  suppVec = model.SVs;
else
  suppVec=suppVec(model.SVs,:);
end
coeff = model.sv_coef;
coeff = repmat(coeff, 1, size(suppVec, 2));
minModel.rho = model.rho;
size(coeff)
size(suppVec)
minModel.w = sum(coeff .* suppVec);
firstLabel = model.Label(1);
minModel.w=minModel.w*firstLabel;
minModel.rho=minModel.rho*firstLabel;
if isfield(model, 'threshold')
  minModel.threshold = model.threshold;
end
if isfield(model, 'info')
  minModel.info = model.info;
end
end
