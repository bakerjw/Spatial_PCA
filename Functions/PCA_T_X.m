%% function PCA_T_to_X = PCA_T_X()
% Function that takes from transformed space back to the original space
% Input:
% coef          Coefficient Matrix (txn)
% mu            Mean of the original space (1xt)
% transformed   Matrix in the transformed space (mxn)
% Output:
% original      Matrix in the original space (mxt)

function original = PCA_T_X(transformed, coef, mu)
    original = bsxfun(@plus,transformed*coef',mu); 
end