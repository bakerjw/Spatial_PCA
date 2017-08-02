function [COV] = getNugCOV(varModel,DIST)

COV = eye(size(DIST))*varModel.Cn;

end