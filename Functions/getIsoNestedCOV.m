function [COV] = getIsoNestedCOV(varModel,DIST)

var = varModel.Cn+varModel.C1+varModel.C2;
COV = var - (varModel.Cn+varModel.C1*(1-exp(-3*DIST/varModel.a1))+varModel.C2*(1-exp(-3*DIST/varModel.a2)));
COV(DIST == 0) = var;

end