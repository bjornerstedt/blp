
% Compile the demo as a mex file
if ismac
    mex -I/opt/local/include  -L/opt/local/lib -larmadillo findDeltaCpp.cpp 
    mex -I/opt/local/include  -L/opt/local/lib -larmadillo shareCalcCpp.cpp
else
    error('armaMex on PC not implemented')
end

% clear
% load testvals
% edel = findDelta2(s, expmu, iweight, edelta)
% assert(max(abs(newedel -  edel))<1e-10)
% 
% % Run the demo using X and Y
% Z = findDeltaCpp(s, expmu, iweight, edelta)

