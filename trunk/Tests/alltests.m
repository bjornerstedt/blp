% Run all tests
clear all

ff = runtests
if all([ff(:).Passed])
    load handel
else
    load gong
end
sound(y,Fs)

% return
% 
% import matlab.unittest.TestSuite
% import matlab.unittest.TestRunner
% import matlab.unittest.plugins.CodeCoveragePlugin
% 
% suite = TestSuite.fromPackage('tests');
% runner = TestRunner.withTextOutput;
% 
% runner.addPlugin(CodeCoveragePlugin.forFolder(pwd))
% result = runner.run(suite);