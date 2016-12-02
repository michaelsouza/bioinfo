%% set parameters
protein = '1gpv'; % PDB ID
modelID = 1;      % model identification

fprintf('Creating an artifical DMDGP instance\n');

%% download pdb data
fprintf('Downloading protein data from PDB ....');
% contatenate file name and pdb http address
file = sprintf('http://www.rcsb.org/pdb/files/%s.pdb', protein);

% try to download data
try
    pdb  = pdbread(file);
catch
    error('\nThe data for protein %s could not be found.\n', protein);
end
fprintf(' done!\n');

%% instance creation
fprintf('Instance creation\n');
% select the model
Atom = pdb.Model(modelID).Atom;

chainID = [Atom.chainID];
element = [Atom.element];
resSeq  = [Atom.resSeq];

resSeqID = unique(resSeq);

for i = 1:length(resSeqID)
    index = resSeq == resSeqID(i);
    fprintf('ResSeq(%d)\n', resSeqID(i))
    fprintf('   ChainID '); disp(chainID(index));
    fprintf('   Element '); disp(element(index));
end