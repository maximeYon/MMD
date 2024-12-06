function atoms = fReadAtoms(fname)

format = '%d%s%d%s%s%d%f%f';

fid = fopen(fname, 'r');
text.atoms = textscan(fid, format, 'headerlines', 1,'CommentStyle',';');
fclose(fid);

atoms.nr = text.atoms{1};
atoms.type = text.atoms{2};
atoms.resnr = text.atoms{3};
atoms.residue = text.atoms{4};
atoms.atom = text.atoms{5};
atoms.cgnr = text.atoms{6};
atoms.charge = text.atoms{7};
atoms.mass = text.atoms{8};

