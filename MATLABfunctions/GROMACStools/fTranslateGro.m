function gro_out = fTranslateGro(gro_in,TranslationVector)

if ~isfield(TranslationVector,'x')
    TranslationVector.x = 0;
end
if ~isfield(TranslationVector,'y')
    TranslationVector.y = 0;
end
if ~isfield(TranslationVector,'z')
    TranslationVector.z = 0;
end

gro = gro_in;

gro.x = gro.x + TranslationVector.x;
gro.y = gro.y + TranslationVector.y;
gro.z = gro.z + TranslationVector.z;

gro_out = gro;