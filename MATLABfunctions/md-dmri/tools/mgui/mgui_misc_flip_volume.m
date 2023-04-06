function [I,p] = mgui_misc_flip_volume(I, ori_init, ori_final)
% function I = mgui_misc_flip_volume(I, ori_init, ori_final)

if (numel(ori_init) ~= 3) || (numel(ori_final) ~= 3)
    error('Ori must have 3 elements');
end

if (strcmp(ori_init, ori_final))
    p = [1 2 3];
    return;
end


p = zeros(1,3);
for c = 1:3
    switch (ori_init(c))
        case {'R', 'L'}
            p(find( (ori_final == 'R') | (ori_final == 'L'), 1)) = c;
        case {'A', 'P'}
            p(find( (ori_final == 'A') | (ori_final == 'P'), 1)) = c;
        case {'S', 'I'}
            p(find( (ori_final == 'S') | (ori_final == 'I'), 1)) = c;
    end
    
end

if (any(p == 0))
    error('Strange orientation string');
end


% Intermediate step, flippling left
I = permute(I, [p 4]);
ori_init = ori_init(p);


for c = 1:3
    if (ori_init(c) ~= ori_final(c))
        I = flipdim(I, c);
        switch (ori_init(c))
            case 'L'
                ori_init(c) = 'R';
            case 'R'
                ori_init(c) = 'L';
            case 'P'
                ori_init(c) = 'A';
            case 'A'
                ori_init(c) = 'P';
            case 'I'
                ori_init(c) = 'S';
            case 'S'
                ori_init(c) = 'I';
        end                
    end    
end


if (~strcmp(ori_init, ori_final))
    error('Could not flip dimensions');
end

