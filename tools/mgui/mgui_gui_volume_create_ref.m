function ref = mgui_gui_volume_create_ref(tags, fn, In, f_callback)
% function ref = mgui_gui_volume_create_ref(tags, fn, In, f_callback)

if (nargin < 4), f_callback = @(x) x; end

ref.tags = tags;
ref.fn = fn;
ref.In = In;
ref.f_callback = f_callback;
