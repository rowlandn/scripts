function h_vect = obj2vect2(obj,field)
% OBJ2VECT converts an object array of handles to a vector of handles.
% (the field is usually a handle, h)

% Get a cell array of handles from the object.
h_cell = get(obj,field);

% Convert the cell array to a vector of handles.
if iscell(h_cell),
	h_vect = [h_cell{:}];
else
    h_vect = h_cell;
end