function validatefields(Struct, Fields)
% validatefields(Struct, Fields)
% test if FIELDS are defined in STRUCT

for i = 1:numel(Fields)
    assert(isfield(Struct,Fields{i}), 'Field "%s" NOT DEFINED!', Fields{i});
end
end