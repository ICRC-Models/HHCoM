% Prepares linear indices used to access model compartments
function newVec = makeVec(varargin)
varray = ['d' , 'v' , 'h' , 's' , 'p' , 'g' , 'a' , 'r'];
vec = zeros(size(varargin));
for i = 1 : length(varargin)
    if isequal(varray(i) , varargin{i})
        varargin{i} = -1;
    end
end
vec = allcomb(varargin{1, :});
d = sym('d' , 'real');
v = sym('v' , 'real');
h = sym('h' , 'real');
s = sym('s' , 'real');
p = sym('p' , 'real');
g = sym('g' , 'real');
a = sym('a' , 'real');
r = sym('r' , 'real');
symVec = [d v h s p g a r];
replace = vec(1 , :) < 0; % index array of negatives
symVec = replace .* symVec;
vars = ones(size(vec , 1) , 1) * symVec;
vec(vec < 0) = 0;
newVec = vec + vars;
% newVec = matlabFunction(newVec);