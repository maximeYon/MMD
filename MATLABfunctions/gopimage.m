function [h,im]=gopimage(x,y,z,s)
    
% GOPIMAGE(Z) displays a colour image of the complex matrix Z.
%   The argument is encoded by the gop color table and the
%   magnitude by the intensity. The magnitude is scaled so that the
%   maximum magnitude corresponds to intensity one.
%   If the global variable IMAGEGAMMA exists in the workspace, the
%   image is gamma corrected by that value.
%
% GOPIMAGE(Z), where Z is an MxNx2 real array displays a vector
%   field analogously.
%
% GOPIMAGE(X,Y) does the same thing, but here the components are
%   split into two matrices.
%
% GOPIMAGE(Z,M) and GOPIMAGE(X,Y,M), where M is a scalar, scales
%   the magnitude so that the value M gives maximum intensity.
%
% GOPIMAGE(...,'invert') blends the color determined by the
%   argument with white instead of black, giving images that are
%   better suited for slides.
%
% H=GOPIMAGE(...) returns a handle to the image object.
%
% [H,IM]=GOPIMAGE(...) returns a handle to the image object and the
%   image itself. 
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden

global IMAGEGAMMA	% Global gamma value, defined in workspace
  
% Check the input parameters.
error(nargchk(1,4,nargin))

% Should the image be inverted?
invert=0;
nin=nargin;
switch nin
 case 2
  if ischar(y)
    invert=1;
  end
 case 3
  if ischar(z)
    invert=1;
  end
 case 4
  if ischar(s)
    invert=1;
  end
end
if invert
  nin=nin-1;
end
    
% Is a max level provided?
M=-1;
switch nin
 case 2
  if prod(size(y))==1
    M=y;
    nin=1;
  end
 case 3
  if prod(size(z))==1
    M=z;
    nin=2;
  else
    error('The third parameter must be a scalar.')
  end
end
M=1.0/M; % Convert to normalization factor

if nin==1
  if ndims(x)==2
    u=real(x);
    v=imag(x);
  elseif ndims(x)==3
    if size(x,3)~=2
      error('Z must have a third dimension of size 2.')
    else
      u=x(:,:,1);
      v=x(:,:,2);
    end
  else
    error('Too many dimensions.')
  end
else % nin==2
  if size(x)~=size(y)
    error('X and Y does not have the same size')
  elseif ndims(x)>2
    error('Too many dimensions.')
  else
    u=x;
    v=y;
  end
end

% Get gamma value
if exist('IMAGEGAMMA','var') & isa(IMAGEGAMMA,'double') & ~ ...
      isempty(IMAGEGAMMA) 
  g = IMAGEGAMMA(1);
else
  g = 1.0;
end

w=u+i*v;
gtab=goptab; % Load gop color table
gtab=[gtab;gtab(1,:)]; % Make it cyclic to allow interpolation for all angles
tabangles=2*pi*(0:256)/256;

magw=abs(w);
if (M < 0)
  M = 1.0/max(max(magw)); % Autonormalization
end
magw=min(1,M*magw).^g; % Gamma correct luminance
argw=pi+angle(-w);
argim=reshape(interp1(tabangles,gtab,argw(:)),[size(argw) 3]);
if invert
  gopim=argim.*repmat(magw,[1 1 3])+repmat(1-magw,[1 1 3]);
else
  gopim=argim.*repmat(magw,[1 1 3]);
end
hh=image(gopim);

if nargout>0
  h=hh;
end

if nargout>1
  im=gopim;
end

function gopcoltable=goptab() 
% 
% Generates a GOP colortable 
% 
% by matsa & knutte 
% Dept. of Biomedical Engineering 
% Linköping University, Sweden 
% matsa@imt.liu.se  knutte@imt.liu.se 
% 
gopcoltable=[ 
93 234 63 
91 233 66 
89 232 69 
86 231 73 
84 230 76 
82 228 80 
80 227 83 
78 226 87 
76 225 91 
74 223 95 
72 222 99 
70 220 103 
69 219 107 
67 217 111 
65 216 115 
64 214 119 
62 212 123 
61 210 128 
59 208 132 
58 206 136 
57 204 141 
56 202 145 
55 200 149 
54 197 154 
54 195 158 
53 193 162 
53 190 166 
53 187 171 
53 185 175 
53 182 179 
53 179 183 
54 176 187 
55 173 191 
55 170 195 
56 167 199 
57 164 202 
59 161 206 
60 158 209 
61 155 212 
63 151 216 
64 148 219 
66 145 222 
68 142 225 
70 138 227 
72 135 230 
73 132 232 
75 128 235 
77 125 237 
79 121 239 
82 118 241 
84 115 243 
86 112 244 
88 108 246 
90 105 247 
92 102 249 
95 99 250 
97 95 251 
99 92 252 
101 89 253 
104 86 253 
106 83 254 
108 80 254 
111 77 255 
113 74 255 
115 71 255 
118 68 255 
120 65 255 
122 63 255 
125 60 254 
127 57 254 
130 54 253 
132 52 252 
135 49 252 
137 47 251 
140 44 250 
142 42 248 
145 39 247 
147 37 246 
150 34 244 
152 32 243 
155 30 241 
158 28 239 
160 26 237 
163 23 235 
165 21 233 
168 19 230 
171 17 228 
173 16 225 
176 14 222 
178 12 220 
181 10 217 
183 9 214 
186 7 210 
188 5 207 
191 4 204 
193 3 200 
195 1 196 
198 0 193 
200 0 189 
202 0 185 
204 0 181 
206 0 177 
208 0 173 
210 0 169 
212 0 164 
214 0 160 
215 0 156 
217 0 151 
218 0 147 
220 0 142 
221 0 138 
222 0 133 
223 0 129 
225 0 124 
226 0 120 
227 0 116 
227 0 111 
228 0 107 
229 0 103 
229 0 99 
230 0 94 
230 0 90 
231 0 86 
231 0 82 
231 0 79 
232 0 75 
232 0 71 
232 0 68 
232 0 64 
232 1 61 
232 2 58 
232 4 55 
232 5 51 
232 7 49 
232 8 46 
232 10 43 
232 11 40 
232 13 38 
232 15 36 
232 16 33 
232 18 31 
232 20 29 
232 22 27 
232 24 25 
232 26 23 
232 28 22 
232 30 20 
232 33 19 
232 35 17 
233 37 16 
233 39 14 
233 42 13 
233 44 12 
234 47 11 
234 50 10 
235 52 9 
235 55 8 
235 58 7 
236 61 7 
237 63 6 
237 66 5 
238 70 5 
238 73 4 
239 76 4 
240 79 3 
241 82 3 
241 86 2 
242 89 2 
243 93 2 
244 96 2 
245 100 1 
245 104 1 
246 108 1 
247 111 1 
248 115 1 
249 119 0 
250 123 0 
250 128 0 
251 132 0 
252 136 0 
252 140 0 
253 144 0 
254 149 0 
254 153 0 
254 157 0 
255 161 0 
255 166 0 
255 170 0 
255 174 0 
255 178 0 
255 183 0 
254 187 0 
254 191 0 
253 195 0 
253 199 0 
252 203 0 
251 206 0 
250 210 0 
249 213 0 
247 217 0 
246 220 0 
244 223 0 
242 226 0 
240 229 0 
238 232 0 
236 234 0 
234 237 0 
231 239 0 
229 241 0 
226 243 1 
224 245 1 
221 246 1 
218 248 1 
215 249 1 
212 250 1 
209 251 2 
206 252 2 
203 253 2 
200 253 3 
197 254 3 
193 254 4 
190 255 4 
187 255 5 
184 255 5 
181 255 6 
177 255 7 
174 255 7 
171 255 8 
168 254 9 
165 254 10 
162 254 11 
159 253 12 
156 253 13 
153 252 14 
150 252 15 
147 251 17 
144 250 18 
141 250 20 
138 249 21 
135 248 23 
132 248 24 
130 247 26 
127 246 28 
124 245 30 
122 245 32 
119 244 34 
117 243 36 
114 242 39 
112 241 41 
109 240 43 
107 240 46 
104 239 48 
102 238 51 
100 237 54 
97 236 57 
95 235 60 
]/255; 

