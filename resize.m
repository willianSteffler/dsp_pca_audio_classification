% Copyright (C) 2022 wiu_h
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {} {@var{retval} =} resize (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: wiu_h <wiu_h@NOTE-WIU>
% Created: 2022-06-28

function [retval,factor] = resize (filePath,sz)
    [y,Fs] = audioread(filePath);

    dim = 1:1:length(y);                        %espaço original
    tgt = 1:length(y)/sz:length(y);             %espaço reduzido
    s = interp1(dim,y(:,1),tgt);                %sinal interpolado
    if length(s) < sz
        s = [s,0];
    end
    retval = s;
    factor = sz/length(y);
end
