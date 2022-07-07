% Copyright (C) 2022 wiu
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
% @deftypefn {} {@var{retval} =} menor_amostra (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: wiu <wiu@wiu-ubuntu>
% Created: 2022-06-29

function [sz,idx] = menor_amostra (files)
  
  Filepath = files;
  Tf = length(Filepath);
  menorAmostra = 1000000000;
  indice_menor_amostra = 0;
  for i = 1 : Tf
      [y,Fs] = audioread(Filepath{i,1});
      %size(y)
      if (size(y,1)<menorAmostra)
          menorAmostra = size(y,1)
          indice_menor_amostra = i;
      end
  end
  sz = menorAmostra;
  idx = indice_menor_amostra;
end
