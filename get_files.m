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
% @deftypefn {} {@var{retval} =} get_files (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: wiu <wiu@wiu-ubuntu>
% Created: 2022-06-29

function retval=get_files (folder)
  samples=dir(folder);
  Tf=numel(samples); %total de arquivos
  Tf=Tf-2; %remove . e ..
  
  Filepath = cell(Tf,1);
  for i = 1 : Tf
      Filepath{i,1}=strcat(folder,'/',samples(i+2).name);
  end
  retval=Filepath;
end
