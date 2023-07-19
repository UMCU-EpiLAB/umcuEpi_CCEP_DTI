function data = position2reslicedImage_mif_nosave(els,fname,number)
% input: 
% els [3x] matrix with x y z coordinates of x electrodes (native space)
% default: 
%     Copyright (C) 2009  D. Hermes, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
%
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%    
%   Version 1.1.0, released 26-11-2009
%   Changed february 2022 by Susanne Jelsma to work for .mif format imaging data and use only to save some coordinates in an empty volume

%% select resliced image for electrodes
if ~exist(fname)
    error('image for electrodes doesnt exists')
else 
    %data.Name=fname;
end
data = read_mrtrix (fname);
transform = data.transform;

% convert electrodes from native 2 indices
els_ind=round((els-repmat(transform(1:3,4),1,size(els,1))')*...
    inv(transform(1:3,1:3)'));
temp.electrode=zeros(size(data.data));

if isempty(els)
output = temp.electrode;
disp('isempty')
else
for elec=1:size(els_ind,1)
temp.electrode(els_ind(elec,1),els_ind(elec,2),els_ind(elec,3)) = number; % give the voxels a chosen value
output=temp.electrode;
end
data.data = output;
%write_mrtrix (data, outname) % no save part
end
