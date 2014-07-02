function voxelise_and_write_vtk(vox_dim, STLin, VTKout)
% voxelise_and_write_vtk Voxelise an stl file and write result into .vtk file
%==========================================================================
% FILENAME:          voxelise_and_write_vtk.m
% AUTHOR:            Nikola Tchipev
% PURPOSE:           Voxelise an stl file and write result into .vtk file
%
% USAGE:             voxelise_and_write_vtk([vox_dimX vox_dimY vox_dimZ], STLin, VTKout)
%
% INPUT PARAMETERS:
%
%     vox_dim - Mandatory - 1x3 array     - List of the number of gridpoints per direction
%
%     STLin   - Mandatory - string        - STL input file name
%
%     VTKout  - Mandatory - string        - VTK output file name

disp('Assuming the number of gridpoints in all directions is >= 2')
[outputGrid, gcoX, gcoY, gcoZ] = VOXELISE(vox_dim(1), vox_dim(2), vox_dim(3), STLin);

% compute dr
dr(1) = gcoX(2) - gcoX(1);
dr(2) = gcoY(2) - gcoY(1);
dr(3) = gcoZ(2) - gcoZ(1);

% compute origin
or(1) = gcoX(1);
or(2) = gcoY(1);
or(3) = gcoZ(1);

% save files
% Note that we take the inverse of the mask ~outputGrid

k = strfind(VTKout, '.vtk');
if(length(k) == 0)
  disp('outputfile should have .vtk extension')
  disp('aborting')
  return
end
disp('writing files:')
disp(VTKout)
savevtk(~outputGrid, or, dr, VTKout)

disp(' ')


% print to screen
%%% NOTE: you can specify these variables via .dat input or read them directly from the vtk file
disp('Printing origin and gridsize')
disp('Origin: ')
disp([or(1) or(2) or(3)])
disp('Gridsize: ')
disp([dr(1) dr(2) dr(3)])

end
