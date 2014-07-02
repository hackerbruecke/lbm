function savevtk(array, origin, spacing, filename)
% original source: http://www.exolete.com/recipes/matlab_vtk
%
% modifications: Nikola Tchipev
%
%  savevtk Save a 3-D scalar array in VTK format.
%  savevtk(array, filename) saves a 3-D array of any size to
%  filename in VTK format.

    [nx, ny, nz] = size(array);
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'file read by LBM and uniformly spaced\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
    fprintf(fid, '\n');
    fprintf(fid, 'ORIGIN    %f   %f   %f\n', origin(1), origin(2), origin(3));
    fprintf(fid, 'SPACING   %f   %f   %f\n', spacing(1), spacing(2), spacing(3));
    fprintf(fid, '\n');
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
    fprintf(fid, 'SCALARS inputfluidMask unsigned_char\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '\n');
    for c=1:nz
        for b=1:ny
            for a=1:nx
                fprintf(fid, '%d ', array(a,b,c));
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
return
