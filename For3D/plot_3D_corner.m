function plot_3D_corner(Ni, Nj, Nk)

%% function plot_3D_corner(Ni, Nj, Nk)
%
% Add borders to help visualizing in 3D the corner removed by remove_3D_corner.
%
% -- inputs --
% Ni, Nj, Nk: position of lines in the 3D axes, corresponding to the size
% of the stack (width, height & number of images)
%
% (c) Arnauld SERGE, 2011
%
% See also remove_3D_corner, rendering_3D


Nx = Nj; Ny = Ni; Nz = Nk;

a = axis;
Ax = a(2); Ay = a(4); Az = a(6);

hold on

plot3([0 Ax Ax Ax], [Ay Ay 0 0], [Az Az Az 0], 'k')

plot3([0 0 Nx/2 Nx/2 Ax], [0 0 0 0 0], [0 Nz/2 Nz/2 Az Az], 'k')
plot3([0 0 0 0], [0 Ny/2 Ny/2 Ay], [Nz/2 Nz/2 Az Az], 'k')
plot3([Nx/2 Nx/2 0], [0 Ny/2 Ny/2], [Az Az Az], 'k')
plot3([Nx/2 Nx/2 0], [0 Ny/2 Ny/2], [Nz/2 Nz/2 Nz/2], 'k')
plot3([Nx/2 Nx/2], [Ny/2 Ny/2], [Nz/2 Az], 'k')

%%%