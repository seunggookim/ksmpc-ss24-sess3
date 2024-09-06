function h = myspm_plotgifti(M)
h_axes = gca;
h = patch('Faces',M.faces, 'Vertices',M.vertices, 'LineStyle','none');
h.Parent = h_axes;
h.FaceVertexCData = M.cdata;
h.FaceColor = 'interp';
axis off image;
camlight; material dull
end