function plotMatrix(A,mesh,cmap)
    % plotMatrix: plots a matrix using image
    % A: matrix
    % mesh: srtruct containing ys and zs
    % cmap: Colormap
    colormap(cmap);
    image(mesh.ys,mesh.zs,flipud(A),'AlphaData',0.5);
end
