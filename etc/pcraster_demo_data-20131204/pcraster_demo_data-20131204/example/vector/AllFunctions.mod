gradx.map = gradx(dem.map);
grady.map = grady(dem.map);
lax.map = lax(dem.map, 0.5);
laplacian.map = laplacian(dem.map);
divergence.map = divergence(gradx.map, grady.map);
