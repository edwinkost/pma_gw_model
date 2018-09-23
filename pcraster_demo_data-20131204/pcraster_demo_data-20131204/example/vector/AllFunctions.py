from pcraster import *
dem = readmap("dem.map")
Gradx = gradx(dem);
report(Gradx, "gradx.map")
Grady = grady(dem);
report(Grady, "grady.map")
Lax = lax(dem, 0.5);
report(Lax, "lax.map")
Laplacian = laplacian(dem);
report(Laplacian, "laplacian.map")
Divergence = divergence(Gradx, Grady);
report(Divergence, "divergence.map")
Diver = diver(Gradx, Grady, spatial(scalar(1)), spatial(scalar(1)));
report(Diver, "diver.map")
