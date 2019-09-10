# Spectral-Based Saliency Estimation (SBSE) for the identification of features in 3D meshes

SBSE is a method for the feature curve extraction on triangle meshes. The pipeline of the proposed method is separated into two basic
steps. 

At the first step, we estimate the saliency of each vertex using spectral
analysis. The magnitude of the estimated saliency identifies if a vertex is feature
or not. Based on the geometry, we can say that the feature vertices represent
the edges of a feature curve (both crests and valleys) or corners. At the second
step, we estimate the mean curvature of the extracted features and we use
it to classify the different feature curves (if exist). Additionally, we use the
information related to the mean curvature and the saliency of each feature
curve in order to find similarities with feature curves of other models.


The code is described in more detail on this paper:

### BibTeX

@inproceedings {or.20191066,
booktitle = {Eurographics Workshop on 3D Object Retrieval},
editor = {Biasotti, Silvia and Lavoué, Guillaume and Veltkamp, Remco},
title = {{Feature Curve Extraction on Triangle Meshes}},
author = {Thompson, E. Moscoso and Arvanitis, G. and Moustakas, K. and Hoang-Xuan, N. and Nguyen, E. R. and Tran, M. and Lejemble, T. and Barthe, L. and Mellado, N. and Romanengo, C. and Biasotti, S. and Falcidieno, B.},
year = {2019},
publisher = {The Eurographics Association},
ISSN = {1997-0471},
ISBN = {978-3-03868-077-2},
DOI = {10.2312/3dor.20191066}
}
