# HMRF-EM
A MATLAB implementation of the HMRF as described in "Segmentation of Brain MR Images Through a Hidden Markov Random Field Model and the Expectation-Maximization Algorithm" (Zhang et al., 2001). The HMRF is applied to segment images from the cross-sectional OASIS-brains dataset but the code provided can be modified for any 3D image segmentation.

What can you test this algorithm on?

Any 3D image, but I've validated the model using the OASIS- Cross-sectional dataset. This dataset consists of 416 normal and early onset Alzheimer Disease subjects from ages 18-96. Ground truth labels are provided in the dataset - which come from the FAST-FSL implementation of the HMRF described in the Zhang et al. paper.

Associated blog post here:
http://tanyanair.github.io/HMRF-EM/
