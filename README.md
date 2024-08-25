# Structured-EDG-Solver
This is a convex algorithm designed for the structured Euclidean Distance Geometry Problem. In this problem, we are given a set of $m$ anchor nodes and $n$ mobile nodes, along with partial distance data between the anchor nodes and between the anchor and mobile nodes. There is no distance information available between the mobile nodes. The objective is to estimate the positions of both the mobile and anchor nodes based on the available partial distance information. These codes were developed as part of a research project on the structured Euclidean Distance Geometry problem by Samuel Lichtenberg and Abiy Tasissa.

## Sampling scheme
We sample the pairwise distance matrix of the anchors using the Bernoulli model, where $\gamma$ represents the probability that any given entry is selected. By design, the diagonal entries of this matrix are zero and are considered known. Additionally, there is a central anchor node from which distances to all other anchors and mobile nodes are known. In addition, for each mobile node, aside from the central node, distance information is available from $\alpha-1$ anchors, chosen uniformly at random.

## Choice of anchors
The selection of the number of anchors and which points are designated as anchors is crucial, as a "poor" choice can result in sub-optimal performance of the algorithm. Ideally, this decision should be guided by domain knowledge. In the absence of such knowledge, a simple approach is to randomly select a subset of the points as anchors.

## MATLAB files description
`demo_structured_edg_protein.m`: This is the main script that loads various proteins from the data folder and extracts their $(x,y,z)$ coordinates. Using these coordinates, it constructs the complete distance matrix for the considered protein. Following this, the script selects specific entries from the distance matrix based on the sampling scheme discussed above. Finally, it calls the main algorithm, `AB_completion.m` to estimate the configuration of the points and calculate the error. 

`demo_structured_edg_synthetic.m`: This file is similar to `demo_structured_edg_protein.m` but it operates on synthetic data instead.

`AB_completion.m`: This is the main algorithm for the structured Euclidean Distance Geometry Problem for the case of exact partial information. 

'pdb2mat.m': This reads .pdb files and converts it to MATLAB files. This script is obtained from the following [link](https://www.mathworks.com/matlabcentral/fileexchange/42957-read-and-write-pdb-files-using-matlab?s_tid=FX_rc2_behav) 

## List of protein data
* `1ptq`
* `1ax8`
* `1w2e`
* `5wov`
* `2lum`

All the protein data are downloaded from the Protein Data Bank [here](https://www.rcsb.org/).

## Dependencies

* MATLAB: 2022a
* CVX: 2.2
* MOSEK: 10.2.1

## Instructions

To test the algorithm, you can start with either of the scripts `demo_structured_edg_synthetic` or `demo_structured_edg_protein`. You can select a protein from the available list in the dataset or use a protein of your choice. Be sure to set the parameters $\alpha$, $\gamma$, $m$, and the selection of anchors accordingly. 

## Citation

If you use our code or find our paper useful and relevant, we would appreciate if you cite our paper. 
> Samuel Lichtenberg and Abiy Tasissa, “Localization from Structured Distance Matrices via Low-Rank Matrix Recovery”, IEEE Transactions on Information Theory, 2024.
[arXiv link](https://arxiv.org/pdf/2311.18076). 

## Feedback

If you have any questions about the code or feedback, email <a href="mailto:abiy19@gmail.com">Abiy Tasissa</a>.

