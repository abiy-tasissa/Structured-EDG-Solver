# Structured-EDG-Solver
This is a convex algorithm for the structured Euclidean Distance Geometry Problem. Formally, we have a set of m anchor nodes and p mobile nodes. We have partial distance information between the anchor nodes, and between the anchor and mobile nodes. There is no distance information between the mobile nodes. The goal is to
estimate the configuration of the mobile nodes and anchor nodes, from the partial distance information. These set of codes were developed as part of a research project on structured Euclidean Distance Geometry Problem. They were written by Abiy Tasissa and Samuel Lichtenberg. 

## MATLAB files description
`demo_structured_edg_protein.m`: This is the main file which loads different kinds of protein from the data folder and reads the (x,y,z) coordinates of the data. Using these coordinates, the full distance matrix for the data is constructed. Then, based on the sampling scheme, some entries of the distance matrix are selected. The script calls the main algorithm 'AB_completion.m' to estimate the configuration of the points.

`demo_structured_edg_synthetic.m`: This is similar file to 'demo_structured_edg_protein.m', except it is based on synthetic data. 

`AB_completion.m`: This is the main algorithm for the structured Euclidean Distance Geometry Problem for the case of exact partial information. 


## List of protein data
* `1k.off`: Sphere with 1002 points
* `cow.off`: Cow with 2601 points
* `ptswiss.mat`: Swiss roll data with 2048 points
* `UScities.mat`: Data of 2920 US cities

The first and second data sets are taken from [here](http://visionair.ge.imati.cnr.it/ontologies/shapes/search.jsp).
The third data was obtained by simply plotting the parametric equations of a Swiss roll. The last data uses Latitude
and Longitude of US cities, in different zip codes, to generate the point coordinates. 

## Instructions

To try the algorithm, you can try either of the starting scripts `demo_structured_edg_synthetic` or demo_structured_edg_protein'. You can choose a protein from the list of proteins in the data, or use a protein of interest. Set the parameters alpha, gamma and m, and choice of anchors.  

## References

Samuel Lichtenberg and Abiy Tasissa, “Localization from Structured Distance Matrices via Low-Rank Matrix Recovery”, IEEE Transactions on Information Theory, 2024.
[arXiv link](https://arxiv.org/pdf/2311.18076). 

## Feedback

Email your feedback to <a href="mailto:abiy19@gmail.com">Abiy Tasissa</a>.

