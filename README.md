# **Clustering of Longitudinal Data via Probabilistic Models: the case of Semi-Markov Mixtures**

This project implements a clustering algorithm for longitudinal data described by categorical variables, using a probabilistic model: Mixtures of Semi-Markov Chains. 

The algorithm has been developed for sociological life trajectories, but can more generally be applied to handle any kind of longitudinal data, as long as the variables are categorical. The goal of the algorithm is to identify clusters of similar sequence observations and to model the underlying distribution of the data within each group.

For more detail, check the file `REPORT_research.pdf`.

### Usage

To use the algorithm, one needs to provide a dataset of longitudinal data, which should be a CSV file (or a txt file) with the following format:

- each row corresponds to an individual
- each column corresponds to an observed state over a time step (the timing can be different from one row to another, as long as the same time steps are used between 2 observations across the whole dataset). One must only have columns associated to the longitudinal observations, so one should get rid of columns such as "ID", for instance.

The algorithm knows how to handle sequences of varying lengths. So be careful to set all the non-existing states to NA (if for instance artificial states have been added previously to make all the sequences have equal lengths — which is a common practice when applying most of the traditional methods for clustering longitudinal data).

N.B: A function is then employed in the code to directly process the data in order to convert it to a more specific format that is adapted to the model. This format corresponds to the modeling of semi-Markov chains, where for each sequence is associated 2 lists: one for the successive visited states and another for the associated sojourn times.

**To run the algorithm, one should refer to the file `main.R`.**  This is an example template of how to use the code.

In output, one should have a list with the recap of the results (including the optimal number of clusters using the BIC criterion, the semi-Markov parameters characterising each cluster, as well as a few more analytical parameters).

## **Acknowledgments**

This project was developed as part of my 3-month research internship at the CEREMADE (*Centre de Recherche en Mathématiques de la Décision*) in june-august 2024. The internship was supervised by researcher Madalina Olteanu.

Link to Madalina Olteanu’s page: https://sites.google.com/view/madalina-olteanu

## References

- Hervé Cardot, Guillaume Lecuelle, Pascal Schlich, and Michel Visalli. Estimating finite
mixtures of semi-markov chains: an application to the segmentation of temporal sensory
data. 2018. [https://inria.hal.science/hal-02789807/]
