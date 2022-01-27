# Things to do or discuss

## To discuss

- Per [AGU style guide](https://www.agu.org/Publish-with-AGU/Publish/Author-Resources/Grammar-Style-Guide#hyphenation), sea level is not hyphenated (as best I can tell)

## To do

A number of action items are noted in the text through yellow (or red) boxes.
More substantive points are listed here.

- What are the priors we actually want to use on SLR and surge?
- Fix the period to T=50?
- At end of methods, recap what we want to achieve. Most readers will have forgotten this after reading the methods section.
- I'm comfortable talking about "multiple models" rather than "multiple priors" for sea level rise in the future, but there needs to be a couple of sentences explaining that these are models based on a subjective prior.

### Figures

- XLRM figure:
    - add Xl, L, R, and M to figure
    - do we really need this figure?
- violin plot of sea level rise in 2100:
    - make the violin plot taller
    - add larger fonts
- figure of historical surges
    - add panel labels (a) and (b)
    - Make them taller -- and panel (a) narrower
    - parallel y axes and same font sizes
    - (b): show 5-95% CI
- scenario maps:
    - It would be better here to have a figure that has $\Delta h$ on the $y$ axis and uses expected lifetime cost to color the points
    - x: $\Delta h$, y: SLR in 2100, color: Expected Lifetime Costs / Damages. 1 dot = 1 SOW. Can actually use a KNN to map averages.
    - I think that changing to four values of $\Delta h$ and showing two or three different values of $h_0$ might be useful (alternatively we could vary the discount rate $\gamma$).
    - Also need to rescale the $y$-axis by house value.
    - It would also be cool to create an interactive online version of this figure that lets you scroll over a particular dot and learn which scenario it's from!
- Trade-offs for RCP scenarios:
    - add panel labels (a) and (b)
    - add succinct summary sentence
- Plot of discretization
    - label axes
    - [??]
- Summary of subjective synthesis
    - Perhaps add a range of "optimal" and "robust" performance for different curves?
    - Perhaps add other metrics on a $y$-axis

### Text

A succinct summary of the multiple PDF problem (to go in the introduction):

> The problem of how to make decisions under deep uncertainty is an area of long-sanding need still highly active research [CITE Lempert, Ellsburg, ...]. While expected utility maximization provides an axiomatically appealing and computationally efficient appraohch in cases of shallower (i.e., a single PDF) uncertainty, this ... [??].. Here we demonstrate a workflow designed to keep [??]

A few sentences here

> Consider, as an example, the task of projecting the contributions of the West Ntarctic Ice Sheet to future sea level rise. These projections hinge on assumptions about the relative importance of physical  uncertainties driving the ice-sheet dynamics and future anthropogenic climate forcings (see discussion in Wang & Keller 2018). It is, thus far, an open challenge to resolve this deep uncertainty surrounding the [???] structure and forcings using the observational record [e.g. Srikrishman et al, 2021, reference WAIS work by Alley et al]