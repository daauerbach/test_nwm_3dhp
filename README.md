# test_nwm_3dhp

scratch space to examine coupling coarser scale flows to high res hydrography, potentially supporting white paper on considering tradeoffs of enhanced 'variable flow hydrography'


## Concept & context

Natural systems exhibit few universally consistent 'bright lines', as even apparently clear distinctions are more accurately understood as artifacts of a particular scale of observation. 

Nonetheless, the legal frameworks that structure natural resource management and regulation are nonetheless formulated around various thresholds and categorical distinctions.

In particular, the idea of 'water typing' or 'stream typing' is a key determinant in _many_ practical decision contexts affecting various land uses.

This construct is premised on a recognition that water is not uniformly present in all geomorphic forms through which it travels, and this empirical observation is often then translated into basic 'perennial', 'intermittent', and 'ephemeral' classes.

We can feel reasonably confident that the mouth of the Columbia River is not 'ephemeral', but even the ingredients of this simple "P-I-E" are elusive: is a single day of no flow in a single year enough to designate a given portion of channel network as 'intermittent'? How to label a longer section consisting of interspersed surface and subsurface flow reaches?

None of this is new, and the point is not an easy deconstruction of unstable terms, but rather to provide the context and motivation for asking how new datasets may or may not facilitate more nuanced natural resource management that reduces tradeoffs between competing values (e.g., different forms of harvest).

A rich body of scientific research has examined the problems of flow classification and presence, from classic streamflow analyses such as [Poff & Ward 1989](https://cdnsciencepub.com/doi/abs/10.1139/f89-228) and [Poff 1996](https://pofflab.colostate.edu/wp-content/uploads/2019/08/Poff_1996_Ahydrogeographyofunregulated.pdf), to extensions that directly integrate geomorphic indicators, such as [Poff et al. 2006](https://bledsoe.engr.uga.edu/wp-content/uploads/2017/11/Poff-etal-2006-RRA-Placing-streamflow-in-geographic.pdf) and [Reidy Liermann et al. 2012](https://onlinelibrary.wiley.com/doi/abs/10.1002/Rra.1541), and more recent efforts to better address non-perennial settings such as [Sando et al. 2022](https://doi.org/10.1016/j.hydroa.2022.100138) and [Jaeger et al. 2023](https://onlinelibrary.wiley.com/doi/pdfdirect/10.1002/hyp.14813).

Seeking to supplement the understanding developed in these studies, this repo contains preliminary exploration of how flow predictions available at a continental scale might be coupled to high resolution hydrography at a 'project planning' scale. 

Such a linkage has the potential to deliver data-driven characterizations at previously unreachable extents and resolutions.

## Technical objective

"Invert" traditional basin modeling approaches that sum over drainage area pixels or reaches/subbasins to reach an outlet flow, instead disaggregating NWM (retrospective) values for a given NHDplus 1:100K COMID as the outlet flow for the 3DHP 'internal drainage' of that catchment and ‘downscaling’ throughout that local contributing network. 

Alternative ideas to explore, in some sort of vaguely increasing order of complexity:

  - A: “Backwards travel time, topology” – a minimalist combo of network position and basic linework geometry, using topological and ground distance with a uniform expectation
    - initial scripting has gotten the 3DHP routed fairly easily into directed graphs that allow all sorts of fun graph theory things, but this seems like it would be mostly useful as a sort of null model baseline to compare against other more nuanced methods
  - B: “Backwards travel time, hydraulic geometry” – same as above, but attempting to add further topographic attributes as (linear?) weights so that the per-3DHP flowline values better reflects slopes and widths (channel, valley, hillslope)
  - C: “Backwards travel, full process” – building on 2, but somehow adding lithology and soils and land cover data to infer that within-NHDplus-catchment variation
  - D: “Nested ensembling (? No idea what to call this)” - something else entirely. it’s likely to be difficult to train any other class of model on observed flows (GAM, hierarchical Bayesian, deep learning, whatever), because of the lack of training data (i.e., multiple surface flow time series at a sub-NHDplus catchment scale), but maybe use a statistical ‘bridge’ model trained on a higher resolution process model
    - not a first priority 

2024-10-29 mvp of approach A ![](f_ins_tt_unif_contrib_1k.png)
