# Raw Data: Tide at Sewell's Point, Norfolk, VA

All raw data for Sewell's Point is copied directly (unchanged in any way) from https://github.com/scrim-network/local-coastal-flood-risk as described in Ruckert et al (2019).

* [`1928-01-01_SewellsPoint_1942-08-25.txt`](./1928-01-01_SewellsPoint_1942-08-25.txt) and [`1943-09-15_SewellsPoint_2015-12-31.txt`](./1943-09-15_SewellsPoint_2015-12-31.txt) come from the NOAA Tides and Currents database station 8638610 (Sewells Point, VA) available at https://tidesandcurrents.noaa.gov/inventory.html?id=8638610. This data is hourly (the time zone is GMT) and is given relative to the current NOAA national tidal datum epoch (NTDE; 1983–2001). Units are meters.
* [`BRICK_NOfastDynamics_SP_20Nov2018.nc`](./BRICK_NOfastDynamics_SP_20Nov2018.nc) gives simulations using the BRICK 0.2 model with slow ice sheet dynamics (Wong et al, 2017) used in Ruckert et al (2019).
* [`BRICK_SewellsPoint_FastDynamics_20Nov2018.nc`](./BRICK_SewellsPoint_FastDynamics_20Nov2018.nc) gives simulations using the BRICK 0.2 model with fast ice sheet dynamics (Wong et al, 2017) used in Ruckert et al (2019).

## References

> Ruckert, K. L., Srikrishnan, V., & Keller, K. (2019). Characterizing the deep uncertainties surrounding coastal flood hazard projections: A case study for Norfolk, VA. Scientific Reports, 9(1), 1–12. https://doi.org/10.1038/s41598-019-47587-6

> Wong, T. E., Bakker, A. M. R., Ruckert, K., Applegate, P., Slangen, A. B. A., & Keller, K. (2017). BRICK v0.2, a simple, accessible, and transparent model framework for climate and regional sea-level projections. Geoscientific Model Development, 10(7), 2741–2760. https://doi.org/10.5194/gmd-10-2741-2017
