# GLM-Project

The data set (RoadKills.txt) is taken from a two-year study on vertebrate road kills in a National Road of southern Portugal (IP2, stretch Portalegre-Monforte, 27 km long). The surveyed road has paved verges with two lanes and a moderate amount of traffic (less than 10,000 vehicles per day). Road surroundings are dominated by cork Quercus suber and holm oak Q. rotundifolia tree stands, named montado and open land, including pastures, meadows, and fallows. The road was inspected for amphibian roadkills every two weeks between March 1995 and March 1997. Surveys were made by a car slowly (10–20 km per hour) driving along the road on the hard-shoulder. Each animal found dead was identified to species level, whenever possible, and its geographic location, on UTM coordinates, was determined with help of detailed cartography (1:2000) of horizontal and vertical road profiles and aerial photographs. All carcasses were removed from the road to avoid double counting. The response variable is the total number of amphibian fatalities per segment (TOT.N). All animals found dead on each segment were allocated to the coordinates of its middle point. A list with all available explanatory variables and the abbreviations used is given as follows.

* OPEN.L Open lands (ha)
* MONT.S Montado with shrubs (ha)
* POLIC Policulture (ha)
* D.PARK Distance to Natural Park (m)
* SHRUB Shrubs (ha)
* WAT.RES Water reservoirs (ha)
* L.WAT.C Length of water courses (km)
* L.P.ROAD Paved road length (km)
* D.WAT.COUR Distance to water courses
* OLIVE Olive grooves (ha)
* MONT Montado without shrubs (ha)
* URBAN Urban (ha)
* L.D.ROAD Dirty road length (m)
* D.WAT.RES Distance to water reservoirs
* N.PATCH Number of habitat Patches
* P.EDGE Edges perimeter
* L.SDI Landscape Shannon diversity index

We start fitting a Poisson regression model in a frequentist manner by relating the count TOT.N to the predictors. After performing model selection, we check the adequacy and show a prediction plot. As overdispersion is found, we fit a negative binomial model and a quasi-Poisson model. 

Finally, we redo the analysis in a Bayesian manner.
