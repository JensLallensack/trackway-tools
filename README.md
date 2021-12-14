# trackway-tools
R scripts to calculate common parameters from trackways of quadrupedal animals

== Dependencies ==

The following R packages need to be installed:

retistruct
DescTools
rlist
zoo
emmeans
circular
BBmisc


== Running the script ==

1) Place gaitcalculator.R and tracklib.R, as well as the trackway coordinate data, into a single folder, and set that folder as working folder in R.

2) Read-in the gaitcalculator script:
source("gaitcalculator.R")

3) Read-in the trackway coordinates:
coords  <- list(read.csv("trackway_coordinates.csv",row.names=1))

4) More trackways can be added to the list if needed, and will then be processed together.

5) Run the script. Output files will appear in the working folder:
gaitcalculator(coords)


== Support ==

Please contact Lallensack (jens.lallensack@gmail.com) for any questions.
