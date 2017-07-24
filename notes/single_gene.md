## Condition specific findings

### Menadione

It's so reassuring to see that the high expression outliers make sense. E.g., RNR3 and TSA2 have menadione as an outlier;  


### YPD

beautifully, RPL18 and RPL20A have YPD as an outlier.


## Kinase specific findings

## Gene specific findings

### AGA1

AGA1 is a reporter for the mating pathway. It seems to turn on in Salt Stress when the kinase activity os repressed for any of the following kinases:
* PBS2
* HOG1
* SNF1

Two of these kinases seem to be in the same salt pathway. I wonder if this is known? do they share a common kinase? (yes they are)

The SNF1 one is novel. according to SGD:
	AMP-activated S/T protein kinase; forms a complex with Snf4p and members of the Sip1p/Sip2p/Gal83p family; required for transcription of glucose-repressed genes, thermotolerance, sporulation, and peroxisome biogenesis; regulates nucleocytoplasmic shuttling of Hxk2p; regulates filamentous growth and acts as a non-canonical GEF, activating Arf3p during invasive growth; SUMOylation by Mms21p inhibits its function and targets Snf1p for destruction via the Slx5-Slx8 Ub ligase

So SNF1 must act directly or indirectly to supress the mating pathway under salt. but it also is requred for a bunch of stresses:
* Glucose repressed genes
* Heat tolerance
* Sporulation (how is this connected with mating if at all?)
* Preoxisone biogenesis

according to david:
AGA1 has salt as an outlier in the hog1 and pbs2 strains (and in snf1 – more new biology!).

### HUG1:
* HUG1 is a reporter for DNA damage
* It looks like HUG1 is under negative control by Hog1 and positive control by Pbs2 – this may be another window into a new way of thinking about the HOG pathway.
* There are 0 papers on pubmed that mention both HUG1 and Hog1...
* Hug1 seems to be regulated by the TF Rfx1 (Jaehnig Cell Rep 2013)
* Rfx1 seems to be regualted by the kinase Dun1 (Jaehnig Cell Rep 2013)
* I found 1 paper that mentions that Hog1 interacts with Dun1, but it is half a sentence [via yeast 2 hybrid] (bilsland-marchesan mcb 2000) but paper mainly about Rck2
* Although DDR2 and RNR3 are also reporters of dna damage, they do not seem to have the same expression signature
Reportedly, RNR3 is under control of Rfx1 as well, which is evidence against the idea that HOG1's effect on hug1 goes through Rfx1 (from figure 3 Jaehnig Cell Rep 2013)


### KAR2
KAR2 has tunicamycin as outlier;




### ACT1
ACT1 goes up in response to AZC? Bizarre.

### Cell Wall Integrity

It’s cool that Kdx1 and Prm5 look very similar and they look different from Gsc2. That’s what Andres saw when inducing with drugs and overexpressing Msn2.

He also tried the promoter for YLR042C, it reports cell wall but not through Slt2, just curious if you get some intuition on what controls it from your data

### SDP1 (paralog of MSG5)
From andres data: Sdp1 goes up with Msn2-5a overexpression.
Negativly regulates SLT2p by direct dephosphorylation
Punctate formaiton after heatshock
paralog of MSG5

Sed1 would turn on with both Msn2-5a and zymolyase.
