# Research

## Radial Velocity

Ch. 5.1 in https://arxiv.org/pdf/1706.01629.pdf gives an expectation formula for the radial velocity in terms of 
* M: the mass of the star
* m: the mass of the exoplanet
* T: the orbital period of the exoplanet around the star
* _e_: the eccentricity of the orbit
* I: the angle between line of sight and the orbital angular momentum
* &omega;: the angle of pericenter
* &tau;: time of passage through pericenter
* v<sub>0</sub>: the mean radial velocity of the center of mass

We will combine m and M into one parameter:
* (2&pi;G)<sup>1/3</sup>m/(M+m)<sup>2/3</sup>
## Priors

The eccentricity of a closed orbit ranges from 0 (circular) to 1, not inclusive.
* _e_: [0, 1)

For the possible masses of stars https://rdcu.be/cJPdQ and planets 
* M: [0.072 M<sub>&#9737;</sub>, 150 M<sub>&#9737;</sub>]
However, not every mass is equally as likely, "for every 1,000 dwarf stars we find just one star weighing 20 solar masses, and more massive stars are rarer still." The border between massive planets and brown dwarfs in contentious, so put the upper limit on planetary mass at 24 M<sub>J</sub> for generality. 
* m: [?, 0.02291 M<sub>&#9737;</sub>]
## Likelihoods
