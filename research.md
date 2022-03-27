# Research

## Radial Velocity

Ch. 5.1 in https://arxiv.org/pdf/1706.01629.pdf gives an expectation formula for the radial velocity in terms of 
* _M_: the mass of the star
* _m_: the mass of the exoplanet
* _T_: the orbital period of the exoplanet around the star
* _e_: the eccentricity of the orbit
* _I_: the angle between line of sight and the orbital angular momentum
* &omega;: the angle of pericenter
* &tau;: a time of passage through pericenter
* v<sub>0</sub>: the mean radial velocity of the center of mass

We will combine m and M into one parameter:
* &mu;=(2&pi;G)<sup>1/3</sup>_m_/(_M_+_m_)<sup>2/3</sup>

## Priors

The eccentricity of a closed orbit ranges from 0 (circular) to 1, not inclusive.
* _e_: [0, 1)

For the possible masses of stars https://rdcu.be/cJPdQ and planets 
* _M_: [0.072 M<sub>&#9737;</sub>, 150 M<sub>&#9737;</sub>]

However, not every mass is equally as likely, "for every 1,000 dwarf stars we find just one star weighing 20 solar masses, and more massive stars are rarer still." The border between massive planets and brown dwarfs in contentious, so put the upper limit on planetary mass at 24 M<sub>J</sub> for generality. We can also consider the orbit of very small bodies such that the masses are almost zero. 
* _m_: [0, 0.02291 M<sub>&#9737;</sub>]

These ranges and the assumption that _M_+_m_&approx;_M_ reduce the possiblities of the mass parameter to:
* &mu;: [0, 1.246059x10<sup>6</sup>]

The angles _I_ and &omega; can be found by consulting a diagram https://en.wikipedia.org/wiki/Argument_of_periapsis#/media/File:Orbit1.svg:
* _I_: (-&pi;, &pi;]
* &omega;: [0, &pi;/2]

&tau; is any time at which the orbiting body passes through the pericenter i.e. it is closest to the central body. for simplicity take the smallest such time. this means that &tau; will have the same range as _T_.

Taking the most extreme known orbital periods of exoplanets https://arxiv.org/pdf/1808.06796.pdf, https://arxiv.org/pdf/2107.02805.pdf: 
* _T_, &tau;: [3282.3503 s, 3.46896x10<sup>13</sup> s]

Tentatively put the range on v<sub>0</sub> to be very wide until more information can be found:
* v<sub>0</sub>: [-1000 m/s, 1000 m/s]

## Likelihoods
