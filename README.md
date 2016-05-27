# Pseudo Detector
Do acceptance and smearing of particles from Monte Carlo.
Provides methods InsideAcceptance(pdgid, fourmomentum) and Smear(pdgid, fourmomentum).
Main class checks several subdetectors all deriving from base_detector and providing methods inside(pdgid, fourmomentum) and Smear(fourmomentum).
