#include "Detector.h"
#include "TRandom.h"
#include "math.h"
#include<iostream>

bool SubDetector::inside(int pdgid, TLorentzVector p)
{
    if (allowid.count(pdgid)==0) return false;
    std::map<float, bool>::iterator eta = --alloweta.upper_bound( fabs(p.Eta()) );
    if ( !(eta->second) ) return false;
    std::map<float, bool>::iterator pt = --allowpt.upper_bound( p.Pt() );
    if ( !(pt->second) ) return false;
    return true;
}

double SubDetector::sigma(TLorentzVector p)
{
    if (type == DetType::Calo)
    {
	return sqrt(siga*siga/p.E()+sigb*sigb)*p.E();
    }
    else if (type == DetType::Tracker) 
    {
	return sqrt(siga*siga*p.Perp2()+sigb*sigb); 
    }
    else if (type == DetType::Muon)
    {
	return siga*p.Pt();
    }
    else
    {
	return 1e10;
    }
}

TLorentzVector SubDetector::Smear(TLorentzVector p)
{
    if (type == DetType::Calo)
    {
	double newE = p.E()+rand.Gaus(0,sigma(p));
	double c = sqrt(newE*newE - p.M2() ) / p.P();
	p = TLorentzVector(c*p.Vect(), newE);
	    
    }
    else if (type == DetType::Tracker || type == DetType::Muon)
    {
	double newPt = p.Pt()+rand.Gaus(0,sigma(p));
	double newEta = acosh(sqrt(1+p.Pz()*p.Pz()/(newPt*newPt)));
	p.SetPtEtaPhiM(newPt, newEta, p.Phi(), p.M());
    }
    return p;
}

bool Detector::InsideAcceptance(int pdgid, TLorentzVector p)
{
    for (std::vector<base_detector*>::iterator it = subsystems.begin(); it != subsystems.end(); ++it)
    {
	if ((*it)->inside(pdgid,p))
	{
	    //std::cout << "Detected by " << (*it)->name << std::endl;
	    return true;
	}
    }
    return false;
}

TLorentzVector Detector::Smear(int pdgid, TLorentzVector p)
{
    std::vector<base_detector*>::iterator min = subsystems.begin();
    double minsig = 1e9;
    double sig;
    for (std::vector<base_detector*>::iterator it = subsystems.begin(); it != subsystems.end(); ++it)
    {
	sig = (*it)->sigma(p);
	if ( (*it)->inside(pdgid,p) && sig < minsig )
	{
	    minsig = sig;
	    min = it;
	}
    }
    //std::cout << "Smeared by " << (*min)->name << std::endl;
    //std::cout << "with sigma " << (*min)->sigma(p) << std::endl;
    return (*min)->Smear(p);
}

void Detector::ATLAS()
{
    
    SubDetector* sub = new SubDetector("Tracker", DetType::Tracker, 0.0005, 0.01);
    // Tracker accepts electrons
    sub->allowid.insert(11);
    sub->allowid.insert(-11);
    // Muons
    sub->allowid.insert(13);
    sub->allowid.insert(-13);
    // Protons
    sub->allowid.insert(2212);
    sub->allowid.insert(-2212);
    // Pions
    sub->allowid.insert(211);
    sub->allowid.insert(-211);
    // Charged Kaons
    sub->allowid.insert(321);
    sub->allowid.insert(-321);

    // Tracker covers |eta| < 2.5
    sub->alloweta[0] = true;
    sub->alloweta[2.5] = false;
    
    sub->allowpt[0] = false;
    sub->allowpt[0.1] = true;
    subsystems.push_back(sub);

    sub = new SubDetector("ECAL",DetType::Calo,0.1,0.007);
    // ECAL sees photons and electrons
    sub->allowid.insert(22);
    sub->allowid.insert(11);
    sub->allowid.insert(-11);
    // ECAL covers |eta| < 3.2
    sub->alloweta[0] = true; // Barrel
    sub->alloweta[1.37] = false;
    sub->alloweta[1.52] = true; // Endcap outer wheel
    sub->alloweta[2.47] = false;
    sub->alloweta[2.5] = true; // Endcap inner wheel
    sub->alloweta[3.2] = false;
    
    sub->allowpt[0] = false;
    sub->allowpt[25] = true;

    subsystems.push_back(sub);
    

    sub = new SubDetector("HCAL",DetType::Calo,0.5,0.03);
    // HCAL sees protons,
    sub->allowid.insert(2212);
    sub->allowid.insert(-2212);
    // Neutrons
    sub->allowid.insert(2112);
    sub->allowid.insert(-2112);
    // Pions
    sub->allowid.insert(211);
    sub->allowid.insert(-211);
    // Charged Kaons
    sub->allowid.insert(321);
    sub->allowid.insert(-321);
    // neutral K_L
    sub->allowid.insert(130);
    
    // HCAL covers |eta| < 3.2
    sub->alloweta[0] = true; // Barrel
    sub->alloweta[1.47] = false;
    sub->alloweta[1.52] = true; // Endcap
    sub->alloweta[3.2] = false;
    
    sub->allowpt[0] = false;
    sub->allowpt[25] = true;
    
    subsystems.push_back(sub);

    sub = new SubDetector("FCAL",DetType::Calo, 1, 0.1);
    // FCAL sees protons,
    sub->allowid.insert(2212);
    sub->allowid.insert(-2212);
    // Neutrons
    sub->allowid.insert(2112);
    sub->allowid.insert(-2112);
    // Pions
    sub->allowid.insert(211);
    sub->allowid.insert(-211);
    // Charged Kaons
    sub->allowid.insert(321);
    sub->allowid.insert(-321);
    // neutral K_L
    sub->allowid.insert(130);
    // FCAL also sees photons and electrons
    sub->allowid.insert(22);
    sub->allowid.insert(11);
    sub->allowid.insert(-11);

    sub->alloweta[3.2] = true; // FCal
    sub->alloweta[4.9] = false;

    sub->allowpt[0] = false;
    sub->allowpt[25] = true;

    subsystems.push_back(sub);
    
    sub = new SubDetector("Muon spectrometer",DetType::Muon,0.1);
    sub->allowid.insert(13);
    sub->allowid.insert(-13);

    sub->alloweta[0] = true;
    sub->alloweta[2.7] = false;

    sub->allowpt[0] = false;
    sub->allowpt[3] = true;
    
    subsystems.push_back(sub);
}

