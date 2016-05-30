#ifndef DETECTOR_H
#define DETECTOR_H

#include "TLorentzVector.h"
#include "TRandom.h"
#include<vector>
#include<map>
#include<set>
#include<string>
#include<math.h>

enum class Geometry {ATLAS, Custom};
enum class DetType {Tracker, Calo, Muon, NONE};

class base_detector
{
 public:
    base_detector(std::string s="DetName") { name = s; }
    virtual ~base_detector() {};
    std::string name;
    virtual bool inside(int, TLorentzVector) =0;
    virtual TLorentzVector Smear(TLorentzVector p) { return p; }
    virtual double sigma(TLorentzVector p) { return 1e10; };
};

class SubDetector: public base_detector
{
 public:
 SubDetector(std::string s="DetName", DetType t = DetType::NONE, double a = 1, double b = 0.1): base_detector(s), type(t), siga(a), sigb(b)
    {
	alloweta[0] = true;
	allowpt[0] = false;
	rand = TRandom(11027);
    }
    DetType type;
    TRandom rand;
    double siga; // Stochastic term in uncertainty
    double sigb; // Constant term in uncertainty
    std::set<int> allowid;
    std::map<float, bool> alloweta;
    std::map<float, bool> allowpt;
    bool inside(int, TLorentzVector);
    double sigma(TLorentzVector);
    TLorentzVector Smear(TLorentzVector);
};


class Detector
{
 public:
    Detector(Geometry geom = Geometry::ATLAS)
    {
	if (geom == Geometry::ATLAS) ATLAS();
	else if (geom == Geometry::Custom) { }
    }
    ~Detector()
    {
	for (std::vector<base_detector*>::iterator it = subsystems.begin(); it != subsystems.end(); ++it)
	{
	    delete *it;
	}
    }
    bool InsideAcceptance(int, TLorentzVector);
    double Efficiency(int pdgid, TLorentzVector p) { return InsideAcceptance(pdgid, p); }
    TLorentzVector Smear(int, TLorentzVector);

    std::vector<base_detector*> subsystems;
 private:
    void ATLAS();
};


#endif
