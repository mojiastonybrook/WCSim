#ifndef WCSimPhotonPropagationBuilder_h
#define WCSimPhotonPropagationBuilder_h 1

#include "G4VPhysicsConstructor.hh"
#include "WCSimPhotonPropagation.hh"

class WCSimPhotonPropagationBuilder : public G4VPhysicsConstructor{
public:
	WCSimPhotonPropagationBuilder() {}
	virtual ~WCSimPhotonPropagationBuilder() {}

	void ConstructParticle();

	void ConstructProcess();
};

#endif // WCSimPhotonPropagationBuiler_h
