#ifndef WCSimSteppingAction_h
#define WCSimSteppingAction_h 1

#include "G4Event.hh"
#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
//#include "G4TrackVector.hh"  //M. Jia: add for photon propagation
//#include "G4OpticalPhoton.hh"  //M. Jia
//#include "G4DynamicParticle.hh"  //M. Jia
//#include "G4ThreeVector.hh"  //M. Jia

#include <vector>

class G4HCofThisEvent;
class G4Event;
//class WCSimPrimaryGeneratorAction;  //M. Jia: add for photon propagation.

class WCSimSteppingAction : public G4UserSteppingAction
{

public:
//  WCSimSteppingAction(WCSimPrimaryGeneratorAction* );
//  {};  //M. Jia: add for photon propagation.
  WCSimSteppingAction()
  { };
  ~WCSimSteppingAction()
  { };

  void UserSteppingAction(const G4Step*);

  G4int G4ThreeVectorToWireTime(G4ThreeVector *pos3d,
				G4ThreeVector lArPos,
				G4ThreeVector start,
				G4int i);
  
  void Distortion(G4double x,
		  G4double y);

  G4double FieldLines(G4double x,
		      G4double y,
		      G4int xy);


  static G4int n_photons_through_mPMTLV;
  static G4int n_photons_through_acrylic;
  static G4int n_photons_through_gel;
  static G4int n_photons_on_blacksheet;
  static G4int n_photons_on_smallPMT;

private:

  G4double ret[2];
//  WCSimPrimaryGeneratorAction* generatorAction; //M. Jia: add for photon propagation. 
//  G4TrackVector* newOpticalPhoton; 
//  std::vector<G4Track* >* newOpticalPhoton;
//  G4int primaryStepCounter;
};


#endif
