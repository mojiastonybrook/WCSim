#ifndef WCSimTrackingAction_h
#define WCSimTrackingAction_h

#include <set>
#include "G4UserTrackingAction.hh"
#include "globals.hh"
//#include "G4TrackVector.hh"  //M. Jia: add for photon propagation
//#include "G4OpticalPhoton.hh"  //M. Jia
//#include "G4DynamicParticle.hh"  //M. Jia
//#include "G4ThreeVector.hh"  //M. Jia

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Track;
class WCSimTrackingMessenger;
//class WCSimPrimaryGeneratorAction;  // M. Jia: add for photon propagation

class WCSimTrackingAction : public G4UserTrackingAction
{
 public:
//   WCSimTrackingAction(WCSimPrimaryGeneratorAction* ); // M. Jia: modified for photon propagation.
  WCSimTrackingAction();
  ~WCSimTrackingAction();

  void PreUserTrackingAction (const G4Track* aTrack);
  void PostUserTrackingAction(const G4Track*);

  void SetFractionChPhotons(G4double fraction){percentageOfCherenkovPhotonsToDraw = fraction;}
  
private:
  std::set<G4String> ProcessList;
  std::set<G4int> ParticleList;
  std::set<G4int> pi0List;

  // TF: define in macro now
  G4float percentageOfCherenkovPhotonsToDraw;

  WCSimTrackingMessenger* messenger;
  //M. Jia: add for photon propagation.
//  WCSimPrimaryGeneratorAction* generatorAction;
//  G4TrackVector* newOpticalPhoton;
};


#endif


