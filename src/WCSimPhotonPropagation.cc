#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Poisson.hh"
#include "G4EmProcessSubType.hh"

#include "G4LossTableManager.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"

#include "G4RunManager.hh"
#include "WCSimPrimaryGeneratorAction.hh"

#include "WCSimPhotonPropagation.hh"

/////////////////////////
// Class Implementation  
/////////////////////////

        //////////////////////
        // static data members
        //////////////////////

//G4bool G4Cerenkov::fTrackSecondariesFirst = false;
//G4double G4Cerenkov::fMaxBetaChange = 0.;
//G4int G4Cerenkov::fMaxPhotons = 0;

        //////////////
        // Operators
        //////////////

// G4Cerenkov::operator=(const G4Cerenkov &right)
// {
// }

        /////////////////
        // Constructors
        /////////////////

WCSimPhotonPropagation::WCSimPhotonPropagation(const G4String& processName, G4ProcessType type)
           : G4VProcess(processName, type) ,
            fTrackSecondariesFirst(false)
{
        SetProcessSubType(fCerenkov);

        thePhysicsTable = NULL;

	if (verboseLevel>0) {
           G4cout << GetProcessName() << " is created " << G4endl;
	}
}

// G4Cerenkov::G4Cerenkov(const G4Cerenkov &right)
// {
// }

        ////////////////
        // Destructors
        ////////////////

WCSimPhotonPropagation::~WCSimPhotonPropagation() 
{
	if (thePhysicsTable != NULL) {
	   thePhysicsTable->clearAndDestroy();
           delete thePhysicsTable;
	}
}

        ////////////
        // Methods
        ////////////

G4bool WCSimPhotonPropagation::IsApplicable(const G4ParticleDefinition& aParticleType)
{
    G4bool result = false;
    if (aParticleType.GetPDGCharge() != 0.0 && 
	aParticleType.GetPDGMass() != 0.0 &&
	aParticleType.GetParticleName() != "chargedgeantino" &&
	!aParticleType.IsShortLived() ) { result = true; }

    return result;
}

void WCSimPhotonPropagation::SetTrackSecondariesFirst(const G4bool state)
{
        fTrackSecondariesFirst = state;
}

//void WCSimPhotonPropagation::SetMaxBetaChangePerStep(const G4double value)
//{
//        fMaxBetaChange = value*CLHEP::perCent;
//}

//void WCSimPhotonPropagation::SetMaxNumPhotonsPerStep(const G4int NumPhotons)
//{
//        fMaxPhotons = NumPhotons;
//}

void WCSimPhotonPropagation::BuildPhysicsTable(const G4ParticleDefinition&)
{
    if (!thePhysicsTable) BuildThePhysicsTable();
}

// PostStepDoIt
// -------------
//
G4VParticleChange*
WCSimPhotonPropagation::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
	//Get the pointer of Primary Generator Action to retrieve data 	
	run = G4RunManager::GetRunManager();
        generatorAction = (WCSimPrimaryGeneratorAction*)run->GetUserPrimaryGeneratorAction();
	if(!generatorAction) {
		std::cout<<"User PrimaryGeneratorAction is not found!!!"<<"\n";
		return 0;
	}
	if(!(generatorAction->IsUsingPhoProEvtGenerator())){
		return 0;
	}
	if(aTrack.GetCreatorProcess()){
		std::cout<<"Not a Primary particle!"<<"\n";
		std::cout<<"Particle name: "<<aTrack.GetDefinition()->GetParticleName()<<"\n";
		return 0;
	}

        std::cout<< "Photon Propagation process starts!"<<"\n";
	
        aParticleChange.Initialize(aTrack);

        G4int NumPhotons = generatorAction->GetMaxPhotonNumber(); 	
	////////////////////////////////////////////////////////////////

	aParticleChange.SetNumberOfSecondaries(NumPhotons);

        if (fTrackSecondariesFirst) {
           if (aTrack.GetTrackStatus() == fAlive )
                   aParticleChange.ProposeTrackStatus(fSuspend);
        }
	
	////////////////////////////////////////////////////////////////

	for (G4int i = 0; i < NumPhotons; i++) {

		// Determine photon energy
		G4double photonEnergy;
		photonEnergy = 1.2389e-3 / (generatorAction->GetPhotonWaveLength())[i];
		
		// Photon position 
		G4ThreeVector aSecondaryPosition = (generatorAction->GetPhotonVtx())[i];

		// Photon momentum 
		G4ParticleMomentum photonMomentum((generatorAction->GetPhotonDir())[i].x(),
                                                  (generatorAction->GetPhotonDir())[i].y(),
                                                  (generatorAction->GetPhotonDir())[i].z()); 
                
		// Photon polarization
		G4double sx, sy, sz;
		sx = (generatorAction->GetPhotonPol())[i].x();
		sy = (generatorAction->GetPhotonPol())[i].y();
		sz = (generatorAction->GetPhotonPol())[i].z();

		G4ThreeVector photonPolarization(sx, sy, sz);  
		std::cout << sx << " " << sy << " " << sz << "\n"; // debug

                // Generate a new photon:
                G4DynamicParticle* aCerenkovPhoton =
                  new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), 
  					                 photonMomentum);
		aCerenkovPhoton->SetPolarization
				     (photonPolarization.x(),
				      photonPolarization.y(),
				      photonPolarization.z());

		std::cout << photonPolarization.x() << " " 
			  << photonPolarization.y() << " "
			  << photonPolarization.z() << "\n";  // debug
		std::cout << (aCerenkovPhoton->GetPolarization()).x() << " "
			  << (aCerenkovPhoton->GetPolarization()).y() << " "
			  << (aCerenkovPhoton->GetPolarization()).z() << "\n ";

		aCerenkovPhoton->SetKineticEnergy(photonEnergy);

                // Generate new G4Track object:
		G4double aSecondaryTime;
                aSecondaryTime = (generatorAction->GetPhotonTime())[i];
		
		G4Track* aSecondaryTrack = 
		new G4Track(aCerenkovPhoton,aSecondaryTime,aSecondaryPosition);

                aSecondaryTrack->SetTouchableHandle(
                                 aStep.GetPreStepPoint()->GetTouchableHandle());

                aSecondaryTrack->SetParentID(aTrack.GetTrackID());

		aParticleChange.AddSecondary(aSecondaryTrack);
	}

	if (verboseLevel>0) {
	   G4cout <<"\n Exiting from G4Cerenkov::DoIt -- NumberOfSecondaries = "
	          << aParticleChange.GetNumberOfSecondaries() << G4endl;
	}

//	Kill primary particle:
	aParticleChange.ProposeTrackStatus(fStopAndKill);
        return pParticleChange;
}

// BuildThePhysicsTable for the Cerenkov process
// ---------------------------------------------
//

void WCSimPhotonPropagation::BuildThePhysicsTable()
{
	if (thePhysicsTable) return;

	const G4MaterialTable* theMaterialTable=
	 		       G4Material::GetMaterialTable();
	G4int numOfMaterials = G4Material::GetNumberOfMaterials();

	// create new physics table
	
	thePhysicsTable = new G4PhysicsTable(numOfMaterials);

	// loop for materials

	for (G4int i=0 ; i < numOfMaterials; i++)
	{
	        G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector = 0;

		// Retrieve vector of refraction indices for the material
		// from the material's optical properties table 

		G4Material* aMaterial = (*theMaterialTable)[i];

		G4MaterialPropertiesTable* aMaterialPropertiesTable =
				aMaterial->GetMaterialPropertiesTable();

		if (aMaterialPropertiesTable) {

		   aPhysicsOrderedFreeVector = new G4PhysicsOrderedFreeVector();
		   G4MaterialPropertyVector* theRefractionIndexVector = 
		    	   aMaterialPropertiesTable->GetProperty("RINDEX");

		   if (theRefractionIndexVector) {
		
		      // Retrieve the first refraction index in vector
		      // of (photon energy, refraction index) pairs 

                      G4double currentRI = (*theRefractionIndexVector)[0];

		      if (currentRI > 1.0) {

			 // Create first (photon energy, Cerenkov Integral)
			 // pair  

                         G4double currentPM = theRefractionIndexVector->
                                                 Energy(0);
			 G4double currentCAI = 0.0;

			 aPhysicsOrderedFreeVector->
			 	 InsertValues(currentPM , currentCAI);

			 // Set previous values to current ones prior to loop

			 G4double prevPM  = currentPM;
			 G4double prevCAI = currentCAI;
                	 G4double prevRI  = currentRI;

			 // loop over all (photon energy, refraction index)
			 // pairs stored for this material  

                         for (size_t ii = 1;
                              ii < theRefractionIndexVector->GetVectorLength();
                              ++ii)
			 {
                                currentRI = (*theRefractionIndexVector)[ii];
                                currentPM = theRefractionIndexVector->Energy(ii);

				currentCAI = 0.5*(1.0/(prevRI*prevRI) +
					          1.0/(currentRI*currentRI));

				currentCAI = prevCAI + 
					     (currentPM - prevPM) * currentCAI;

				aPhysicsOrderedFreeVector->
				    InsertValues(currentPM, currentCAI);

				prevPM  = currentPM;
				prevCAI = currentCAI;
				prevRI  = currentRI;
			 }

		      }
		   }
		}

	// The Cerenkov integral for a given material
	// will be inserted in thePhysicsTable
	// according to the position of the material in
	// the material table. 

	thePhysicsTable->insertAt(i,aPhysicsOrderedFreeVector); 

	}
}

// GetMeanFreePath
// ---------------
//

G4double WCSimPhotonPropagation::GetMeanFreePath(const G4Track&,
                                           G4double,
                                           G4ForceCondition*)
{
        return 1.01*mm ;
}

G4double WCSimPhotonPropagation::PostStepGetPhysicalInteractionLength(
                                           const G4Track& aTrack,
                                           G4double,
                                           G4ForceCondition* condition)
{
	std::cout<<"From Photon propagation: Interaction length is returned!"<<"\n";
//        *condition = NotForced;
//        G4double StepLimit = DBL_MAX;

//        const G4Material* aMaterial = aTrack.GetMaterial();
//	G4int materialIndex = aMaterial->GetIndex();

	// If Physics Vector is not defined no Cerenkov photons
	//    this check avoid string comparison below
//	if(!(*thePhysicsTable)[materialIndex]) { return StepLimit; }

//        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
//        const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();

//        G4double kineticEnergy = aParticle->GetKineticEnergy();
//        const G4ParticleDefinition* particleType = aParticle->GetDefinition();
//        G4double mass = particleType->GetPDGMass();

        // particle beta
//        G4double beta = aParticle->GetTotalMomentum() /
//	                aParticle->GetTotalEnergy();
        // particle gamma
//        G4double gamma = aParticle->GetTotalEnergy()/mass;
//
//        G4MaterialPropertiesTable* aMaterialPropertiesTable =
//                            aMaterial->GetMaterialPropertiesTable();

//        G4MaterialPropertyVector* Rindex = NULL;

//        if (aMaterialPropertiesTable)
//                     Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");

//        G4double nMax;
//        if (Rindex) {
//           nMax = Rindex->GetMaxValue();
//        } else {
//           return StepLimit;
//        }

//        G4double BetaMin = 1./nMax;
//        if ( BetaMin >= 1. ) return StepLimit;

//        G4double GammaMin = 1./std::sqrt(1.-BetaMin*BetaMin);

//        if (gamma < GammaMin ) return StepLimit;

//        G4double kinEmin = mass*(GammaMin-1.);

//        G4double RangeMin = G4LossTableManager::Instance()->
//                                                   GetRange(particleType,
//                                                            kinEmin,
//                                                            couple);
//        G4double Range    = G4LossTableManager::Instance()->
//                                                   GetRange(particleType,
//                                                            kineticEnergy,
//                                                            couple);

//        G4double Step = Range - RangeMin;
//        if (Step < 1.*um ) return StepLimit;

//        if (Step > 0. && Step < StepLimit) StepLimit = Step; 

        // If user has defined an average maximum number of photons to
        // be generated in a Step, then calculate the Step length for
        // that number of photons. 
 
//        if (fMaxPhotons > 0) {

           // particle charge
//           const G4double charge = aParticle->
//                                   GetDefinition()->GetPDGCharge();

//	   G4double MeanNumberOfPhotons = 
//                    GetAverageNumberOfPhotons(charge,beta,aMaterial,Rindex);

//           Step = 0.;
//           if (MeanNumberOfPhotons > 0.0) Step = fMaxPhotons /
//                                                 MeanNumberOfPhotons;

//           if (Step > 0. && Step < StepLimit) StepLimit = Step;
//        }

        // If user has defined an maximum allowed change in beta per step
//       if (fMaxBetaChange > 0.) {

//          G4double dedx = G4LossTableManager::Instance()->
//                                                   GetDEDX(particleType,
//                                                           kineticEnergy,
//                                                           couple);

//           G4double deltaGamma = gamma - 
//                                 1./std::sqrt(1.-beta*beta*
//                                                 (1.-fMaxBetaChange)*
//                                                 (1.-fMaxBetaChange));

//           Step = mass * deltaGamma / dedx;

//           if (Step > 0. && Step < StepLimit) StepLimit = Step;

//        }

        *condition = StronglyForced;

	G4double StepLimit = 1.01*nm ;
        return StepLimit;
}

// GetAverageNumberOfPhotons
// -------------------------
// This routine computes the number of Cerenkov photons produced per
// GEANT-unit (millimeter) in the current medium. 
//             ^^^^^^^^^^

//G4double 
//WCSimPhotonPropagation::GetAverageNumberOfPhotons(const G4double charge,
//                              const G4double beta, 
//			      const G4Material* aMaterial,
//			      G4MaterialPropertyVector* Rindex) const
//{
//	const G4double Rfact = 369.81/(eV * cm);
//
//        if(beta <= 0.0)return 0.0;
//
//        G4double BetaInverse = 1./beta;
//
	// Vectors used in computation of Cerenkov Angle Integral:
	// 	- Refraction Indices for the current material
	//	- new G4PhysicsOrderedFreeVector allocated to hold CAI's
// 
//	G4int materialIndex = aMaterial->GetIndex();
//
	// Retrieve the Cerenkov Angle Integrals for this material  
//
//	G4PhysicsOrderedFreeVector* CerenkovAngleIntegrals =
//	(G4PhysicsOrderedFreeVector*)((*thePhysicsTable)(materialIndex));
//
//        if(!(CerenkovAngleIntegrals->IsFilledVectorExist()))return 0.0;
//
	// Min and Max photon energies 
//	G4double Pmin = Rindex->GetMinLowEdgeEnergy();
//	G4double Pmax = Rindex->GetMaxLowEdgeEnergy();
//
	// Min and Max Refraction Indices 
//	G4double nMin = Rindex->GetMinValue();	
//	G4double nMax = Rindex->GetMaxValue();
//
	// Max Cerenkov Angle Integral 
//	G4double CAImax = CerenkovAngleIntegrals->GetMaxValue();
//
//	G4double dp, ge;
//
	// If n(Pmax) < 1/Beta -- no photons generated 
//
//	if (nMax < BetaInverse) {
//		dp = 0;
//		ge = 0;
//	} 
//
	// otherwise if n(Pmin) >= 1/Beta -- photons generated  
//
//	else if (nMin > BetaInverse) {
//		dp = Pmax - Pmin;	
//		ge = CAImax; 
//	} 
//
	// If n(Pmin) < 1/Beta, and n(Pmax) >= 1/Beta, then
	// we need to find a P such that the value of n(P) == 1/Beta.
	// Interpolation is performed by the GetEnergy() and
	// Value() methods of the G4MaterialPropertiesTable and
	// the GetValue() method of G4PhysicsVector.  
//
//	else {
//		Pmin = Rindex->GetEnergy(BetaInverse);
//		dp = Pmax - Pmin;
//
		// need boolean for current implementation of G4PhysicsVector
		// ==> being phased out
//		G4bool isOutRange;
//		G4double CAImin = CerenkovAngleIntegrals->
//                                  GetValue(Pmin, isOutRange);
//		ge = CAImax - CAImin;
//
//		if (verboseLevel>0) {
//			G4cout << "CAImin = " << CAImin << G4endl;
//			G4cout << "ge = " << ge << G4endl;
//		}
//	}
	
	// Calculate number of photons 
//	G4double NumPhotons = Rfact * charge/eplus * charge/eplus *
//                                 (dp - ge * BetaInverse*BetaInverse);
//
//	return NumPhotons;		
//}
