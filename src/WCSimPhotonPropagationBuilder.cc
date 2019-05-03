#include "G4ProcessManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleDefinition.hh"

#include "WCSimPhotonPropagationBuilder.hh"
#include "WCSimPhotonPropagation.hh"

#include "G4PhysicsConstructorFactory.hh"
G4_DECLARE_PHYSCONSTR_FACTORY(WCSimPhotonPropagationBuilder);

void WCSimPhotonPropagationBuilder::ConstructParticle(){
	G4OpticalPhoton::OpticalPhotonDefinition();
}

void WCSimPhotonPropagationBuilder::ConstructProcess(){

	WCSimPhotonPropagation* photonPropagationProcess = new WCSimPhotonPropagation();
	G4ProcessManager * pManager = 0;

	aParticleIterator->reset();
	while( (*aParticleIterator)() ){

		G4ParticleDefinition* particle = aParticleIterator->value();
		G4String particleName = particle->GetParticleName();
		pManager = particle->GetProcessManager();
		if(!pManager) {
                  std::ostringstream o;
                  o << "Particle " << particleName << "without a Process Manager";
                  G4Exception("WCSimPhotonPropagationBuilder::ConstructProcess()","",
                       FatalException,o.str().c_str());
                   return;                 // else coverity complains for pManager use below        
		}
		if( photonPropagationProcess->IsApplicable(*particle) ){
		  pManager->AddProcess(photonPropagationProcess);
		  pManager->SetProcessOrdering(photonPropagationProcess,idxPostStep); 
                  std::cout<<"Particle "<< particleName << " now has photonpropagation process."<<"\n";
		}	
	}
}
