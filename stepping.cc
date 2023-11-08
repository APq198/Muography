#include "stepping.hh"

MySteppingAction::MySteppingAction(MyEventAction * eventAction)
{
	fEventAction = eventAction;
}

MySteppingAction::~MySteppingAction()
{}

void MySteppingAction::UserSteppingAction(const G4Step * step)
{
	//G4double edep = step->GetTotalEnergyDeposit();
	//fEventAction->AddEdep(edep);
	//G4cout << "something2_step" << G4endl;
	// static G4int cntr = 0;
	// if (cntr<100000000)
	// {
	// 	cntr += 1;
	// 	// G4cout << "u";
	// }
	//if (cntr < 1000)
	//	G4cout << "Print many times but no that many (<1000)" << G4endl;
	//G4cout << cntr << ") ..." << G4endl;
	// cntr ++;
	#ifndef ONE_THREADED
		//G4cout << "a";
		//std::cout << std::flush;
		// G4cout << "uiop";
		//std::cout << "a";
	#endif

	G4Track * track = step->GetTrack();
	G4String particleName = track->GetParticleDefinition()->GetParticleName();
	if (particleName.compare(0, 2, "nu")==0 || 
		particleName.compare(0, 7, "anti_nu")==0 //|| 
		//particleName.compare(0, 5, "gamma")==0
	) {
		track->SetTrackStatus(fStopAndKill);
		//G4cout << "stopped :" << particleName << G4endl;
	} 

	#ifdef DEBUG_MODE
		else if (	particleName.compare(0, 6, "proton")!=0 &&
					particleName.compare(0, 2, "e-")!=0
		) {
			G4cout << "? --- " << particleName << G4endl;	// N14, O16, neutron, B11, alpha, C12, deuteron, C13, e+, mu+, C14, N15, N13, pi+, triton, Be8, mu-,        
		}
	#endif


	if (particleName.compare(0, 2, "mu")==0)
	{
		G4ThreeVector pos = track->GetPosition();
		if ( pos[1] < -Y_WORLD_VAL+1000*m )
		{
			//fEventAction->numIncidentMuons += 1;
			G4double energy = track->GetTotalEnergy();
			G4cout 	<< particleName << ", "
					<< energy 
					<< G4endl;
			// G4cout << particleName + ", ";
			// for (int i=0; i<10; i++)
			// 	G4cout << "u";
			track->SetTrackStatus(fStopAndKill);
		}
		G4cout << "mu!" << G4endl;
		
	}

}