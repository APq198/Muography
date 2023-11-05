#include <iostream>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4VPhysicalVolume.hh"
//#include "QGSP_BERT.hh"
//#include "QGSP_BIC.hh"
#include "FTFP_BERT_HP.hh" 		// <- mars article: "FTFP_BERT_HP is recommended for all cosmic ray applications by GEANT4 authors" or smth like that
//#include "QGSP_BERT_HP.hh"
//#include "FTFP_BERT_EMZ.hh"		// see choosingPhysicsLists.key

#include "G4MTRunManager.hh"


#include "constructor.hh"
#include "physics.hh"
#include "action.hh"
#include "globals.hh"


int main(int argc, char** argv)
{
	#ifdef ONE_THREADED
		G4RunManager * runManager = new G4RunManager();
	#else
		#ifdef G4MULTITHREADED
			G4MTRunManager * runManager = new G4MTRunManager();
		#else
			G4RunManager * runManager = new G4RunManager();
		#endif
	#endif
	//G4RunManager * runManager = new G4RunManager();

	runManager->SetUserInitialization(new DetectorConstruction());
	runManager->SetUserInitialization(new PhysicsList());
	runManager->SetUserInitialization(new ActionInitialization());

	G4VModularPhysicsList* physics = new FTFP_BERT_HP();
    //physics->RegisterPhysics(new G4DecayPhysics());		//
    runManager->SetUserInitialization(physics);

	//runManager->Initialize(); // /run/initialize instead
	

	G4UIExecutive *ui = new G4UIExecutive(argc, argv);

	G4UImanager *UImanager = G4UImanager::GetUIpointer();

	#ifdef ONE_THREADED
		G4VisManager *visManager = new G4VisExecutive();
		visManager->Initialize();
		UImanager -> ApplyCommand("/control/execute vis.mac");
	#else
		#ifdef G4MULTITHREADED
			#ifdef SEARCHING_FOR_WINDOW
				// cut from here 
				UImanager -> ApplyCommand("/control/execute search_for_window.mac");
				G4int N = 100;
				G4double lgE_start = 10;
				G4double lgE_end = 14;
				G4double lgE = 0;
				G4double E_MeV = 0;
				for (int i=0; i<N; i++)
				{
					lgE = ( (N-i)*lgE_start + (i)*lgE_end ) / N ;
					E_MeV = pow(10, lgE-6);	// 6 -> Mega (eV)
					UImanager -> ApplyCommand("/generator/setParticleEnergy " + std::to_string(E_MeV));	//round()
					UImanager -> ApplyCommand("/run/beamOn 4");
				}
				// to here
			#else
				UImanager -> ApplyCommand("/control/execute mt_run.mac");
			#endif
		#endif
	#endif
	ui->SessionStart();
	return 0;
}