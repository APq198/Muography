#include <iostream>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4VPhysicalVolume.hh"
#include "QGSP_BERT.hh"
//#include "QGSP_BIC.hh"
//#include "FTFP_BERT.hh"
//#include "QGSP_BERT_HP.hh"

#include "G4MTRunManager.hh"


#include "constructor.hh"
#include "physics.hh"
#include "action.hh"

int main(int argc, char** argv)
{
	#ifdef G4MULTITHREADED
		G4MTRunManager * runManager = new G4MTRunManager();
	#else
		G4RunManager * runManager = new G4RunManager();
	#endif
	//G4RunManager * runManager = new G4RunManager();

	runManager->SetUserInitialization(new DetectorConstruction());
	runManager->SetUserInitialization(new PhysicsList());
	runManager->SetUserInitialization(new ActionInitialization());

	G4VModularPhysicsList* physics = new QGSP_BERT();
    //physics->RegisterPhysics(new G4DecayPhysics());
    runManager->SetUserInitialization(physics);

	//runManager->Initialize(); //run/initialize
	

	G4UIExecutive *ui = new G4UIExecutive(argc, argv);

	G4VisManager *visManager = new G4VisExecutive();
	visManager->Initialize();

	G4UImanager *UImanager = G4UImanager::GetUIpointer();

	UImanager -> ApplyCommand("/control/execute run.mac");
	ui->SessionStart();
	return 0;
}