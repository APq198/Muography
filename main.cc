#include <iostream>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4VPhysicalVolume.hh"
//#include "QGSP_BERT.hh"
//#include "QGSP_BIC.hh"
#include "FTFP_BERT.hh"


#include "constructor.hh"
#include "physics.hh"
#include "action.hh"

int main(int argc, char** argv)
{
	G4RunManager * runManager = new G4RunManager();
	runManager->SetUserInitialization(new DetectorConstruction());
	runManager->SetUserInitialization(new PhysicsList());
	runManager->SetUserInitialization(new ActionInitialization());

	G4VModularPhysicsList* physics = new FTFP_BERT();
    //physics->RegisterPhysics(new G4DecayPhysics());
    runManager->SetUserInitialization(physics);

	runManager->Initialize();
	

	G4UIExecutive *ui = new G4UIExecutive(argc, argv);

	G4VisManager *visManager = new G4VisExecutive();
	visManager->Initialize();

	G4UImanager *UImanager = G4UImanager::GetUIpointer();

	//UImanager -> ApplyCommand("/control/execute init_vis.mac");
	UImanager->ApplyCommand("/vis/open OGL 600x600-0+0");
	UImanager->ApplyCommand("/vis/drawVolume");
	UImanager->ApplyCommand("/vis/viewer/set/autoRefresh true");
	UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
	UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
	ui->SessionStart();
	return 0;
}