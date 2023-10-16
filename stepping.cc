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
	static G4int cntr = 0;
	//G4cout << cntr << ") ..." << G4endl;
	cntr ++;
}

