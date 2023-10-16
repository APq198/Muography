#include "event.hh"

MyEventAction::MyEventAction(MyRunAction * myrunaction)
{
	//fEdep = 0;
}

MyEventAction::~MyEventAction()
{}

void MyEventAction::BeginOfEventAction(const G4Event * event)
{
	//fEdep = 0.;
	/* 
	 * event 
	 *		-> GetPrimaryVertex();
	 * 				-> GetPosition();
	 * 				   і все?
	 */
	G4cout << "begin of event; " << G4endl;
}


void MyEventAction::EndOfEventAction(const G4Event * event)
{
	G4cout << "end of event; " << G4endl;
}
