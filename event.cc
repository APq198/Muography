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
	G4cout << "something01" << G4endl;
}


void MyEventAction::EndOfEventAction(const G4Event * event)
{
	//G4cout << "Energy deposition: " << fEdep << G4endl;
}