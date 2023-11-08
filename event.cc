#include "event.hh"

MyEventAction::MyEventAction(MyRunAction * myrunaction)
{
	//fEdep = 0;
	// implement flux here // no, maybe in stepping? or detector?
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
	// G4cout << "begin of event; " << G4endl;
	numIncidentMuons = 0;
}


void MyEventAction::EndOfEventAction(const G4Event * event)
{
	// G4cout << G4endl << "end of event; " << G4endl;
	#ifdef SEARCHING_FOR_WINDOW
		G4cout << G4endl << numIncidentMuons << " # counted Muons" << G4endl;
	#endif
}
