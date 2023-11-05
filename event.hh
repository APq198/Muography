#ifndef MY_EVENT_HH
#define MY_EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"

#include "G4UserRunAction.hh"
#include "run.hh"
#include "globals.hh"

class MyEventAction : public G4UserEventAction
{
	public:
		MyEventAction(MyRunAction *);
		~MyEventAction();

		virtual void BeginOfEventAction(const G4Event*);
		virtual void EndOfEventAction(const G4Event*);
		G4int numIncidentMuons = 0;

};

#endif