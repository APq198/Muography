#ifndef MY_STEPPING_ACTION_HH
#define MY_STEPPING_ACTION_HH

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"

#include "constructor.hh"
#include "event.hh"

class MySteppingAction : public G4UserSteppingAction
{
	public: 
		MySteppingAction(MyEventAction*);
		~MySteppingAction();
		virtual void UserSteppingAction(const G4Step*);
	private:
		MyEventAction *fEventAction;
};


#endif