#include "action.hh"


ActionInitialization::ActionInitialization()
{}

ActionInitialization::~ActionInitialization()
{}

void ActionInitialization::BuildForMaster() const
{
	MyRunAction *runAction = new MyRunAction();
	SetUserAction(runAction);
}

void ActionInitialization::Build() const
{
	PrimaryGenerator * generator  = new PrimaryGenerator();
	SetUserAction(generator);

	MyRunAction *runAction = new MyRunAction();
	SetUserAction(runAction);

	MyEventAction *eventAction = new MyEventAction(runAction);
	SetUserAction(eventAction);

	MySteppingAction * steppingAction = new MySteppingAction(eventAction);
	SetUserAction(steppingAction);
}