#ifndef ACTION_HH
#define ACTION_HH 1

#include "G4VUserActionInitialization.hh"

#include "generator.hh"

class ActionInitialization : public G4VUserActionInitialization
{
	public:
		ActionInitialization();
		~ActionInitialization();

		virtual void Build() const;
};

#endif