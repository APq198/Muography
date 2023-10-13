#ifndef DETECTOR_HH
#define DETECTOR_HH 1

#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include <fstream>

//#include "globals.hh"

class SensitiveDetector : public G4VSensitiveDetector
{
	public:
		SensitiveDetector(G4String);
		~SensitiveDetector();
	private:
		virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
		virtual G4bool ProcessHits_Old(G4Step *, G4TouchableHistory *);
};

#endif