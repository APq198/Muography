#ifndef GENERATOR_HH
#define GENERATOR_HH 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"
#include "math.h"
#include "G4GenericMessenger.hh"
#include "constructor.hh"


class PrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
	public:
		PrimaryGenerator();
		~PrimaryGenerator();

		virtual void GeneratePrimaries(G4Event * anEvent);
	private:
		void MyGeneratePrimaries_Muons(G4Event * anEvent);
		void MyGeneratePrimaries_CosmicRays_Asteroid(G4Event * anEvent);
		void MyGeneratePrimaries_CosmicRays_Surface(G4Event * anEvent);
		G4ParticleGun * fParticleGun;
		void GenerateMuonFlux();

		G4GenericMessenger * fMessenger;

		G4double momentum;
		G4double energy;
		G4String particleName;
};

#endif