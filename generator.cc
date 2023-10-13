#include "generator.hh"

#define PI 3.1415926535897
#define R 4.7*m
#define H 4.*m


/*
 * Particles:
 * GenericIon, He3, alpha, anti_neutron
 * anti_nu_e, anti_nu_mu, anti_nu_tau, anti_proton
 * chargedgeantino, deuteron, e+, e-
 * eta, eta_prime, gamma, geantino
 * mu+, mu-, neutron, nu_e
 * nu_mu, nu_tau, opticalphoton, pi+
 * pi-, pi0, proton, tau+
 * tau-, triton, 
*/

PrimaryGenerator::PrimaryGenerator()
{
	fMessenger = new G4GenericMessenger(this, "/generator/", "Particle generator");
	fMessenger->DeclareProperty("setParticleMomentum", momentum, "Change momentum of the particle");
	fMessenger->DeclareProperty("setParticleEnergy", energy, "Change momentum of the particle");

	momentum = 1*GeV;
	energy = 1*GeV;

	fParticleGun = new G4ParticleGun(1);
}

PrimaryGenerator::~PrimaryGenerator()
{
	delete fParticleGun;
}

void PrimaryGenerator::GenerateMuonFlux()
{
	
}

void PrimaryGenerator::GeneratePrimaries(G4Event * anEvent)
{
	//MyGeneratePrimaries_CosmicRays_Asteroid(anEvent);
	MyGeneratePrimaries_CosmicRays_Surface(anEvent);
}

void PrimaryGenerator::MyGeneratePrimaries_Muons(G4Event * anEvent)
{
	// for random see Randomize.hh or CLHEP random
	// remember: if f(x) is probability density ( works on [a,b] ), then we can get random
	//	variable with this PDF by: f( G4UniformRand * (b-a) + a )
	// 	*(b-a) to stretch, +a to shift

	// defining particle
	G4ParticleTable * particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "mu-";
	G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
	fParticleGun->SetParticleDefinition(particle);


	// random position
	G4double r = sqrt(G4UniformRand()) * R;
	G4double angle = G4UniformRand()*2*PI;
	G4ThreeVector pos(r*cos(angle), H, r*sin(angle)); // x = r * cos(t) ; y = r * sin(t)
	fParticleGun->SetParticlePosition(pos);


	// random direction
	G4double theta = G4UniformRand()*PI/2;
	G4double phi = G4UniformRand()*2*PI;
	G4ThreeVector direction = G4ThreeVector(  sin(theta)*cos(phi),  -1*cos(theta),  sin(theta)*sin(phi)  ).unit();
	fParticleGun->SetParticleMomentumDirection(direction);


	// momentum 
	// simple model: p = p_0 * cos3(theta)
	G4double p0 = 1*GeV;
	G4double p = p0 * pow(cos(theta), 3);
	fParticleGun->SetParticleMomentum(p);
	//fParticleGun->SetParticleEnergy(50 * MeV);


	fParticleGun->GeneratePrimaryVertex(anEvent); 
}

void PrimaryGenerator::MyGeneratePrimaries_CosmicRays_Surface(G4Event * anEvent)
{
	// defining particle
	//std::cout << "Launching a proton, momentum = " << momentum << std::endl;
	G4cout << "Launching a proton, momentum = " << momentum << G4endl;
	G4ParticleTable * particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "proton";
	G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
	fParticleGun->SetParticleDefinition(particle);

	// position: 
	G4ThreeVector pos(0,0.98*40*km,0);
	fParticleGun->SetParticlePosition(pos);


	// direction
	// make it random
	//G4ThreeVector direction = G4ThreeVector(0, -1, 0).unit();
	//fParticleGun->SetParticleMomentumDirection(direction);
	G4double theta_max = PI/2/10;
	G4double theta = G4UniformRand()*theta_max;		// angle between momentum and OY axis
	G4double phi = G4UniformRand()*2*PI;
	G4ThreeVector direction = G4ThreeVector(  sin(theta)*cos(phi),  -1*cos(theta),  sin(theta)*sin(phi)  ).unit();
	fParticleGun->SetParticleMomentumDirection(direction);
	

	// energy
	//G4double momentum = 1*TeV;
	fParticleGun->SetParticleMomentum(momentum);



	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGenerator::MyGeneratePrimaries_CosmicRays_Asteroid(G4Event * anEvent)
{
	// defining particle
	G4ParticleTable * particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "proton";
	G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
	fParticleGun->SetParticleDefinition(particle);

	// position: y=2*radius
	G4double radius = 100*m;	// same as in constructor.cc	// do better
	G4double y = 2*radius;
	G4ThreeVector pos(0, y, 0);
	fParticleGun->SetParticlePosition(pos);

	// direction
	// make it random
	//G4ThreeVector direction = G4ThreeVector(0, -1, 0).unit();
	//fParticleGun->SetParticleMomentumDirection(direction);
	G4double theta_max = PI/2/10;
	G4double theta = G4UniformRand()*theta_max;		// angle between momentum and OY axis
	G4double phi = G4UniformRand()*2*PI;
	G4ThreeVector direction = G4ThreeVector(  sin(theta)*cos(phi),  -1*cos(theta),  sin(theta)*sin(phi)  ).unit();
	fParticleGun->SetParticleMomentumDirection(direction);
	

	// energy
	G4double momentum = 1*MeV;
	fParticleGun->SetParticleMomentum(momentum);



	fParticleGun->GeneratePrimaryVertex(anEvent); 
}