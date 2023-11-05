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
	//fMessenger->DeclareProperty("setParticleMomentum", momentum, "Change momentum of the particle (doesn't work, use energy)");
	fMessenger->DeclareProperty("setParticleEnergy", energy, "Change momentum of the particle");
	fMessenger->DeclareProperty("setParticleName", particleName, "Change name of the primary particle");
	fMessenger->DeclareProperty("useDistribution", useDistribution, "Whether to use distribution based on PCR fluxes (1-use, 0-use specified momentum)");
	fMessenger->DeclareProperty("launchVertically", launchVertically, "Whether to launch PCRs vertically (only [0,-1,0])");

	momentum = 1*GeV;
	energy = 1*GeV;
	particleName = "proton";
	useDistribution = 0;
	launchVertically = 1;

	fParticleGun = new G4ParticleGun(1);
}

PrimaryGenerator::~PrimaryGenerator()
{
	delete fParticleGun;
}

G4double PrimaryGenerator::my_lerp(G4double x1, G4double x2, G4double y1, G4double y2, G4double x)
{	return y1 + (x-x1)*(y2-y1)/(x2-x1);		}

G4double PrimaryGenerator::phi_interpolated(G4double lgE)
{
	G4int l = 29;
	for(int i=0; i<l-1; i++){
		if ( lgE < 8.451046876324488 )
			return 0.0;
		if ( log10(energies0[i])<=lgE && lgE<log10(energies0[i+1]) ) {
			return my_lerp(log10(energies0[i]), log10(energies0[i+1]), log10(fluxes0[i]), log10(fluxes0[i+1]), lgE);
		}
	}
	return 0.0;
}

G4double PrimaryGenerator::Psi(G4double E) { return pow(10, b) * pow(E, k) / NC_Psi; }

G4double PrimaryGenerator::inverse_CDF(G4double y) {	return pow((y * (pow(E_max, (k+1)) - pow(E_min, (k+1))) + pow(E_min, (k+1))),  (1/(k+1)));	}

G4double PrimaryGenerator::generate_accurate_E()
{
	while(1)
	{
		G4double X3 = inverse_CDF(G4UniformRand());
		if(G4UniformRand() < pow(10, phi_interpolated(log10(X3))) / (NC_Psi*Psi(X3)))
			return X3;
	}
}

void PrimaryGenerator::GenerateMuonFlux() {}

void PrimaryGenerator::GeneratePrimaries(G4Event * anEvent)
{
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
	G4ParticleTable * particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
	fParticleGun->SetParticleDefinition(particle);

	// position: 
	G4double yWorld = Y_WORLD_VAL;
	G4ThreeVector pos(0,0.99*yWorld,0);
	fParticleGun->SetParticlePosition(pos);

	// direction
	if (launchVertically) {
		G4ThreeVector direction = G4ThreeVector( 0, -1, 0 ).unit();
		fParticleGun->SetParticleMomentumDirection(direction);
	} else {
		// random
		G4double theta_max = PI/2/10;
		G4double theta = G4UniformRand()*theta_max;		// angle between momentum and OY axis
		G4double phi = G4UniformRand()*2*PI;
		G4ThreeVector direction = G4ThreeVector(  sin(theta)*cos(phi),  -1*cos(theta),  sin(theta)*sin(phi)  ).unit();
		fParticleGun->SetParticleMomentumDirection(direction);
	}
	

	// energy
	if (useDistribution) {
		G4double E = generate_accurate_E() * eV;
		fParticleGun->SetParticleEnergy(E);
		G4cout << "Launching a " << particleName << " with distributed energy, energy = " << E << G4endl;
	} else {
		//fParticleGun->SetParticleMomentum(momentum);
		//G4cout << "Launching a " << particleName << ", momentum = " << momentum << G4endl;
		fParticleGun->SetParticleEnergy(energy);
		G4cout << "Launching a " << particleName << " with predefined(!) energy, energy = " << energy << G4endl;
	}
	#ifndef DONT_LAUNCH
		fParticleGun->GeneratePrimaryVertex(anEvent);
	#endif
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