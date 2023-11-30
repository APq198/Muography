#include "generator.hh"

#define PI 3.1415926535897

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


G4double my_lerp(G4double x1, G4double x2, G4double y1, G4double y2, G4double x)
{	return y1 + (x-x1)*(y2-y1)/(x2-x1);		}


PrimaryGenerator::PrimaryGenerator()
{
	fMessenger = new G4GenericMessenger(this, "/generator/", "Particle generator");
	//fMessenger->DeclareProperty("setParticleMomentum", momentum, "Change momentum of the particle (doesn't work, use energy)");
	fMessenger->DeclareProperty("setParticleEnergy", energy, "Change momentum of the particle");
	fMessenger->DeclareProperty("setParticleName", particleName, "Change name of the primary particle");
	fMessenger->DeclareProperty("useDistribution", useDistribution, "Whether to use distribution based on PCR fluxes (1-use, 0-use specified momentum)");
	fMessenger->DeclareProperty("launchVertically", launchVertically, "Whether to launch PCRs vertically (only [0,-1,0])");
	fMessenger->DeclareProperty("launchCos3Theta", cos3DistributionOverAzimuthAngle, "Whether to launch PCRs with Phi(theta)~cos3(theta)");
	fMessenger->DeclareProperty("sin2aDistribution", sin2aDistribution, "Phi~sin(2*theta)");
	
	momentum = 1*GeV;
	energy = 1*GeV;
	particleName = "proton";
	useDistribution = 1;
	launchVertically = 0;
	cos3DistributionOverAzimuthAngle = 0;
	sin2aDistribution = 1;

	fParticleGun = new G4ParticleGun(1);
}

PrimaryGenerator::~PrimaryGenerator()
{
	delete fParticleGun;
}

void PrimaryGenerator::GeneratePrimaries(G4Event * anEvent)
{
	MyGeneratePrimaries_CosmicRays_Surface(anEvent);
}


void PrimaryGenerator::MyGeneratePrimaries_CosmicRays_Surface(G4Event * anEvent)
{
	// defining particle
	G4ParticleTable * particleTable = G4ParticleTable::GetParticleTable();
	if (INCLUDE_HELIUM == 0 ) {
		G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
		fParticleGun->SetParticleDefinition(particle);
	}
	G4ParticleDefinition *proton_particle = particleTable->FindParticle("proton");
	G4ParticleDefinition *alpha_particle = particleTable->FindParticle("alpha");


	// direction
	G4double zenith_angle, azimuth_angle;
	G4ThreeVector direction;
	if (sin2aDistribution) {
		while (1) {
			zenith_angle = G4UniformRand() * PI / 2;	// theta
			if (G4UniformRand() < sin(2*zenith_angle) && zenith_angle < THETA_MAX)
				break;
		}
		azimuth_angle = 0;//G4UniformRand() * PI * 2;		// phi
		direction = G4ThreeVector( sin(zenith_angle)*cos(azimuth_angle), -cos(zenith_angle), sin(zenith_angle)*sin(azimuth_angle) ).unit();
		fParticleGun->SetParticleMomentumDirection(direction);
	} else if (cos3DistributionOverAzimuthAngle) {		//тут треба щось робити з dOmega - тілесним кутом? На кшталт ~sin(2theta)
		while (1) {
			zenith_angle = G4UniformRand() * PI / 2;	// theta
			if (G4UniformRand() < pow(cos(zenith_angle), 1) && zenith_angle < THETA_MAX)
				break;
		}
		azimuth_angle = 0;//G4UniformRand() * PI * 2;		// phi
		direction = G4ThreeVector( sin(zenith_angle)*cos(azimuth_angle), -cos(zenith_angle), sin(zenith_angle)*sin(azimuth_angle) ).unit();
		fParticleGun->SetParticleMomentumDirection(direction);
	} else if (launchVertically) {
		direction = G4ThreeVector( 0, -1, 0 ).unit();
		fParticleGun->SetParticleMomentumDirection(direction);
	} else {
		// random
		G4double theta_max = PI/2/10;
		G4double theta = G4UniformRand()*theta_max;		// angle between momentum and OY axis
		G4double phi = 0;//G4UniformRand()*2*PI;
		direction = G4ThreeVector(  sin(theta)*cos(phi),  -1*cos(theta),  sin(theta)*sin(phi)  ).unit();
		fParticleGun->SetParticleMomentumDirection(direction);
	}


	// position: 
	G4double yWorld = Y_WORLD_VAL;
	G4double x_launch = X_WORLD_VAL - R - Y_WORLD_VAL*tan(zenith_angle);
	G4ThreeVector pos(x_launch,0,0);	// .501*yWorld
	fParticleGun->SetParticlePosition(pos);


	// energy
	if (true) {
		G4double E_kin=0;
		if (G4UniformRand() < 0.1) {	// 10% всіх частинок - Гелій
			fParticleGun->SetParticleDefinition(alpha_particle);
			E_kin = AlphaGenerator->generate_accurate_E() * eV;
			G4cout << G4endl << "Launching an alpha with distributed energy, energy = " << E_kin;// << G4endl;
			fParticleGun->SetParticleEnergy(E_kin); 	// kinetic energy is actually being set	
		} else {					// інші 90% - протони
			fParticleGun->SetParticleDefinition(proton_particle);
			G4double E_kin = ProtonGenerator->generate_accurate_E() * eV;
			G4cout << G4endl << "Launching a proton with distributed energy, energy = " << E_kin;// << G4endl;
			fParticleGun->SetParticleEnergy(E_kin); 
		}
		G4cout <<  "," <<
				direction[0] << "," <<
				direction[1] << "," <<
				direction[2] <<
				G4endl;
	}

	#ifndef DONT_LAUNCH
		fParticleGun->GeneratePrimaryVertex(anEvent);
	#endif
}



// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //



AccurateGenerator_Protons::AccurateGenerator_Protons(G4double E_min, G4double E_max)
{
	//E_min = E_min;
	//E_max = E_max;
}

G4double AccurateGenerator_Protons::phi_interpolated(G4double lgE)
{
	if ( lgE < log10(energies0[0]) ){
		return 0.0;
	}	
	for(int i=0; i<length-1; i++) {
		if ( log10(energies0[i])<=lgE && lgE<log10(energies0[i+1]) ) {
			return my_lerp(log10(energies0[i]), log10(energies0[i+1]), log10(fluxes0[i]), log10(fluxes0[i+1]), lgE);
		}
	}
	return 0.0;
}

G4double AccurateGenerator_Protons::generate_accurate_E() 
{
	while(1)
	{
		G4double X3 = inverse_CDF(G4UniformRand());
		if(G4UniformRand() < pow(10, phi_interpolated(log10(X3))) / (NC_Psi*Psi(X3)))
			return X3;
	}
}



AccurateGenerator_Alpha::AccurateGenerator_Alpha(G4double E_min, G4double E_max)
{
	//E_min = E_min;
	//E_max = E_max;
}

G4double AccurateGenerator_Alpha::phi_interpolated(G4double lgE)
{
	if ( lgE < log10(energies0[0]) ){
		return 0.0;
	}	
	for(int i=0; i<length-1; i++) {
		if ( log10(energies0[i])<=lgE && lgE<log10(energies0[i+1]) ) {
			return my_lerp(log10(energies0[i]), log10(energies0[i+1]), log10(fluxes0[i]), log10(fluxes0[i+1]), lgE);
		}
	}
	return 0.0;
}

G4double AccurateGenerator_Alpha::generate_accurate_E() 
{
	while(1)
	{
		G4double X3 = inverse_CDF(G4UniformRand());
		if(G4UniformRand() < pow(10, phi_interpolated(log10(X3))) / (NC_Psi*Psi(X3)))
			return X3;
	}
}