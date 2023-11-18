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
	//fMessenger->DeclareProperty("setMinEnergyForDistribution", E_min_mes, "Change E_min for distributing PCR");
	//fMessenger->DeclareProperty("setMaxEnergyForDistribution", E_max_mes, "Change E_max for distributing PCR");
	fMessenger->DeclareProperty("sin2aDistribution", sin2aDistribution, "Phi~sin(2*theta)");
	
	momentum = 1*GeV;
	energy = 1*GeV;
	particleName = "proton";
	useDistribution = 1;
	launchVertically = 0;
	cos3DistributionOverAzimuthAngle = 0;
	sin2aDistribution = 1;
	//E_min_mes = E_min;
	//E_max_mes = E_max;

	// G4double proton_energies0[] = {199921544.76111484, 240485460.88902575, 302953070.2457401, 381647033.5962363, 459082901.149991, 578332568.9900707, 728558087.2579212, 876381920.7579796, 1054199087.8478093, 1268095211.1129405, 1456544582.5280664, 1752076273.9864798, 2012449209.4022143, 2207189093.8715105, 2655025936.0148873, 3193728503.1375237, 3841733375.705851, 4412646532.26677, 5068402076.484144, 5821608284.521383, 6686747126.012868, 7680452710.314457, 8821831112.155489, 10132827726.011518, 11638649212.356182, 13368248148.594358, 15354879703.107552, 17636741035.637882, 22217995325.46859, 25519772050.509686, 29312219935.676872, 34454764327.90557, 40499518204.23605, 48716836961.05355, 59970402372.01342, 77312932946.73355, 101999017517.5311, 128493903143.87724, 154565209847.23392, 185926363123.77448, 223650668469.3146, 245292807379.31125, 295062515177.4548, 354930455536.01013, 407676038561.3912, 490393169700.62427, 563269625479.7191, 743122091047.5924, 936152722988.0698, 1179324274322.8462, 1354581554882.7917, 1960034024119.4404, 3110546343581.526, 4936393162825.432, 7833986305419.286, 11871307848915.094, 18839598072570.695, 29898176346968.05, 57075140782810.05, 108955531520501.83, 250196253568066.0, 601685613520239.0};
	// G4double proton_fluxes0[] = {1492.49554505182, 1805.1850955720192, 1805.1850955720192, 1805.1850955720192, 1805.1850955720192, 1492.49554505182, 1492.49554505182, 1233.9692796398103, 1020.2242734613477, 843.5036311953877, 697.3940871117293, 576.5932649858903, 476.71725265692675, 394.1415080287651, 325.86932292753795, 222.75429519995427, 152.2680183094802, 125.89254117941584, 104.08575681597432, 71.14983758398104, 58.82544448529416, 40.21127336696097, 27.48719571845632, 18.78940568655321, 15.534752835119006, 10.619081562525311, 7.258879135601047, 4.9619476030028675, 2.804310428135217, 1.9169407765334439, 1.191480548089105, 0.8957233857950082, 0.6122890909287914, 0.346043291889239, 0.19557095110072745, 0.1005018168811834, 0.046961194001077754, 0.024132854564823503, 0.013639002490846559, 0.007708262959346138, 0.0052691326305264845, 0.0036018177927383435, 0.0024620924014946105, 0.0013914841406961597, 0.0009511759691218655, 0.0006501947796417433, 0.00036746619407366577, 0.00020767838809992954, 9.704128120312554e-05, 5.484416576120982e-05, 3.09959069042689e-05, 1.7517747448309716e-05, 4.626125296946266e-06, 1.477628358301306e-06, 3.9021534863390524e-07, 1.5075123182938445e-07, 3.981071705534953e-08, 8.692214496115054e-09, 1.2973091412086619e-09, 1.3235462760443923e-10, 1.975386920666182e-11, 2.9482561788479336e-12};
	// ProtonGenerator = new AccurateGenerator(E_min, E_max, -2.6, 27.7, "proton", 62, proton_energies0, proton_fluxes0);

	//AlphaGenerator = new AccurateGenerator(E_min, E_max, -2.7, 26.5, "alphaa", 57, alpha_energies0, alpha_fluxes0);

	fParticleGun = new G4ParticleGun(1);
}

PrimaryGenerator::~PrimaryGenerator()
{
	delete fParticleGun;
}


// G4double PrimaryGenerator::phi_interpolated(G4double lgE)
// {
// 	G4int l = 29;
// 	for(int i=0; i<l-1; i++){
// 		if ( lgE < 8.451046876324488 )	// min value
// 			return 0.0;
// 		if ( log10(energies0[i])<=lgE && lgE<log10(energies0[i+1]) ) {
// 			return my_lerp(log10(energies0[i]), log10(energies0[i+1]), log10(fluxes0[i]), log10(fluxes0[i+1]), lgE);
// 		}
// 	}
// 	return 0.0;
// }

// G4double PrimaryGenerator::Psi(G4double E) { return pow(10, b) * pow(E, k) / NC_Psi; }

// G4double PrimaryGenerator::inverse_CDF(G4double y) {	return pow((y * (pow(E_max, (k+1)) - pow(E_min, (k+1))) + pow(E_min, (k+1))),  (1/(k+1)));	}

// G4double PrimaryGenerator::generate_accurate_E()
// {
// 	while(1)
// 	{
// 		G4double X3 = inverse_CDF(G4UniformRand());
// 		if(G4UniformRand() < pow(10, phi_interpolated(log10(X3))) / (NC_Psi*Psi(X3)))
// 			return X3;
// 	}
// }

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
	if (INCLUDE_HELIUM == 0 ) {
		G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
		fParticleGun->SetParticleDefinition(particle);
	}
	G4ParticleDefinition *proton_particle = particleTable->FindParticle("proton");
	G4ParticleDefinition *alpha_particle = particleTable->FindParticle("alpha");

	// position: 
	G4double yWorld = Y_WORLD_VAL;
	//G4ThreeVector pos(0,0.99*yWorld,0);
	G4ThreeVector pos(0,0.51*yWorld,0);
	fParticleGun->SetParticlePosition(pos);

	// direction
	G4double zenith_angle, azimuth_angle;
	G4ThreeVector direction;
	if (sin2aDistribution) {
		while (1) {
			zenith_angle = G4UniformRand() * PI / 2;	// theta
			if (G4UniformRand() < sin(2*zenith_angle))
				break;
		}
		azimuth_angle = G4UniformRand() * PI / 2;		// phi
		direction = G4ThreeVector( sin(zenith_angle)*cos(azimuth_angle), -cos(zenith_angle), sin(zenith_angle)*sin(azimuth_angle) ).unit();
		fParticleGun->SetParticleMomentumDirection(direction);
	} else if (cos3DistributionOverAzimuthAngle) {		//тут треба щось робити з dOmega - тілесним кутом? На кшталт ~sin(2theta)
		while (1) {
			zenith_angle = G4UniformRand() * PI / 2;	// theta
			if (G4UniformRand() < pow(cos(zenith_angle), 3))
				break;
		}
		azimuth_angle = G4UniformRand() * PI / 2;		// phi
		direction = G4ThreeVector( sin(zenith_angle)*cos(azimuth_angle), -cos(zenith_angle), sin(zenith_angle)*sin(azimuth_angle) ).unit();
		fParticleGun->SetParticleMomentumDirection(direction);
	} else if (launchVertically) {
		direction = G4ThreeVector( 0, -1, 0 ).unit();
		fParticleGun->SetParticleMomentumDirection(direction);
	} else {
		// random
		G4double theta_max = PI/2/10;
		G4double theta = G4UniformRand()*theta_max;		// angle between momentum and OY axis
		G4double phi = G4UniformRand()*2*PI;
		direction = G4ThreeVector(  sin(theta)*cos(phi),  -1*cos(theta),  sin(theta)*sin(phi)  ).unit();
		fParticleGun->SetParticleMomentumDirection(direction);
	}

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
	// } else if (useDistribution) {
	// 	G4double E = generate_accurate_E() * eV;
	// 	fParticleGun->SetParticleEnergy(E);
	// 	G4cout << G4endl << "Launching a " << particleName << " with distributed energy, energy = " << E << G4endl;
	// } else {
	// 	//fParticleGun->SetParticleMomentum(momentum);
	// 	//G4cout << "Launching a " << particleName << ", momentum = " << momentum << G4endl;
	// 	fParticleGun->SetParticleEnergy(energy);
	// 	G4cout << G4endl << "Launching a " << particleName << " with predefined(!) energy, energy = " << energy << G4endl;
	// }
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

// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //

// AccurateGenerator::AccurateGenerator(G4double E_min_, G4double E_max_, G4double k_, G4double b_, G4String particleName_, int lenght_, G4double energies0_[], G4double fluxes0_[])
// {
// 	E_min = E_min_;
// 	E_max = E_max_;
// 	k = k_;
// 	b = b_;
// 	particleName = particleName_;
// 	NC_Psi = pow(10, b) * (pow(E_max, (k+1)) - pow(E_min, (k+1))) / (k+1);
// 	length = lenght_;
// 	energies0 = energies0_;
// 	fluxes0 = fluxes0_;
// 	// G4double *lg_energies0_ = new G4double[length];
// 	// for(int i=0; i<length; i++) {
// 	// 	lg_energies0_[i] = log10(energies0[i]);
// 	// }
// 	// lg_energies0 = lg_energies0_;
// }


// AccurateGenerator::~AccurateGenerator() {}


// G4double AccurateGenerator::phi_interpolated(G4double lgE)
// {
// 	//eV! lgE -= 9;	//in GeV
// 	if ( lgE < log10(energies0[0]) ){//8.451046876324488 )	// min value
// 		return 0.0;
// 	}	
// 	for(int i=0; i<length-1; i++) {
// 		if ( log10(energies0[i])<=lgE && lgE<log10(energies0[i+1]) ) {
// 			G4cout << G4endl << "pAp " << i << G4endl;
// 			return my_lerp11(log10(energies0[i]), log10(energies0[i+1]), log10(fluxes0[i]), log10(fluxes0[i+1]), lgE);
// 		}
// 	}
// 	return 0.0;
// }

// // G4double AccurateGenerator::phi_interpolated(G4double lgE)
// // {
// // 	if ( lgE < log10(energies0[0]) )//8.451046876324488 )	// min value
// // 		return 0.0;
// // 	for(int i=0; i<length-1; i++) {
// // 		if ( lg_energies0[i]<=lgE && lgE<lg_energies0[i+1] ) {
// // 			return my_lerp(log10(energies0[i]), log10(energies0[i+1]), log10(fluxes0[i]), log10(fluxes0[i+1]), lgE);
// // 		}
// // 	}
// // 	return 0.0;
// // }


// G4double AccurateGenerator::Psi(G4double E) { return pow(10, b) * pow(E, k) / NC_Psi; }


// G4double AccurateGenerator::inverse_CDF(G4double y) {	return pow((y * (pow(E_max, (k+1)) - pow(E_min, (k+1))) + pow(E_min, (k+1))),  (1/(k+1)));	}


// G4double AccurateGenerator::generate_accurate_E()
// {
// 	int cntr = 0;
// 	while(1)
// 	{
// 		G4double X3 = inverse_CDF(G4UniformRand());
// 		if(G4UniformRand() < pow(10, phi_interpolated(log10(X3))) / (NC_Psi*Psi(X3)))
// 			return X3;
// 		cntr += 1;
// 		if (cntr > 100)
// 			G4cout << G4endl << X3 << G4endl;
// 			G4cout << G4endl << log10(X3) << "!" << G4endl;
// 			G4cout << G4endl << phi_interpolated(log10(X3)) << "?" << G4endl;
// 			G4cout << G4endl << pow(10, phi_interpolated(log10(X3))) / (NC_Psi*Psi(X3)) << "???" << G4endl;
// 			G4cout << G4endl << pow(10, phi_interpolated(log10(X3))) << ",,, " << (NC_Psi*Psi(X3)) << G4endl;
// 			cntr = cntr % 100;
// 	}
// }



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