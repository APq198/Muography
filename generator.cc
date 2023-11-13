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


G4double /*PrimaryGenerator::*/my_lerp(G4double x1, G4double x2, G4double y1, G4double y2, G4double x)
{	return y1 + (x-x1)*(y2-y1)/(x2-x1);		}


PrimaryGenerator::PrimaryGenerator()
{
	fMessenger = new G4GenericMessenger(this, "/generator/", "Particle generator");
	//fMessenger->DeclareProperty("setParticleMomentum", momentum, "Change momentum of the particle (doesn't work, use energy)");
	fMessenger->DeclareProperty("setParticleEnergy", energy, "Change momentum of the particle");
	fMessenger->DeclareProperty("setParticleName", particleName, "Change name of the primary particle");
	fMessenger->DeclareProperty("useDistribution", useDistribution, "Whether to use distribution based on PCR fluxes (1-use, 0-use specified momentum)");
	fMessenger->DeclareProperty("launchVertically", launchVertically, "Whether to launch PCRs vertically (only [0,-1,0])");
	//fMessenger->DeclareProperty("setMinEnergyForDistribution", E_min_mes, "Change E_min for distributing PCR");
	//fMessenger->DeclareProperty("setMaxEnergyForDistribution", E_max_mes, "Change E_max for distributing PCR");

	momentum = 1*GeV;
	energy = 1*GeV;
	particleName = "proton";
	useDistribution = 1;
	launchVertically = 1;
	//E_min_mes = E_min;
	//E_max_mes = E_max;

	G4double proton_energies0[] = {0.19992154476111484, 0.24048546088902575, 0.3029530702457401, 0.3816470335962363, 0.45908290114999095, 0.5783325689900707, 0.7285580872579213, 0.8763819207579796, 1.0541990878478094, 1.2680952111129404, 1.4565445825280665, 1.7520762739864797, 2.0124492094022144, 2.2071890938715106, 2.655025936014887, 3.1937285031375238, 3.841733375705851, 4.412646532266771, 5.068402076484144, 5.821608284521384, 6.686747126012868, 7.680452710314457, 8.82183111215549, 10.132827726011518, 11.638649212356182, 13.368248148594358, 15.354879703107551, 17.636741035637883, 22.21799532546859, 25.519772050509687, 29.312219935676872, 34.45476432790557, 40.499518204236054, 48.71683696105355, 59.97040237201342, 77.31293294673355, 101.9990175175311, 128.49390314387725, 154.56520984723392, 185.92636312377448, 223.6506684693146, 245.29280737931126, 295.06251517745477, 354.93045553601013, 407.6760385613912, 490.39316970062424, 563.2696254797191, 743.1220910475924, 936.1527229880699, 1179.3242743228461, 1354.5815548827918, 1960.0340241194403, 3110.546343581526, 4936.393162825431, 7833.9863054192865, 11871.307848915094, 18839.598072570694, 29898.17634696805, 57075.14078281005, 108955.53152050183, 250196.253568066, 601685.613520239, };
	G4double proton_fluxes0[] = {1492.49554505182, 1805.1850955720192, 1805.1850955720192, 1805.1850955720192, 1805.1850955720192, 1492.49554505182, 1492.49554505182, 1233.9692796398103, 1020.2242734613477, 843.5036311953877, 697.3940871117293, 576.5932649858903, 476.71725265692675, 394.1415080287651, 325.86932292753795, 222.75429519995427, 152.2680183094802, 125.89254117941584, 104.08575681597432, 71.14983758398104, 58.82544448529416, 40.21127336696097, 27.48719571845632, 18.78940568655321, 15.534752835119006, 10.619081562525311, 7.258879135601047, 4.9619476030028675, 2.804310428135217, 1.9169407765334439, 1.191480548089105, 0.8957233857950082, 0.6122890909287914, 0.346043291889239, 0.19557095110072745, 0.1005018168811834, 0.046961194001077754, 0.024132854564823503, 0.013639002490846559, 0.007708262959346138, 0.0052691326305264845, 0.0036018177927383435, 0.0024620924014946105, 0.0013914841406961597, 0.0009511759691218655, 0.0006501947796417433, 0.00036746619407366577, 0.00020767838809992954, 9.704128120312554e-05, 5.484416576120982e-05, 3.09959069042689e-05, 1.7517747448309716e-05, 4.626125296946266e-06, 1.477628358301306e-06, 3.9021534863390524e-07, 1.5075123182938445e-07, 3.981071705534953e-08, 8.692214496115054e-09, 1.2973091412086619e-09, 1.3235462760443923e-10, 1.975386920666182e-11, 2.9482561788479336e-12, };
	ProtonGenerator = new AccurateGenerator(E_min, E_max, -2.6, 27.7, "proton", 62, proton_energies0, proton_fluxes0);

	G4double alpha_energies0[] = {2.0140290107538013, 2.5371856997767424, 3.196235625594224, 3.8447491912360126, 4.624845623126356, 5.563222976028452, 7.008305095945919, 8.430284405131696, 10.140782140394839, 12.774914012550095, 15.366933502823404, 18.484871604474925, 22.235436768900176, 26.746988504051473, 32.17393035591118, 38.70199422227376, 46.554596849424726, 58.647445868231635, 70.54696414255153, 88.87198130187436, 102.07908820203875, 122.79085079465382, 141.03858049585872, 161.99807281694356, 186.07231797241715, 223.8262374360644, 269.240395942248, 309.25173178693325, 371.9986524632697, 427.28070984321295, 513.9756125755048, 590.3566132004903, 678.0884583280767, 778.8579767490382, 936.8876162366342, 1076.1168167092542, 1236.0366207593947, 1486.8274291008397, 1707.7822059672492, 1961.572678802497, 2253.0785019188766, 3112.9881668557823, 3575.603795497607, 4301.090849134134, 5674.432777858529, 7840.136096477761, 11344.413665925818, 18854.38742594138, 31335.945221677146, 39475.65386641065, 47485.231395660376, 78920.3370142788, 119592.7053878722, 181225.9770178466, 287603.0687416792, 549029.6614584042, 1320336.5117757008, };
	G4double alpha_fluxes0[] = {50.62299388226923, 74.05684692262398, 89.57233857950138, 89.57233857950138, 89.57233857950138, 74.05684692262398, 61.22890909287977, 50.62299388226923, 41.854208209341984, 34.60432918892454, 23.654460291311693, 19.557095110072783, 16.169465057924743, 11.05295141126013, 9.138387404950647, 6.246721929411474, 4.270067916167639, 2.918887732593866, 1.9952623149688666, 1.363900249084656, 0.9323204650824172, 0.6373057341948934, 0.4356426937402771, 0.2977920116300827, 0.2035615045653933, 0.13914841406961598, 0.09511759691218655, 0.053756970097504006, 0.03674661940736658, 0.025118864315095617, 0.017170486827250735, 0.011737219254280154, 0.008023203838600897, 0.004534419947930674, 0.00309959069042689, 0.002118785326128759, 0.0014483367988178699, 0.000990038706112178, 0.0006767601571680691, 0.00046261252969462664, 0.00031622776601683664, 0.0001221677348996783, 8.351012436296131e-05, 5.7084965001963284e-05, 2.667392689200085e-05, 1.2463848857844373e-05, 3.981071705534953e-06, 1.2715921145526467e-06, 4.0615859883769385e-07, 1.8978455660928563e-07, 1.297309141208662e-07, 3.425962529608272e-08, 4.227532514733751e-09, 1.116416089898579e-09, 3.56593902717452e-10, 6.437179983574068e-11, 6.567367271881885e-12, };
	AlphaGenerator = new AccurateGenerator(E_min, E_max, -2.7, 26.5, "alphaa", 57, alpha_energies0, alpha_fluxes0);

	fParticleGun = new G4ParticleGun(1);
}

PrimaryGenerator::~PrimaryGenerator()
{
	delete fParticleGun;
}


G4double PrimaryGenerator::phi_interpolated(G4double lgE)
{
	G4int l = 29;
	for(int i=0; i<l-1; i++){
		if ( lgE < 8.451046876324488 )	// min value
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
	if (! INCLUDE_HELIUM ) {
		G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
		fParticleGun->SetParticleDefinition(particle);
	}

	// position: 
	G4double yWorld = Y_WORLD_VAL;
	//G4ThreeVector pos(0,0.99*yWorld,0);
	G4ThreeVector pos(0,0.51*yWorld,0);
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
	if (INCLUDE_HELIUM) {
		G4double E_kin;
		if (G4UniformRand() < 0.1) {	// 10% всіх частинок - Гелій
			fParticleGun->SetParticleDefinition(particleTable->FindParticle("alpha"));
			E_kin = AlphaGenerator->generate_accurate_E() * eV;
			G4cout << G4endl << "Launching an alpha with distributed energy, energy = " << E_kin << G4endl;
		} else {						// інші 90% - протони
			fParticleGun->SetParticleDefinition(particleTable->FindParticle("proton"));
			E_kin = ProtonGenerator->generate_accurate_E() * eV;
			G4cout << G4endl << "Launching a proton with distributed energy, energy = " << E_kin << G4endl;
		}
		fParticleGun->SetParticleEnergy(E_kin); 	// kinetic energy is actually being set	
	} else if (useDistribution) {
		G4double E = generate_accurate_E() * eV;
		fParticleGun->SetParticleEnergy(E);
		G4cout << G4endl << "Launching a " << particleName << " with distributed energy, energy = " << E << G4endl;
	} else {
		//fParticleGun->SetParticleMomentum(momentum);
		//G4cout << "Launching a " << particleName << ", momentum = " << momentum << G4endl;
		fParticleGun->SetParticleEnergy(energy);
		G4cout << G4endl << "Launching a " << particleName << " with predefined(!) energy, energy = " << energy << G4endl;
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

// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //
// ----- ===== |||| ===== ----- //

AccurateGenerator::AccurateGenerator(G4double E_min_, G4double E_max_, G4double k_, G4double b_, G4String particleName_, int lenght_, G4double * energies0_, G4double * fluxes0_)
{
	E_min = E_min_;
	E_max = E_max_;
	k = k_;
	b = b_;
	particleName = particleName_;
	NC_Psi = pow(10, b) * (pow(E_max, (k+1)) - pow(E_min, (k+1))) / (k+1);
	length = lenght_;
	energies0 = energies0_;
	fluxes0 = fluxes0_;
}


AccurateGenerator::~AccurateGenerator() {}


G4double AccurateGenerator::phi_interpolated(G4double lgE)
{
	if ( lgE < log10(energies0[0]))//8.451046876324488 )	// min value
		return 0.0;
	for(int i=0; i<length-1; i++){
		if ( log10(energies0[i])<=lgE && lgE<log10(energies0[i+1]) ) {
			return my_lerp(log10(energies0[i]), log10(energies0[i+1]), log10(fluxes0[i]), log10(fluxes0[i+1]), lgE);
		}
	}
	return 0.0;
}


G4double AccurateGenerator::Psi(G4double E) { return pow(10, b) * pow(E, k) / NC_Psi; }


G4double AccurateGenerator::inverse_CDF(G4double y) {	return pow((y * (pow(E_max, (k+1)) - pow(E_min, (k+1))) + pow(E_min, (k+1))),  (1/(k+1)));	}


G4double AccurateGenerator::generate_accurate_E()
{
	while(1)
	{
		G4double X3 = inverse_CDF(G4UniformRand());
		if(G4UniformRand() < pow(10, phi_interpolated(log10(X3))) / (NC_Psi*Psi(X3)))
			return X3;
	}
}