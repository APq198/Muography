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

#include "globals.hh"
#ifdef ONE_THREADED
	#include <fstream>
#endif

#define INCLUDE_HELIUM 1


class AccurateGenerator
{
	public:
		AccurateGenerator(G4double E_min, G4double E_max, G4double k, G4double b, G4String particleName, int length, G4double * energies0_, G4double * fluxes0_);
		~AccurateGenerator();
		G4double generate_accurate_E();
	private:
		G4double *energies0, *fluxes0;
		G4double *lg_energies0, *lg_fluxes0;
		G4String particleName;
		G4double E_min, E_max;		// should be same for every object of this class
		G4double k, b;
		G4double NC_Psi;
		int length;
		G4double phi_interpolated(G4double lgE);
		G4double Psi(G4double);
		G4double inverse_CDF(G4double);
};	


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
		G4bool useDistribution;
		G4bool launchVertically;
		static const unsigned int numOfPoints = 29;
		//G4double energies0[29];
		//G4double fluxes0[29];
		G4double energies0[29] = {
			282518489.9667483,
			393149690.8141331,
			547102879.5511398,
			761342480.502315,
			1059682816.9221555,
			1475218971.877678,
			2230891386.7257385,
			3105696577.6373997,
			3980916300.305746,
			5102782642.757679,
			6542080141.015659,
			9109220692.97067,
			11680859003.840343,
			16267653489.799633,
			20856128956.930878,
			26741445395.695816,
			32884848574.952053,
			40455386058.19259,
			51866284781.70766,
			75284564450.58363,
			941103449842.5316,
			3543047381596.9365,
			17091140465325.418,
			36002103821286.086,
			50158847134420.58,
			89558429206492.34,
			181058861289951.53,
			297603041534441.9,
			450224202888678.4	//29
		};
		G4double fluxes0[29] = {
			4277.654199066912,
			3633.153005079787,
			3085.7568527160474,
			2620.8352196482674,
			1890.5829836907567,
			1158.323286254708,
			602.7564188857535,
			369.29709087504904,
			226.26111818250686,
			138.62577005384486,
			72.13665969982199,
			37.53773681925604,
			16.590449618912096,
			7.332435087468232,
			3.8155775399500476,
			1.8298325446794785,
			1.319980075832034,
			0.6868777248548796,
			0.3574304018221044,
			0.14558631963371407,
			0.0002495531870178878,
			0.000004954906469324441,
			8.355757272293554e-8,
			1.632172128563277e-8,
			5.203697672647922e-9,
			1.1967796184265081e-9,
			1.9855111964445867e-10,
			5.3764604311936354e-11,
			2.0182015668528597e-11
		};
		G4double phi_interpolated(G4double E);
		///// G4double my_lerp(G4double x1, G4double x2, G4double y1, G4double y2, G4double x);
		G4double E_min_mes, E_max_mes;
		G4double E_min = 100e9;//pow(10, 10);
		G4double E_max = pow(10, 14);
		G4double k = -2.6;	// from PCR_Fluxes.ipynb - manually
		G4double b = 27.7;
		G4double NC_Psi = pow(10, b) * (pow(E_max, (k+1)) - pow(E_min, (k+1))) / (k+1); //normalizing constant for Psi(x)
		G4double Psi(G4double);
		G4double inverse_CDF(G4double);
		G4double generate_accurate_E();
		AccurateGenerator * ProtonGenerator;
		AccurateGenerator * AlphaGenerator;
};

G4double my_lerp(G4double x1, G4double x2, G4double y1, G4double y2, G4double x);

#endif