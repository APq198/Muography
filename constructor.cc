#include "constructor.hh"
//#include "CADMesh.hh"		// Easy 3d models --- https://github.com/christopherpoole/CADMesh

DetectorConstruction::DetectorConstruction()
{}

DetectorConstruction::~DetectorConstruction()
{}


G4VPhysicalVolume * DetectorConstruction::ConstructSurfaceScene()
{
	nist = G4NistManager::Instance();
	G4Material *worldMat = nist->FindOrBuildMaterial("G4_Galactic");

	yWorld = Y_WORLD_VAL;
	xWorld = X_WORLD_VAL;
	zWorld = R;

	G4int altitude_stand[] = {
		0,
		1000,
		2000,
		3000,
		4000,
		5000,
		6000,
		7000,
		8000,
		9000,
		10000,
		15000,
		20000,
		25000,
		30000,
		40000,
		50000,
		60000,
		70000,
		80000
	};
	G4double density_stand[] = {
		1.225,
		1.112,
		1.007,
		0.9093,
		0.8194,
		0.7364,
		0.6601,
		0.5900,
		0.5258,
		0.4671,
		0.4135,
		0.1948,
		0.08891,
		0.04008,
		0.01841,
		0.003996,
		0.001027,
		0.0003097,
		0.00008283,
		0.00001846
	};
	
	G4Box * solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);	// +detectorHeight
	logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
	physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicWorld, "physWorld", 0, false, 0, true);
	
	G4Box * solidLayer[15];
	G4LogicalVolume * logicLayer[15];
	G4VPhysicalVolume * physLayer[15];
	G4Material * Air[15];

	G4double aN = 14.01 *g/mole;
	G4double aO = 16 *g/mole;
	G4Element * elN = new G4Element("Nitrogen", "N", 7, aN);
	G4Element * elO = new G4Element("Oxygen", "O", 8, aO);

	for (G4int i=0; i < 15; i++)	// 15 - lower than 40 km (see table)
	{
		G4double deltaY = altitude_stand[i+1]*m - altitude_stand[i]*m ;
		G4double y = (altitude_stand[i+1]*m + altitude_stand[i]*m) / 2.0 ;

		std::stringstream stri;
		stri << i;
		G4double density = density_stand[i] * kg/m3;
		Air[i] = new G4Material("STANDARD_AIR_" + stri.str(), density, 2 );
		Air[i]->AddElement(elN, 79*perCent);
		Air[i]->AddElement(elO, 21*perCent);
		
		solidLayer[i] = new G4Box("solidLayer", xWorld, deltaY/2.0 , zWorld);
		logicLayer[i] = new G4LogicalVolume(solidLayer[i], Air[i], "logicLayer");
		physLayer[i] = new G4PVPlacement(0, G4ThreeVector(0., y-yWorld, 0.), logicLayer[i], "physLayer", logicWorld, false, i, true);
	}
	return physWorld;
}


G4VPhysicalVolume * DetectorConstruction::ConstructSurfaceScene_MartianAtmosphere()
{
	nist = G4NistManager::Instance();
	G4Material *worldMat = nist->FindOrBuildMaterial("G4_Galactic");

	yWorld = Y_WORLD_VAL;
	xWorld = X_WORLD_VAL;
	zWorld = R;

	G4int altitude_mars[] = {
		0,
		1000,
		2000,
		3000,
		4000,
		5000,
		6000,
		7000,
		8000,
		9000,
		10000,
		15000,
		20000,
		25000,
		30000,
		40000
	};
	G4double density_mars[] = {
		0.01740355, 	//i=0
		0.01573928, 
		0.01423417, 
		0.01287299, 
		0.01164197,
		0.01052867, 
		0.00952184, 
		0.00861128, 
		0.0077878, 
		0.00704307,
		0.00636956, 
		0.00385341, 
		0.00233121, 
		0.00141032, 
		0.0008532,		//i=14
		0.00031227		//i=15
	};
	
	G4Box * solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);	//+detectorHeight
	logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
	physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicWorld, "physWorld", 0, false, 0, true);
	
	G4Box * solidLayer[15];
	G4LogicalVolume * logicLayer[15];
	G4VPhysicalVolume * physLayer[15];
	G4Material * Air[15];

	G4double aN = 14.01 *g/mole;
	G4double aO = 16 *g/mole;
	G4double aC = 12 *g/mole;
	G4double aAr = 39.79 *g/mole;
	//G4Element * elN = new G4Element("Nitrogen", "N", 7, aN);
	G4Element * elO = new G4Element("Oxygen", "O", 8, aO);
	G4Element * elC = nist->FindOrBuildElement("C");//new G4Element("Carbon", "C", 6, aN);
	//G4Element * elAr = new G4Element("Argon", "Ar", 18, aAr);
	//G4Material * CO2 = new G4Material("Carbon Dioxide CO2", )

	for (G4int i=0; i < 15; i++)	// 15 - lower than 40 km (see table)
	{
		G4double deltaY = altitude_mars[i+1]*m - altitude_mars[i]*m ;
		G4double y = (altitude_mars[i+1]*m + altitude_mars[i]*m) / 2.0 ;

		std::stringstream stri;
		stri << i;
		G4double density = density_mars[i] * kg/m3;
		Air[i] = new G4Material("MARTIAN_AIR_" + stri.str(), density, 2 );
		//Air[i]->AddElement(elN, 79*perCent);
		//Air[i]->AddElement(elO, 21*perCent);
		Air[i]->AddElement(elC, 27.3*perCent);
		Air[i]->AddElement(elO, 72.7*perCent);
		//Air[i]->AddElement(elN, 3*perCent);
		//Air[i]->AddElement(elAr, 2*perCent);
		
		solidLayer[i] = new G4Box("solidLayer", xWorld, deltaY/2.0 , zWorld);
		logicLayer[i] = new G4LogicalVolume(solidLayer[i], Air[i], "logicLayer");
		physLayer[i] = new G4PVPlacement(0, G4ThreeVector(0., y-yWorld, 0.), logicLayer[i], "physLayer", logicWorld, false, i, true);
	}

	return physWorld;
}


G4VPhysicalVolume * DetectorConstruction::Construct()
{
	if (MARTIAN_SURFACE) 
		G4VPhysicalVolume * physWorld = ConstructSurfaceScene_MartianAtmosphere();
	else
		G4VPhysicalVolume * physWorld = ConstructSurfaceScene();
	return physWorld;
}


void DetectorConstruction::ConstructSDandField()	// SD = Sensitive Detector
{
	SensitiveDetector *sensDet = new SensitiveDetector("SensitiveDetector");

	if (logicDetector != nullptr)
	logicDetector->SetSensitiveDetector(sensDet);		// logical Volume of an object accepts 
														// sensitive detector to become sensitive
}