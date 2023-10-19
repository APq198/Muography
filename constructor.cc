#include "constructor.hh"
#include "CADMesh.hh"		// Easy 3d models --- https://github.com/christopherpoole/CADMesh


DetectorConstruction::DetectorConstruction()
{}

DetectorConstruction::~DetectorConstruction()
{}

void DetectorConstruction::ConstructCOW(G4NistManager *nist, G4LogicalVolume* logicWorld)
{
	auto mesh = CADMesh::TessellatedMesh::FromPLY("cow.ply");	// the Cow
    G4VSolid* solidCow = mesh->GetSolid();
	G4LogicalVolume * logicCow = new G4LogicalVolume(solidCow, nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE"), "logicCow");
		// rotate around OX +90 degrees
		G4ThreeVector AxisOfRotation = G4ThreeVector(1,0,0).unit();
		G4RotationMatrix * RotMat = new G4RotationMatrix();
		RotMat -> rotate(+90.0 *deg, AxisOfRotation);
    	//RotMat -> invert();
	G4VPhysicalVolume * physCow = new G4PVPlacement(/*RotMat*/0, G4ThreeVector(0.,0.,0.), logicCow, "physCow", logicWorld, false, 0, true);
}

void DetectorConstruction::ConstructMountainFuji(G4NistManager *nist, G4LogicalVolume* logicWorld)
{
	auto mesh = CADMesh::TessellatedMesh::FromOBJ("Mount_Fuji.obj");	// the Cow
    G4VSolid* solidMountain = mesh->GetSolid();
	G4LogicalVolume * logicMountain = new G4LogicalVolume(solidMountain, nist->FindOrBuildMaterial("G4_AIR"), "logicMountain");
		// rotate around OX +90 degrees
		G4ThreeVector AxisOfRotation = G4ThreeVector(1,0,0).unit();
		G4RotationMatrix * RotMat = new G4RotationMatrix();
		RotMat -> rotate(-90.0 *deg, AxisOfRotation);
    	RotMat -> invert();
	G4VPhysicalVolume * physMaountain = new G4PVPlacement(RotMat, G4ThreeVector(), logicMountain, "physMountain", logicWorld, false, 0, true);
}

void DetectorConstruction::MyConstruct3Detectors(G4NistManager *nist, G4LogicalVolume* logicWorld)
{
	G4Box * solidDetector = new G4Box("solidDetector", 4.*m, 4.*m, .0005*m);
	logicDetector = new G4LogicalVolume(solidDetector, nist->FindOrBuildMaterial("G4_AIR"), "logicDetector");
	G4VPhysicalVolume * physDetector0 = new G4PVPlacement(0, G4ThreeVector(0.,0.,4.7*m), logicDetector, "physDetector0", logicWorld, false, 0, true);
	G4VPhysicalVolume * physDetector1 = new G4PVPlacement(0, G4ThreeVector(0.,0.,3.7*m), logicDetector, "physDetector1", logicWorld, false, 1, true);
	G4VPhysicalVolume * physDetector2 = new G4PVPlacement(0, G4ThreeVector(0.,0.,2.7*m), logicDetector, "physDetector2", logicWorld, false, 2, true);
}

void DetectorConstruction::ConstructObstacle(G4NistManager *nist, G4LogicalVolume* logicWorld)
{
	G4Material *obstacleMat = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
	G4Box * solidObstacle = new G4Box("solidObstacle", 1.5*m, 1*m, 0.75*m);
	G4LogicalVolume * logicObstacle = new G4LogicalVolume(solidObstacle, obstacleMat, "logicObstacle");
	G4VPhysicalVolume * physWorld = new G4PVPlacement(0, G4ThreeVector(3.*m, -3*m, 1*m), logicObstacle, "physObstacle", logicWorld, false, 0, true);
}

void DetectorConstruction::ConstructOmniDetector()
{
	G4Box * solidDetector = new G4Box("solidDetector", xWorld, detectorHeight/2.0, zWorld);
	logicDetector = new G4LogicalVolume(solidDetector, nist->FindOrBuildMaterial("G4_WATER"), "logicDetector");
	#ifndef USE_STANDARD_ATMOSPHERE
		G4VPhysicalVolume * physDetector0 = new G4PVPlacement(0, G4ThreeVector(0., detectorHeight/2.0 - yWorld, 0.), logicDetector, "physDetector0", logicWorld, false, 0, true);
	#else
		G4VPhysicalVolume * physDetector0 = new G4PVPlacement(0, G4ThreeVector(0., -detectorHeight -yWorld, 0.), logicDetector, "physDetector0", logicWorld, false, 32, true);
	#endif
}

void DetectorConstruction::ConstructMountain(G4NistManager *nist, G4LogicalVolume* logicWorld)
{
	int N=2;
	G4double length = 4*km;
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
		{
			G4Box * solidMountaian = new G4Box("solidMountain", 1*m, /*length/(2*N)*/(i+1)*km, 1*m);
			G4LogicalVolume * logicMountain = new G4LogicalVolume(solidMountaian, nist->FindOrBuildMaterial("G4_WATER"), "logicMountain");
			G4VPhysicalVolume * physicalMountain = new G4PVPlacement(0, G4ThreeVector((i/N)*length/2, 0., (j/N)*length/2), logicMountain, "physMountain", logicWorld, false, j+i*N, true);
		}
}


void DetectorConstruction::ConstructAsteriod_Sphere(G4NistManager *nist, G4LogicalVolume* logicWorld)
{
	G4Material * asteriodMat = nist->FindOrBuildMaterial("G4_WATER");
	G4double innerRadius = 0*m;
	G4double radius = 100*m;
	G4Sphere * solidSphere = new G4Sphere("solidSphere", innerRadius, radius, 0, 360*degree, 0, 180*degree);
	G4LogicalVolume * logicSphere = new G4LogicalVolume(solidSphere, asteriodMat, "logicSphere");
	G4VPhysicalVolume * physSphere = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicSphere, "physSphere", logicWorld, false, 0, true);
}



G4VPhysicalVolume * DetectorConstruction::ConstructAsteroidScene()
{
	G4NistManager *nist = G4NistManager::Instance();
	G4Material *worldMat = nist->FindOrBuildMaterial("G4_Galactic");
	G4Box * solidWorld = new G4Box("solidWorld", 250*m, 250*m, 250*m);
	G4LogicalVolume * logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
	G4VPhysicalVolume * physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicWorld, "physWorld", 0, false, 0, true);
	

	//ConstructObstacle(nist, logicWorld);
	//ConstructOmniDetector(nist, logicWorld);
	//ConstructCOW(nist, logicWorld);
	//ConstructMountain(nist, logicWorld);
	ConstructAsteriod_Sphere(nist, logicWorld);

	return physWorld;
}



G4VPhysicalVolume * DetectorConstruction::ConstructSurfaceScene()
{
	nist = G4NistManager::Instance();
	G4Material *worldMat = nist->FindOrBuildMaterial("G4_Galactic");

	yWorld = Y_WORLD_VAL;
	G4double xWorld = 20*km;
	G4double zWorld = 20*km;

	#ifdef USE_STANDARD_ATMOSPHERE
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
		
		

		G4Box * solidWorld = new G4Box("solidWorld", xWorld, yWorld+4*detectorHeight, zWorld);
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
			Air[i] = new G4Material("G4_AIR_" + stri.str(), density, 2 );
			Air[i]->AddElement(elN, 70*perCent);
			Air[i]->AddElement(elO, 30*perCent);
			
			solidLayer[i] = new G4Box("solidLayer", xWorld, deltaY/2.0 , zWorld);
			logicLayer[i] = new G4LogicalVolume(solidLayer[i], Air[i], "logicLayer");
			physLayer[i] = new G4PVPlacement(0, G4ThreeVector(0., y-yWorld, 0.), logicLayer[i], "physLayer", logicWorld, false, i, true);
		}
	//} else {
	#else
		//xWorld = 20*km;
		//yWorld = 1*km;
		//yWorld = Y_WORLD_VAL;
		//zWorld = 20*km;
		G4double world_height = 2 * yWorld;
		G4Box * solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);
		logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
		physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicWorld, "physWorld", 0, false, 0, true);
		
		
		const G4int numOfLayers = 10;

		G4double H = 8400*m;
		G4double density0 = 1.29 * kg/m3;
		G4double aN = 14.01 *g/mole;
		G4double aO = 16 *g/mole;
		G4Element * elN = new G4Element("Nitrogen", "N", 7, aN);
		G4Element * elO = new G4Element("Oxygen", "O", 8, aO);
		G4double degOfFreedom = 3;
		G4double R = 8.3144626181532;
		G4double gamma = (degOfFreedom + 2) / degOfFreedom;
		G4double T = 293.15;
		G4double M = (0.3*aO + 0.7*aN) / 1000.0;
		G4double g0 = 9.8;
		G4double e = 2.718281828459045;
		G4double M_air = 29 / 1000;

		G4Material * Air[numOfLayers];

		G4double zeta = 0.8;
		G4Box * solidAtmosphere = new G4Box("soliAir", xWorld, zeta*yWorld/numOfLayers, zWorld);
		G4LogicalVolume * logicAtmosphere[numOfLayers];
		G4VPhysicalVolume * physAtmosphere[numOfLayers];

		for (int layer=0; layer<numOfLayers; layer++)
		{
			G4double y = (float) yWorld * (2*zeta*layer - numOfLayers + zeta) / (numOfLayers) + detectorHeight;
			//y = yWorld/10. *2*layer - yWorld + yWorld/10. ; 

			std::stringstream stri;
			stri << layer;

			//G4double deltaH = 40e3 / (1.0*numOfLayers) * layer;
			//G4double density = density0 * pow( ( 1 - (gamma-1)/gamma * M * g0 * deltaH / (R*T)  ), (1 / (gamma - 1) ) );
			//G4double density = density0 * pow(e, (-1)*(M_air*g0)/(R*T));	// my
			G4double density = density0 * pow(e,  (-1.0) * y / H  );		// from paper (thesis)
			Air[layer] = new G4Material("G4_AIR_" + stri.str(), density, 2 );
			Air[layer]->AddElement(elN, 70*perCent);
			Air[layer]->AddElement(elO, 30*perCent);

			logicAtmosphere[layer] = new G4LogicalVolume(solidAtmosphere, Air[layer], "logicAtmosphere");
			physAtmosphere[layer] = new G4PVPlacement(0, G4ThreeVector(0., y, 0.), logicAtmosphere[layer], "physAtmosphere", logicWorld, false, layer, true);
		}
	#endif
	//}


	/*G4Material *atmosphereMaterial = nist->FindOrBuildMaterial("G4_AIR");
	G4double atmosphere_height = world_height * 0.9;
	G4Box * solidAtmosphere = new G4Box("Solid Atmosphere", xWorld, atmosphere_height/2, zWorld);
	G4LogicalVolume* logicAtmosphere = new G4LogicalVolume(solidAtmosphere, atmosphereMaterial, "logicAtmosphere");
	G4VPhysicalVolume * physAtmosphere = new G4PVPlacement(0, G4ThreeVector(0.,(atmosphere_height - world_height)/2 ,0.), logicAtmosphere, "physAtmosphere", logicWorld, false, 0, true);
	*/

	return physWorld;
}



G4VPhysicalVolume * DetectorConstruction::Construct()
{
	//G4VPhysicalVolume * physWorld = ConstructAsteroidScene();
	detectorHeight = 10*m;
	G4VPhysicalVolume * physWorld = ConstructSurfaceScene();
	ConstructOmniDetector();
	return physWorld;
}


void DetectorConstruction::ConstructSDandField()
{
	SensitiveDetector *sensDet = new SensitiveDetector("SensitiveDetector");

	if (logicDetector != nullptr)
	logicDetector->SetSensitiveDetector(sensDet);		// logical Volume of an object accepts 
														// sensitive detector to become sensitive
}