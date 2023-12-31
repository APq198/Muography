#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include <G4Sphere.hh>

#include "detector.hh"

#define MARTIAN_SURFACE 0	//if

//#define USE_STANDARD_ATMOSPHERE 1	// always

#ifndef USE_STANDARD_ATMOSPHERE
	#define Y_WORLD_VAL 40*km		// for usual (not standard atmosphere)
#else
	#define Y_WORLD_VAL 40*km		// for standard atmosphere 
#endif
//#define LONG_BOX 1			//if	// always
#define R 40*km
#define X_WORLD_VAL 0.5*(Y_WORLD_VAL*tan(THETA_MAX)+R)


class DetectorConstruction : public G4VUserDetectorConstruction
{
	public:
		DetectorConstruction();
		~DetectorConstruction();

		virtual G4VPhysicalVolume* Construct();
		G4double yWorld;
	private:
		void ConstructMountainFuji(G4NistManager *, G4LogicalVolume*);
		void ConstructCOW(G4NistManager *, G4LogicalVolume*);
		void MyConstruct3Detectors(G4NistManager *, G4LogicalVolume*);
		void ConstructObstacle(G4NistManager *, G4LogicalVolume*);
		void ConstructOmniDetector();
		void ConstructMountain(G4NistManager *, G4LogicalVolume*);
		void ConstructAsteriod_Sphere(G4NistManager *, G4LogicalVolume*);
		void ConstructSimpleAtmosphere(G4NistManager *, G4LogicalVolume*);
		G4VPhysicalVolume * ConstructAsteroidScene();
		G4VPhysicalVolume * ConstructSurfaceScene();
		G4VPhysicalVolume * ConstructSurfaceScene_MartianAtmosphere();
		G4LogicalVolume *logicWorld;
		G4VPhysicalVolume * physWorld;
		G4LogicalVolume *logicDetector = nullptr;
		virtual void ConstructSDandField();
		//G4Material Air[10];
		G4double xWorld, zWorld;
		G4NistManager * nist;
};

#endif