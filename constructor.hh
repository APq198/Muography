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


class DetectorConstruction : public G4VUserDetectorConstruction
{
	public:
		DetectorConstruction();
		~DetectorConstruction();

		virtual G4VPhysicalVolume* Construct();
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
		G4LogicalVolume *logicWorld;
		G4LogicalVolume *logicDetector = nullptr;
		virtual void ConstructSDandField();
		//G4Material Air[10];
		G4double xWorld, yWorld, zWorld, detectorHeight;
		G4NistManager * nist;
};

#endif