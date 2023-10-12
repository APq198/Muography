#include "detector.hh"

SensitiveDetector::SensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{}

SensitiveDetector::~SensitiveDetector()
{}

G4bool SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory * ROhist)
{

	G4String filename_;
	//if (filename != "")
	//	filename_ = filename;
	//else 
		filename_ = "noname_output01.csv";
	
	G4Track * track = aStep->GetTrack();
	track->SetTrackStatus(fStopAndKill);
	G4bool isMuon = aStep->GetTrack()->GetParticleDefinition()->GetParticleName()=="mu-";

	G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
	G4ThreeVector direction = aStep->GetPreStepPoint()->GetMomentumDirection();
	
	G4bool writePositions = false;
	G4bool writeDirections = true;

	if(writePositions && isMuon) {
		G4cout << "muon" << G4endl;
		std::ofstream myfile;
		myfile.open(filename_, std::ios::app);
		myfile << pos[0] << ", " << pos[1] << ", " << pos[2] << std::endl; 
		myfile.close();
	} 
	else if(writeDirections && isMuon) {
		G4cout << "Muon - writing direction" << G4endl;
		std::ofstream myfile;
		myfile.open(filename_, std::ios::app);
		myfile 	<< pos[0] << ", " 
				<< pos[1] << ", " 
				<< pos[2] << ", " 
				<< direction[0] << ", " 
				<< direction[1] << ", " 
				<< direction[2] 
				<< std::endl; 
		myfile.close();
	}


	/*
	 * aStep -> GetTrack()
	 * 				-> GetParticleDefinition() -> GetParticleName()
	 * 				-> GetParticleMomentum()
	 * 				-> SetTrackStatus(fStopAndKill)  			--- detector absorbs this(?) particle
	 * 		 -> GetPreStepPoint()
	 * 				-> GetLocalTime()							
	 * 				-> GetGlobalTime()							--- seems to work (identical?)
	 * 				-> GetPosition()							--- particle's position
	 * 				-> GetTouchable()							--- detector
	 * 					-> GetTranslation 						--- detector's position
	*/
	// go to my_proj for interesting things


	return true;
}