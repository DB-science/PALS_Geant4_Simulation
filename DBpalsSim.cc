#include "DG4LTSettingsManager.hh"
#include "DG4LTLifetimeSimulatorEngine.hh"
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"
#include "Randomize.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4SteppingVerbose.hh"
#include "G4MTRunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //
  CLHEP::HepRandom::setTheSeed((unsigned)clock());
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) { ui = new G4UIExecutive(argc, argv); }

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);
  //use G4SteppingVerboseWithUnits
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  // Construct the default run manager
  //
  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(new DetectorConstruction());
 
 
 
  G4int verbose = 0;
  
  G4VModularPhysicsList* physicsList = new FTFP_BERT(verbose);
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4(verbose));
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics(verbose);

  physicsList->RegisterPhysics(opticalPhysics);
  runManager->SetUserInitialization(physicsList);

  // Set user action classes
  runManager->SetUserInitialization(new ActionInitialization());

 
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  // Process macro or start UI session
  //
  UImanager->ApplyCommand("/tracking/verbose 0");
  UImanager->ApplyCommand("/tracking/storeTrajectory 0");
  
  
  if ( ! ui ) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }


	

 // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
  //
   if (!DG4LTSimulationParameterManager::sharedInstance()->registerFile("input.dg4lt")); // place here filename (example can be found in the folder "input.dg4lt")
		return -1;




  delete visManager;
  delete runManager;
  

 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
