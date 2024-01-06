//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

/// \file avion_esferico_glare.cc
/// \brief Main program of the avion_esferico_glare example

#include "G4Types.hh"

#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"    

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"


//
int main(int argc,char** argv) {

  // Instantiate G4UIExecutive if there are no arguments (interactive mode)
  // G4UIExecutive* ui = nullptr;
  // if ( argc == 1 ) {
    // ui = new G4UIExecutive(argc, argv);
  //}

  // Run manager
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  //runManager->SetNumberOfThreads(4);
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  // UserInitialization classes (mandatory)
  
  //  Create Detector
  DetectorConstruction* detector = new DetectorConstruction;
  detector->SetPhantomSize(G4ThreeVector(1*m, 1*m, 1*m)); //Default
  detector->SetNumberOfSegmentsInPhantom(100, 100, 100); //Default
  //  If your machine does not have enough memory,
  //  please try to reduce Number of Segements in phantom.
  //  detector->SetNumberOfSegmentsInPhantom(1,1,100);  // For small memory size.
  

  runManager->SetUserInitialization(detector);
  
  runManager->SetUserInitialization(new FTFP_BERT());
  

  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // UserAction classes
  runManager->SetUserInitialization(new ActionInitialization);

  //Initialize G4 kernel
  runManager->Initialize();
      
  //get the pointer to the User Interface manager 
  
  G4UImanager * pUI = G4UImanager::GetUIpointer();  


  // if(ui)
  // Define (G)UI terminal for interactive mode  
  //{ 
    // pUI->ApplyCommand("/control/execute vis.mac");
    // ui->SessionStart();
    // delete ui;
  //}
  
  //else
  // Batch mode
  //{ 
  G4String command = "/control/execute vis.mac ";    
  pUI->ApplyCommand(command);
  //}

  delete visManager;
  delete runManager;

  return 0;
}



