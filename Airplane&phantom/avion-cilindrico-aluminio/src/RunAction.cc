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

/// \file src/RunAction.cc
/// \brief Implementation of the RunAction class


#include "RunAction.hh"
#include "Run.hh"

//-- In order to obtain detector information.
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "G4THitsMap.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"    

#include <fstream>

//=======================================================================
// Constructor
RunAction::RunAction()
  : G4UserRunAction(),
    fNx(0), fNy(0), fNz(0)
{
  // - Prepare data member for Run.
  //   vector represents a list of MultiFunctionalDetector names.
  fSDName.push_back(G4String("PhantomSD"));
}

//=======================================================================
// Destructor.
RunAction::~RunAction()
{
  fSDName.clear();
}

//=======================================================================
G4Run* RunAction::GenerateRun()
{
  // Generate new RUN object, which is specially
  // dedicated for MultiFunctionalDetector scheme.
  //  Detail description can be found in Run.hh/cc.
  return new Run(fSDName);
}

//=======================================================================
void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

//=======================================================================
void RunAction::EndOfRunAction(const G4Run* aRun)
{
  if(!IsMaster()) return;

  //- Run object.
  Run* re02Run = (Run*)aRun;
  //--- Dump all socred quantities involved.
  re02Run->DumpAllScorer();
  //---

  //
  //- water phantom (Detector) Information.
  //-- Number of segments in the water phantom.
  const DetectorConstruction* detector =
      (const DetectorConstruction*)
      (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  detector->GetNumberOfSegmentsInPhantom(fNx,fNy,fNz); //Fill fNx,y,z.

  //---------------------------------------------
  // Dump accumulated quantities for this RUN.
  //  (Display only central region of x-y plane)
  //---------------------------------------------
  G4THitsMap<G4double>* totalEdep = re02Run->GetHitsMap("PhantomSD/totalEDep");
  G4THitsMap<G4double>* protonEdep = re02Run->GetHitsMap("PhantomSD/protonEDep");
  G4THitsMap<G4double>* neutronEdep = re02Run->GetHitsMap("PhantomSD/neutronEDep");
  G4THitsMap<G4double>* electronEdep = re02Run->GetHitsMap("PhantomSD/electronEDep");
  G4THitsMap<G4double>* positronEdep = re02Run->GetHitsMap("PhantomSD/positronEDep");
  G4THitsMap<G4double>* gammaEdep = re02Run->GetHitsMap("PhantomSD/gammaEDep");
  G4THitsMap<G4double>* muonpEdep = re02Run->GetHitsMap("PhantomSD/muonpEDep");
  G4THitsMap<G4double>* muonnEdep = re02Run->GetHitsMap("PhantomSD/muonnEDep");
  G4THitsMap<G4double>* pionpEdep = re02Run->GetHitsMap("PhantomSD/pionpEDep");
  G4THitsMap<G4double>* pionnEdep = re02Run->GetHitsMap("PhantomSD/pionnEDep");

  G4cout << "============================================================="
      << G4endl;
  G4cout << " Number of event processed : " << aRun->GetNumberOfEvent() << G4endl;
  G4cout << "============================================================="
      << G4endl;
  G4cout << std::setw(8) << "#Z Cell#";
  G4cout << std::setw(16) << totalEdep->GetName();
  G4cout << std::setw(16) << protonEdep->GetName();
  G4cout << std::setw(16) << neutronEdep->GetName();
  G4cout << std::setw(16) << electronEdep->GetName();
  G4cout << std::setw(16) << positronEdep->GetName();
  G4cout << std::setw(16) << gammaEdep->GetName();
  G4cout << std::setw(16) << muonpEdep->GetName();
  G4cout << std::setw(16) << muonnEdep->GetName();
  G4cout << std::setw(16) << pionpEdep->GetName();
  G4cout << std::setw(16) << pionnEdep->GetName()
      << G4endl;
  G4int ix = fNx / 2;
  G4int iy = fNy / 2;
  G4int iz;

  for (iz = 0; iz < fNz; iz++) {
      G4double* totalED = (*totalEdep)[CopyNo(ix, iy, iz)];
      G4double* protonED = (*protonEdep)[CopyNo(ix, iy, iz)];
      G4double* neutronED = (*neutronEdep)[CopyNo(ix, iy, iz)];
      G4double* electronED = (*electronEdep)[CopyNo(ix, iy, iz)];
      G4double* positronED = (*positronEdep)[CopyNo(ix, iy, iz)];
      G4double* gammaED = (*gammaEdep)[CopyNo(ix, iy, iz)];
      G4double* muonpED = (*muonpEdep)[CopyNo(ix, iy, iz)];
      G4double* muonnED = (*muonnEdep)[CopyNo(ix, iy, iz)];
      G4double* pionpED = (*pionpEdep)[CopyNo(ix, iy, iz)];
      G4double* pionnED = (*pionnEdep)[CopyNo(ix, iy, iz)];
      if (!totalED) totalED = new G4double(0.0);
      if (!protonED) protonED = new G4double(0.0);
      if (!neutronED) neutronED = new G4double(0.0);
      if (!electronED) electronED = new G4double(0.0);
      if (!positronED) positronED = new G4double(0.0);
      if (!gammaED) gammaED = new G4double(0.0);
      if (!muonpED) muonpED = new G4double(0.0);
      if (!muonnED) muonnED = new G4double(0.0);
      if (!pionpED) pionpED = new G4double(0.0);
      if (!pionnED) pionnED = new G4double(0.0);
      G4cout << std::setw(6) << iz << "  "
          << std::setw(12) << G4BestUnit(*totalED, "Energy")
          << std::setw(12) << G4BestUnit(*protonED, "Energy")
          << std::setw(12) << G4BestUnit(*neutronED, "Energy")
          << std::setw(12) << G4BestUnit(*electronED, "Energy")
          << std::setw(12) << G4BestUnit(*positronED, "Energy")
          << std::setw(12) << G4BestUnit(*gammaED, "Energy")
          << std::setw(12) << G4BestUnit(*muonpED, "Energy")
          << std::setw(12) << G4BestUnit(*muonnED, "Energy")
          << std::setw(12) << G4BestUnit(*pionpED, "Energy")
          << std::setw(12) << G4BestUnit(*pionnED, "Energy")
          << G4endl;
  }
  G4cout << "=============================================" << G4endl;

  std::ofstream  file("totalED.txt");
  for (iz = 0; iz < fNz; iz++) {
      for (iy = 0; iy < fNy; iy++) {
          for (ix = 0; ix < fNx; ix++) {
              G4double* totalED = (*totalEdep)[CopyNo(ix, iy, iz)];
              if (!totalED) totalED = new G4double(0.0);
              file << ix << " " << iy << " " << iz << " " << *totalED / MeV << G4endl;

          }
      }
  }
  file.close();

  std::ofstream  file1("protonED.txt");
  for (iz = 0; iz < fNz; iz++) {
      for (iy = 0; iy < fNy; iy++) {
          for (ix = 0; ix < fNx; ix++) {
              G4double* protonED = (*protonEdep)[CopyNo(ix, iy, iz)];
              if (!protonED) protonED = new G4double(0.0);
              file1 << ix << " " << iy << " " << iz << " " << *protonED / MeV << G4endl;

          }
      }
  }
  file1.close();

  std::ofstream  file2("neutronED.txt");
  for (iz = 0; iz < fNz; iz++) {
      for (iy = 0; iy < fNy; iy++) {
          for (ix = 0; ix < fNx; ix++) {
              G4double* neutronED = (*neutronEdep)[CopyNo(ix, iy, iz)];
              if (!neutronED) neutronED = new G4double(0.0);
              file2 << ix << " " << iy << " " << iz << " " << *neutronED / MeV << G4endl;

          }
      }
  }
  file2.close();

  std::ofstream  file3("electronED.txt");
  for (iz = 0; iz < fNz; iz++) {
      for (iy = 0; iy < fNy; iy++) {
          for (ix = 0; ix < fNx; ix++) {
              G4double* electronED = (*electronEdep)[CopyNo(ix, iy, iz)];
              if (!electronED) electronED = new G4double(0.0);
              file3 << ix << " " << iy << " " << iz << " " << *electronED / MeV << G4endl;

          }
      }
  }
  file3.close();

  std::ofstream  file4("positronED.txt");
  for (iz = 0; iz < fNz; iz++) {
      for (iy = 0; iy < fNy; iy++) {
          for (ix = 0; ix < fNx; ix++) {
              G4double* positronED = (*positronEdep)[CopyNo(ix, iy, iz)];
              if (!positronED) positronED = new G4double(0.0);
              file4 << ix << " " << iy << " " << iz << " " << *positronED / MeV << G4endl;

          }
      }
  }
  file4.close();

  std::ofstream  file5("gammaED.txt");
  for (iz = 0; iz < fNz; iz++) {
      for (iy = 0; iy < fNy; iy++) {
          for (ix = 0; ix < fNx; ix++) {
              G4double* gammaED = (*gammaEdep)[CopyNo(ix, iy, iz)];
              if (!gammaED) gammaED = new G4double(0.0);
              file5 << ix << " " << iy << " " << iz << " " << *gammaED / MeV << G4endl;

          }
      }
  }
  file5.close();

  std::ofstream  file6("muonpED.txt");
  for (iz = 0; iz < fNz; iz++) {
      for (iy = 0; iy < fNy; iy++) {
          for (ix = 0; ix < fNx; ix++) {
              G4double* muonpED = (*muonpEdep)[CopyNo(ix, iy, iz)];
              if (!muonpED) muonpED = new G4double(0.0);
              file6 << ix << " " << iy << " " << iz << " " << *muonpED / MeV << G4endl;

          }
      }
  }
  file6.close();

  std::ofstream  file7("muonnED.txt");
  for (iz = 0; iz < fNz; iz++) {
      for (iy = 0; iy < fNy; iy++) {
          for (ix = 0; ix < fNx; ix++) {
              G4double* muonnED = (*muonnEdep)[CopyNo(ix, iy, iz)];
              if (!muonnED) muonnED = new G4double(0.0);
              file7 << ix << " " << iy << " " << iz << " " << *muonnED / MeV << G4endl;

          }
      }
  }
  file7.close();

  std::ofstream  file8("pionpED.txt");
  for (iz = 0; iz < fNz; iz++) {
      for (iy = 0; iy < fNy; iy++) {
          for (ix = 0; ix < fNx; ix++) {
              G4double* pionpED = (*pionpEdep)[CopyNo(ix, iy, iz)];
              if (!pionpED) pionpED = new G4double(0.0);
              file8 << ix << " " << iy << " " << iz << " " << *pionpED / MeV << G4endl;

          }
      }
  }
  file8.close();

  std::ofstream  file9("pionnED.txt");
  for (iz = 0; iz < fNz; iz++) {
      for (iy = 0; iy < fNy; iy++) {
          for (ix = 0; ix < fNx; ix++) {
              G4double* pionnED = (*pionnEdep)[CopyNo(ix, iy, iz)];
              if (!pionnED) pionnED = new G4double(0.0);
              file9 << ix << " " << iy << " " << iz << " " << *pionnED / MeV << G4endl;

          }
      }
  }
  file9.close();



}

