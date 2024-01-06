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

/// \file src/PSEnergyDeposit.cc
/// \brief Implementation of the PSEnergyDeposit class


#include "PSEnergyDeposit.hh"

//==========================================================================
// (Description)
//   This is a primitive scorer class for scoring cell charge.
//   The Cell Charge is defined by  a sum of charge inside the cell
//   which calculates the deposited charge in the cell.  
//==========================================================================

//==========================================================================
PSEnergyDeposit::PSEnergyDeposit(G4String name,
                                         G4int nx, G4int ny, G4int nz)
  :G4PSEnergyDeposit(name),fNx(nx),fNy(ny),fNz(nz)
{;}

//==========================================================================
PSEnergyDeposit::~PSEnergyDeposit()
{;}

//==========================================================================
G4int PSEnergyDeposit::GetIndex(G4Step* aStep)
{
  const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int ix = touchable->GetReplicaNumber(1);
  G4int iy = touchable->GetReplicaNumber(2);
  G4int iz = touchable->GetReplicaNumber(0);

  G4int tmp = fNy;
  if (tmp) return iy*fNx*fNz+ix*fNz+iz;
  else return iy*fNx*fNz+ix*fNz+iz;
}
