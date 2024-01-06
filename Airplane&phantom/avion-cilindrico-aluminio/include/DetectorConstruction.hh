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

/// \file include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class


#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4MultiFunctionalDetector.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

//=======================================================================================================================
// User detector construction class
//   
//   [Geometry] 
//   The world volume is defined as 45 m x 45 m x 45 m box with Air.
//   The airplane is defined as a cylindrical section and two spherical crust sections. The cylindrical section has
//   radius of 7 and 6.996 m and a length of 70 m. The spherical sections represent the front and back of the plane,
//   they are two spherical crusts with radii 7 and 6.996 m and are cut in half. The material of the plane is Aluminium.
//   The interior of the plane is defined as the geometry of the plane but with a cylinder of radius 6.996 m and 
//   length  of 70 m. The spherical regions have a radius of 6.996 m with with Air.
//   Water phantom is defined as  1 m x 1 m x 1 m box with Water.
//   The water phantom is divided into 100 segments in x,y plane using replication, and 
//   then divided into 100 segments perpendicular to z axis using nested parameterised volume.  
//   These values are defined at constructor,
//   e.g. the size of water phantom (fPhantomSize), and number of segmentation of water phantom (fNx, fNy, fNz).
//
//   NIST database is used for materials.
//
//  [Scorer]
//  Assignment of G4MultiFunctionalDetector and G4PrimitiveScorer 
//  is demonstrated in this example.
//      -------------------------------------------------
//      The collection names of defined Primitives are
//       0       PhantomSD/totalEDep 
//       1       PhantomSD/protonEDep 
//       2       PhantomSD/neutronEDep
//       3       PhantomSD /electronEDep
//       4       PhantomSD/positronEDep 
//       5       PhantomSD/gammaEDep 
//       6       PhantomSD/pionpEDep
//       7       PhantomSD /pionnEDep
//       8       PhantomSD/muonpEDep
//       9       PhantomSD/muonnEDep
//     -------------------------------------------------
//     Please see README for detail description.
//
//
// - G4VPhysicalVolume* Construct()
//     retrieves material from NIST database,
//     constructs a water phantom "phantom", plane and air sphere in the world volume "world" and
//     sets detector sensitivities with G4MultiFunctionalDetector
//
// - void SetPhantomSize(G4ThreeVector size)
//     sets the water phantom size which is defined in G4Box
//
// - const G4ThreeVector& GetPhantomSize() const
//     gets the water phantom size
//
// - void SetNumberOfSegmentsInPhantom(G4int nx, G4int ny, G4int nz) 
//     sets the number of segments of the water phantom
//
// - void GetNumberOfSegmentsInPhantom(G4int& nx, G4int& ny, G4int& nz) 
//     gets the number of segments of the water phantom
//=======================================================================================================================


class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  // constructor and destructor.
  DetectorConstruction();
  virtual ~DetectorConstruction();

public:
  // virtual method from G4VUserDetectorConstruction.
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();

public:
  // Get/Set Access methods for data members
  // Size of Whater Phantom
  void SetPhantomSize(G4ThreeVector size) { fPhantomSize=size; }
  const G4ThreeVector& GetPhantomSize() const { return fPhantomSize; }
  // Number of segments of water phantom
  void SetNumberOfSegmentsInPhantom(G4int nx, G4int ny, G4int nz) 
      { fNx=nx; fNy=ny; fNz=nz; }
  void GetNumberOfSegmentsInPhantom(G4int& nx, G4int& ny, G4int& nz) 
     const{ nx=fNx; ny = fNy; nz = fNz; }

private:
  // Data members
  G4ThreeVector fPhantomSize;   // Size of Phantom
  G4int         fNx,fNy,fNz;    // Number of segmentation of phantom.
  G4LogicalVolume* fLVPhantomSens;
  G4LogicalVolume* logicAirplaneTube;
  G4LogicalVolume* logicAirplaneSphere1;
  G4LogicalVolume* logicAirplaneSphere2;
  G4LogicalVolume* logicAirTube;
  G4LogicalVolume* logicAirSphere1;
  G4LogicalVolume* logicAirSphere2;

};
#endif
