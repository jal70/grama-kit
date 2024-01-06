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
//
/// \file src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

 
#include "DetectorConstruction.hh"

#include "G4PSEnergyDeposit3D.hh"

#include "G4SDParticleWithEnergyFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4SDChargedFilter.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"

#include "G4PVParameterised.hh"
#include "NestedPhantomParameterisation.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"    
#include "G4ios.hh"

//====================================================================================================================
//  DetectorConstruction
//
//    
//   [Geometry] 
//   The world volume is defined as 30 m x 30 m x 30 m box with Air.
//   The plane is defined as a sphericall shell of radius 7 and 6.99855 m with GLARE.
//   Glare is an optimised Fibre Metal Laminate for aircraft and it consists of alternating layers of aluminium and 
//   glass fibre pregreg layers. A Glare laminate with fibre orientation according to the Glare 3A-3 definition:
//   - 3 layers of aluminium and 2 fibre layers.
//   - An aluminium layer thickness of 0.4 mm.
//   - An fibre  layer thickness of 0.125mm.
//   More information about Glare: http://www.fmlc.nl/research-development/results-cases/
//   The interior of the plane is defined as a sphere of radius 6.99855 m with Air.
//   Water phantom is defined as  1 m x 1 m x 1 m box with Water.
//   The water phantom is divided into 100 segments in x,y plane using replication, and 
//   then divided into 100 segments perpendicular to z axis using nested parameterised volume.  
//   These values are defined at constructor,
//   e.g. the size of water phantom (fPhantomSize), and number of segmentation of water phantom (fNx, fNy, fNz).
//
//   NIST database is used for materials.
//
//
//   [Scorer]
//   Assignment of G4MultiFunctionalDetector and G4PrimitiveScorer 
//   is demonstrated in this example.
//       -------------------------------------------------
//       The collection names of defined Primitives are
//        0       PhantomSD/totalEDep 
//        1       PhantomSD/protonEDep
//        2       PhantomSD/neutronEDep
//        3       PhantomSD/electonEDep
//        4       PhantomSD/positronEDep
//        5       PhantomSD/gammaEDep
//        6       PhantomSD/muonpEDep
//        7       PhantomSD/muonnEDep
//        8       PhantomSD/pionpEDep
//        9       PhantomSD/pionnEDep       
//      -------------------------------------------------
//      Please see README for detail description.
//
//====================================================================================================================

//====================================================================================================================
DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction() 
{
  // Default size of water phantom,and segmentation.
    fPhantomSize.setX(1.*m);
    fPhantomSize.setY(1.*m);
    fPhantomSize.setZ(1.*m);
    fNx = fNy = fNz = 100;
}

//====================================================================================================================
DetectorConstruction::~DetectorConstruction()
{;}

//====================================================================================================================
G4VPhysicalVolume* DetectorConstruction::Construct()
{

  //==================================
  //       Material Definitions       
  //==================================
  
  //----------------------------------------------- NIST Materials --------------------------------------------------
  //  Material Information imported from NIST database.
  
  G4NistManager* NISTman = G4NistManager::Instance();
  
  G4Material* air  = NISTman->FindOrBuildMaterial("G4_AIR");
  G4Material* water  = NISTman->FindOrBuildMaterial("G4_WATER");
  
  // -- Elements --
  G4Material* Al = NISTman->FindOrBuildMaterial("G4_Al");
  G4Material* Cu = NISTman->FindOrBuildMaterial("G4_Cu");
  G4Material* Mg = NISTman->FindOrBuildMaterial("G4_Mg");
  G4Material* Mn = NISTman->FindOrBuildMaterial("G4_Mn");
  G4Material* Zn = NISTman->FindOrBuildMaterial("G4_Zn");
  G4Material* He = NISTman->FindOrBuildMaterial("G4_He");
  G4Material* Ti = NISTman->FindOrBuildMaterial("G4_Ti");
  G4Material* Cr = NISTman->FindOrBuildMaterial("G4_Cr");
  G4Material* Si = NISTman->FindOrBuildMaterial("G4_Si");
  G4Material* O = NISTman->FindOrBuildMaterial("G4_O");

  // -------- Define a mixture by fractional mass -------
  G4double fractionmass, density;
  G4String name, symbol;
  G4int ncomponents;

  // -- Material: Aluminium2024 -- 
  density = 2.77 * g / cm3;
  G4Material* Aluminium2024 = new G4Material(name = "Aluminium2024", density, ncomponents = 9);
  Aluminium2024->AddMaterial(Al, fractionmass = 91.55 * perCent);
  Aluminium2024->AddMaterial(Cu, fractionmass = 4.3 * perCent);
  Aluminium2024->AddMaterial(Mg, fractionmass = 1.5 * perCent);
  Aluminium2024->AddMaterial(Mn, fractionmass = 0.5 * perCent);
  Aluminium2024->AddMaterial(Zn, fractionmass = 0.5 * perCent);
  Aluminium2024->AddMaterial(He, fractionmass = 0.5 * perCent);
  Aluminium2024->AddMaterial(Ti, fractionmass = 0.15 * perCent);
  Aluminium2024->AddMaterial(Si, fractionmass = 0.5 * perCent);
  Aluminium2024->AddMaterial(Cr, fractionmass = 0.5 * perCent);

  // -- Material: GLARE --
  density = 2.46 * g / cm3;
  G4Material* S2Glass  = new G4Material(name = "S2Glass", density, ncomponents = 4);
  S2Glass->AddMaterial(Si, fractionmass = 30.37 * perCent);
  S2Glass->AddMaterial(O, fractionmass = 50.36 * perCent);
  S2Glass->AddMaterial(Al, fractionmass = 13.23 * perCent);
  S2Glass->AddMaterial(Mg, fractionmass = 6.04 * perCent);

  // Print all the materials defined.
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;


  //============================================================================
  //         Definitions of Solids, Logical Volumes, Physical Volumes            
  //============================================================================


  //=====================
  //     World Volume 
  //=====================

  G4double worldSize = 30.0 * m;

  G4Box * solidWorld = 
      new G4Box("world", worldSize/2., worldSize/2., worldSize/2.);

  G4LogicalVolume * logicWorld = 
      new G4LogicalVolume(solidWorld, air, "World", 0, 0, 0);

  
  //  Must place the World Physical volume unrotated at (0,0,0).
  G4VPhysicalVolume * physiWorld
    = new G4PVPlacement(0,               // no rotation
                        G4ThreeVector(), // at (0,0,0)
                        logicWorld,      // its logical volume
                        "World",         // its name
                        0,               // its mother  volume
                        false,           // no boolean operations
                        0);              // copy number
     


  //==================
  //     Airplane
  //==================


  G4Sphere* solidAirplane =
      new G4Sphere("AirplaneSolid",
          6.9996 * m, 
          7.0 * m,
          0.0 * deg,
          360.0 * deg,
          0.0 * deg,
          360.0 * deg);


  logicAirplane =
      new G4LogicalVolume(solidAirplane, Aluminium2024, "AirplaneLogical");

  G4ThreeVector positionAirplane = G4ThreeVector(0, 0, 0);

  new G4PVPlacement(0,
      positionAirplane,
      logicAirplane,
      "Airplane",
      logicWorld,
      false,
      0);

  // -- 
  G4Sphere* solidAirplane1 = 
      new G4Sphere("AirplaneSolid1", 
          6.999475 * m,
          6.9996 * m, 
          0.0 * deg,
          360.0 * deg,
          0.0 * deg,
          360.0 * deg);


  logicAirplane1 =
      new G4LogicalVolume(solidAirplane1, S2Glass, "AirplaneLogical1");

  G4ThreeVector positionAirplane1 = G4ThreeVector(0, 0, 0);

  new G4PVPlacement(0,
      positionAirplane1,
      logicAirplane1,
      "Airplane1",
      logicWorld,
      false,
      0);

  // -- 
  G4Sphere* solidAirplane2 =
      new G4Sphere("AirplaneSolid2",
          6.999075 * m,
          6.999475 * m,
          0.0 * deg,
          360.0 * deg,
          0.0 * deg,
          360.0 * deg);


  logicAirplane2 =
      new G4LogicalVolume(solidAirplane2, Aluminium2024, "AirplaneLogical2");

  G4ThreeVector positionAirplane2 = G4ThreeVector(0, 0, 0);

  new G4PVPlacement(0,
      positionAirplane2,
      logicAirplane2,
      "Airplane2",
      logicWorld,
      false,
      0);

  // -- 
  G4Sphere* solidAirplane3 =
      new G4Sphere("AirplaneSolid3",
          6.99895 * m,
          6.999075 * m,
          0.0 * deg,
          360.0 * deg,
          0.0 * deg,
          360.0 * deg);


  logicAirplane3 =
      new G4LogicalVolume(solidAirplane3, S2Glass, "AirplaneLogical3");

  G4ThreeVector positionAirplane3 = G4ThreeVector(0, 0, 0);

  new G4PVPlacement(0,
      positionAirplane3,
      logicAirplane3,
      "Airplane3",
      logicWorld,
      false,
      0);

  // -- 
  G4Sphere* solidAirplane4 =
      new G4Sphere("AirplaneSolid4",
          6.99855 * m,
          6.99895 * m,
          0.0 * deg,
          360.0 * deg,
          0.0 * deg,
          360.0 * deg);


  logicAirplane4 =
      new G4LogicalVolume(solidAirplane4, Aluminium2024, "AirplaneLogical4");

  G4ThreeVector positionAirplane4 = G4ThreeVector(0, 0, 0);

  new G4PVPlacement(0,
      positionAirplane4,
      logicAirplane4,
      "Airplane4",
      logicWorld,
      false,
      0);



  //====================
  //     Air Sphere
  //====================


  G4Sphere* solidAirSphere =
      new G4Sphere("AirSolid",
          0.0 * m,
          6.99855 * m,
          0.0 * deg,
          360.0 * deg,
          0.0 * deg,
          360.0 * deg);

  logicAirSphere =
      new G4LogicalVolume(solidAirSphere, air, "AirLogical");

  G4ThreeVector positionAirSphere = G4ThreeVector(0, 0, 0);

  G4VPhysicalVolume* physAir =
      new G4PVPlacement(0,
          positionAirSphere,
          logicAirSphere,
          "air",
          logicWorld,
          false,
          0);


  //=================
  //     Phantom
  //=================

  //-------------------------
  // Mother Volume of Phantom
  //-------------------------

  

  //--  Default size of water phantom is defined at constructor.
  G4ThreeVector phantomSize = fPhantomSize; 
  
  G4Box * solidPhantom
    = new G4Box("phantom",
                phantomSize.x()/2., phantomSize.y()/2., phantomSize.z()/2.);

  G4LogicalVolume * logicPhantom
    = new G4LogicalVolume(solidPhantom, water, "Phantom", 0, 0, 0);  

  G4ThreeVector positionPhantom = G4ThreeVector(0, 0, 0);


  //G4VPhysicalVolume * physiPhantom =
  new G4PVPlacement(0,             
                    positionPhantom, 
                    logicPhantom,    
                    "Phantom",       
                    logicAirSphere,      
                    false,           
                    0);             

  
  //=====================================================================
  //                            Visualiation
  //=====================================================================

  G4VisAttributes* boxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  G4VisAttributes* airplaneVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  G4VisAttributes* airVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
  G4VisAttributes* phantomVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));

  airplaneVisAtt->SetForceWireframe(true);
  airVisAtt->SetForceWireframe(true);
  phantomVisAtt->SetForceWireframe(true);

  logicWorld->SetVisAttributes(boxVisAtt);
  logicAirplane->SetVisAttributes(airplaneVisAtt);
  logicAirSphere->SetVisAttributes(airVisAtt);
  logicPhantom->SetVisAttributes(phantomVisAtt);
  
                                     
  //=====================================================================
  //             Phantom segmentation using Parameterisation         
  //=====================================================================
  
  G4cout << "<-- DetectorConstruction::Construct-------" <<G4endl;
  G4cout << "  Water Phantom Size " << fPhantomSize/mm       << G4endl;
  G4cout << "  Segmentation  ("<< fNx<<","<<fNy<<","<<fNz<<")"<< G4endl;
  G4cout << "<---------------------------------------------"<< G4endl;
  
  //---------------------- Number of segmentation ----------------------
  // - Default number of segmentation is defined at constructor.  

  G4int nxCells = fNx;
  G4int nyCells = fNy;
  G4int nzCells = fNz;

  G4ThreeVector sensSize;
  sensSize.setX(phantomSize.x()/(G4double)nxCells);
  sensSize.setY(phantomSize.y()/(G4double)nyCells);
  sensSize.setZ(phantomSize.z()/(G4double)nzCells);
  
   
  //====================================================================
  //                    Replication of Phantom Volume
  //====================================================================

  //---------------------------- Y Slice -------------------------------

  G4String yRepName("RepY");
  
  G4VSolid* solYRep =
    new G4Box(yRepName,phantomSize.x()/2.,sensSize.y()/2.,phantomSize.z()/2.);
  
  G4LogicalVolume* logYRep =
    new G4LogicalVolume(solYRep,water,yRepName);
  
  //G4PVReplica* yReplica =
  new G4PVReplica(yRepName,logYRep,logicPhantom,kYAxis,fNy,sensSize.y());
  
  
  //--------------------------- X Slice --------------------------------
  G4String xRepName("RepX");
  
  G4VSolid* solXRep =
    new G4Box(xRepName,sensSize.x()/2.,sensSize.y()/2.,phantomSize.z()/2.);
  
  G4LogicalVolume* logXRep =
    new G4LogicalVolume(solXRep,water,xRepName);
  
  //G4PVReplica* xReplica =
  new G4PVReplica(xRepName,logXRep,logYRep,kXAxis,fNx,sensSize.x());

  
  //====================================================================
  //                  Voxel solid and logical volumes
  //====================================================================

  //--------------------------- Z Slice --------------------------------

  G4String zVoxName("phantomSens");
  
  G4VSolid* solVoxel = 
    new G4Box(zVoxName,sensSize.x()/2.,sensSize.y()/2.,sensSize.z()/2.);
  
  fLVPhantomSens = new G4LogicalVolume(solVoxel,water,zVoxName);
  
  
  std::vector<G4Material*> phantomMat(2,water);
  
  //- Parameterisation for transformation of voxels.
  // (voxel size is fixed. 
  // e.g. nested parameterisation handles material and transfomation of voxels.)

  NestedPhantomParameterisation* paramPhantom
    = new NestedPhantomParameterisation(sensSize/2.,nzCells,phantomMat);

  //G4VPhysicalVolume * physiPhantomSens =
    new G4PVParameterised("PhantomSens",     // their name
                          fLVPhantomSens,    // their logical volume
                          logXRep,           // Mother logical volume
                          kUndefined,        // Are placed along this axis 
                          nzCells,           // Number of cells
                          paramPhantom);     // Parameterisation.

  //   Optimization flag is avaiable for,
  //   kUndefined, kXAxis, kYAxis, kZAxis.
 

  // Replica
  G4VisAttributes* yRepVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  logYRep->SetVisAttributes(yRepVisAtt);

  G4VisAttributes* xRepVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  logXRep->SetVisAttributes(xRepVisAtt);
  
  // Skip the visualization for those voxels.
  fLVPhantomSens->SetVisAttributes(G4VisAttributes::GetInvisible());

  
  return physiWorld;
}

void DetectorConstruction::ConstructSDandField() {

  //====================================================================
  //         Sensitive detectors : MultiFunctionalDetector
  //====================================================================
  
  // Sensitive Detector Manager.
  G4SDManager* pSDman = G4SDManager::GetSDMpointer();
  
  // Sensitive Detector Name
  G4String phantomSDname = "PhantomSD";
  G4String airplaneSDname = "AirplaneSD";
  
  //--------------------------------------------------------------------
  //                      MultiFunctionalDetector
  //--------------------------------------------------------------------
  
  // Define MultiFunctionalDetector with name.
  G4MultiFunctionalDetector* mFDet
  = new G4MultiFunctionalDetector(phantomSDname);
  pSDman->AddNewDetector( mFDet );                // Register SD to SDManager.
  fLVPhantomSens->SetSensitiveDetector(mFDet);    // Assign SD to the logical volume.

  G4MultiFunctionalDetector* mFDet1
      = new G4MultiFunctionalDetector(airplaneSDname);
  pSDman->AddNewDetector(mFDet1);                       // Register SD to SDManager.
  logicAirplane->SetSensitiveDetector(mFDet1);          // Assign SD to the logical volume.
  logicAirplane1->SetSensitiveDetector(mFDet1);         // Assign SD to the logical volume.
  logicAirplane2->SetSensitiveDetector(mFDet1);         // Assign SD to the logical volume.
  logicAirplane3->SetSensitiveDetector(mFDet1);         // Assign SD to the logical volume.
  logicAirplane4->SetSensitiveDetector(mFDet1);         // Assign SD to the logical volume.
  logicAirSphere->SetSensitiveDetector(mFDet1);         // Assign SD to the logical volume.

  
  //--------------------------------------------------------------------
  //              SDFilter : Sensitive Detector Filters
  //--------------------------------------------------------------------
  
  // Particle Filter for Primitive Scorer with filter name(fltName)
  // and particle name(particleName),
  // or particle names are given by add("particle name"); method.
  
  G4String fltName,particleName;
  
  //-- proton filter
  G4SDParticleFilter* protonFilter =
  new G4SDParticleFilter(fltName="protonFilter", particleName="proton");

  
  //-- neutron filter
  G4SDParticleFilter* neutronFilter =
  new G4SDParticleFilter(fltName = "neutronFilter", particleName = "neutron");

  
  //-- electron filter
  G4SDParticleFilter* electronFilter =
      new G4SDParticleFilter(fltName = "electronFilter", particleName = "e-");


  //-- positron filter
  G4SDParticleFilter* positronFilter =
      new G4SDParticleFilter(fltName = "positronFilter", particleName = "e+");


  //-- gamma filter
  G4SDParticleFilter* gammaFilter =
      new G4SDParticleFilter(fltName = "gammaFilter", particleName = "gamma");


  //-- mu+ filter
  G4SDParticleFilter* muonpFilter =
      new G4SDParticleFilter(fltName = "muonpFilter", particleName = "mu+");


  //-- mu- filter
  G4SDParticleFilter* muonnFilter =
      new G4SDParticleFilter(fltName = "muonnFilter", particleName = "mu-");


  //-- pi+ filter
  G4SDParticleFilter* pionpFilter =
      new G4SDParticleFilter(fltName = "pionpFilter", particleName = "pi+");


  //-- pi- filter
  G4SDParticleFilter* pionnFilter =
      new G4SDParticleFilter(fltName = "pionnFilter", particleName = "pi-");


  //====================================================================
  //                     PS : Primitive Scorers
  //====================================================================

  // Primitive Scorers are used with SDFilters according to your purpose.
  
  //-- Primitive Scorer for Energy Deposit.
 
  G4String psName;
  G4PSEnergyDeposit3D * scorer0 = new G4PSEnergyDeposit3D(psName="totalEDep",
                                                          fNx,fNy,fNz);

  G4PSEnergyDeposit3D * scorer1 = new G4PSEnergyDeposit3D(psName="protonEDep",
                                                          fNx,fNy,fNz);
  scorer1->SetFilter(protonFilter);

  G4PSEnergyDeposit3D * scorer2 = new G4PSEnergyDeposit3D(psName = "neutronEDep",
                                                         fNx, fNy, fNz);
  scorer2->SetFilter(neutronFilter);

  G4PSEnergyDeposit3D* scorer3 = new G4PSEnergyDeposit3D(psName = "electronEDep",
      fNx, fNy, fNz);
  scorer3->SetFilter(electronFilter);

  G4PSEnergyDeposit3D* scorer4 = new G4PSEnergyDeposit3D(psName = "positronEDep",
      fNx, fNy, fNz);
  scorer4->SetFilter(positronFilter);

  G4PSEnergyDeposit3D* scorer5 = new G4PSEnergyDeposit3D(psName = "gammaEDep",
      fNx, fNy, fNz);
  scorer5->SetFilter(gammaFilter);

  G4PSEnergyDeposit3D* scorer6 = new G4PSEnergyDeposit3D(psName = "muonpEDep",
      fNx, fNy, fNz);
  scorer6->SetFilter(muonpFilter);

  G4PSEnergyDeposit3D* scorer7 = new G4PSEnergyDeposit3D(psName = "muonnEDep",
      fNx, fNy, fNz);
  scorer7->SetFilter(muonnFilter);

  G4PSEnergyDeposit3D* scorer8 = new G4PSEnergyDeposit3D(psName = "pionpEDep",
      fNx, fNy, fNz);
  scorer8->SetFilter(pionpFilter);

  G4PSEnergyDeposit3D* scorer9 = new G4PSEnergyDeposit3D(psName = "pionnEDep",
      fNx, fNy, fNz);
  scorer9->SetFilter(pionnFilter);
  
  //--------------------------------------------------------------------
  //      Register primitive scorers to MultiFunctionalDetector
  //--------------------------------------------------------------------

  mFDet->RegisterPrimitive(scorer0);
  mFDet->RegisterPrimitive(scorer1);
  mFDet->RegisterPrimitive(scorer2);
  mFDet->RegisterPrimitive(scorer3);
  mFDet->RegisterPrimitive(scorer4);
  mFDet->RegisterPrimitive(scorer5);
  mFDet->RegisterPrimitive(scorer6);
  mFDet->RegisterPrimitive(scorer7);
  mFDet->RegisterPrimitive(scorer8);
  mFDet->RegisterPrimitive(scorer9);
 }

