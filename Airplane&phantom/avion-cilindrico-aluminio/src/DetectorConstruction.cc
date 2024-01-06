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
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"

#include "G4PVParameterised.hh"
#include "NestedPhantomParameterisation.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"    
#include "G4ios.hh"

//=========================================================================================================================
//  DetectorConstruction
//
//    
//   [Geometry] 
//   The world volume is defined as 45 m x 45 m x 45 m box with Air.
//   The airplane is defined as a cylindrical section and two spherical crust sections. The cylindrical section has radius
//   of 7 and 6.996 m and a length of 70 m. The spherical sections represent the front and back of the plane, they are 
//   two spherical crusts with radii 7 and 6.996 m and are cut in half. The material of the plane is Aluminium.
//   The interior of the plane is defined as the geometry of the plane but with a cylinder of radius 6.996 m and length  of 
//   70 m. The spherical regions have a radius of 6.996 m with with Air.
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
//       
//      -------------------------------------------------
//      Please see README for detail description.
//
//=========================================================================================================================

//=========================================================================================================================
DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction() 
{
  // Default size of water phantom,and segmentation.
    fPhantomSize.setX(1.*m);
    fPhantomSize.setY(1.*m);
    fPhantomSize.setZ(1.*m);
    fNx = fNy = fNz = 100;
}

//=========================================================================================================================
DetectorConstruction::~DetectorConstruction()
{;}

//=========================================================================================================================
G4VPhysicalVolume* DetectorConstruction::Construct()
{

    //==================================
    //       Material Definitions       
    //==================================

    //----------------------------------------------- NIST Materials --------------------------------------------------
    //  Material Information imported from NIST database.

    G4NistManager* NISTman = G4NistManager::Instance();

    // -- Material: Air --
    G4Material* air = NISTman->FindOrBuildMaterial("G4_AIR");
    // -- Material: Water --
    G4Material* water = NISTman->FindOrBuildMaterial("G4_WATER");

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


    // -- Print all the materials defined --
    G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;


    //============================================================================
    //         Definitions of Solids, Logical Volumes, Physical Volumes            
    //============================================================================


    //=====================
    //     World Volume 
    //=====================

    //G4ThreeVector worldSize = G4ThreeVector(6000*cm, 6000*cm, 6000*cm);

    G4double worldSize = 90.0 * m;

    G4Box* solidWorld =
        new G4Box("world", worldSize / 2., worldSize / 2., worldSize / 2.);

    G4LogicalVolume* logicWorld =
        new G4LogicalVolume(solidWorld, air, "World", 0, 0, 0);


    //  Must place the World Physical volume unrotated at (0,0,0).
    G4VPhysicalVolume* physiWorld
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

    // -- Tube --
    G4Tubs* solidAirplaneTube =
        new G4Tubs("AirplaneSolidTube", // name
            6.996 * m,                // inner radius
            7.0 * m,                // outer radius
            35.0 * m,               // Half length i
            0.0 * deg,              // Starting phi angle
            360.0 * deg);           // Angle of the segment


    logicAirplaneTube =
        new G4LogicalVolume(solidAirplaneTube, Aluminium2024, "AirplaneLogicalTube");

    G4ThreeVector positionAirplaneTube = G4ThreeVector(0, 0, 0);

    G4RotationMatrix* rotAirplaneTube = new G4RotationMatrix();
    rotAirplaneTube->rotateX(180. * deg);

    new G4PVPlacement(rotAirplaneTube,
        positionAirplaneTube,
        logicAirplaneTube,
        "AirplaneTube",
        logicWorld,
        false,
        0);

    // -- Sphere --
    // Curvature right side
    G4Sphere* solidAirplaneSphere1 =
        new G4Sphere("AirplaneSolidSphere1", // name
            6.996 * m,                       // inner radius
            7.0 * m,                         // outer radius
            90.0 * deg,                      // starting Phi angle of the segment
            270.0 * deg,                     // delta Phi angle of the segment
            90.0 * deg,                      // starting Theta angle of the segment
            270.0 * deg);                    // delta Theta of the segment


    logicAirplaneSphere1 =
        new G4LogicalVolume(solidAirplaneSphere1, Aluminium2024, "AirplaneLogicalSphere1");

    G4ThreeVector positionAirplaneSphere1 = G4ThreeVector(0, 0, -35.0 * m);

    new G4PVPlacement(0,
        positionAirplaneSphere1,
        logicAirplaneSphere1,
        "AirplaneSphere1",
        logicWorld,
        false,
        0);

    // -- Sphere --
    // Curvature left side
    G4Sphere* solidAirplaneSphere2 =
        new G4Sphere("AirplaneSolidSphere2", // name
            6.996 * m,                     // inner radius
            7.0 * m,                       // outer radius
            90.0 * deg,                     // starting Phi angle of the segment
            270.0 * deg,                   // delta Phi angle of the segment
            90.0 * deg,                     // starting Theta angle of the segment
            270.0 * deg);                  // delta Theta of the segment


    logicAirplaneSphere2 =
        new G4LogicalVolume(solidAirplaneSphere2, Aluminium2024, "AirplaneLogicalSphere2");

    G4ThreeVector positionAirplaneSphere2 = G4ThreeVector(0, 0, 35.0 * m);

    G4RotationMatrix* rotAirplaneSphere2 = new G4RotationMatrix();
    rotAirplaneSphere2->rotateX(180. * deg);

    new G4PVPlacement(rotAirplaneSphere2,
        positionAirplaneSphere2,
        logicAirplaneSphere2,
        "AirplaneSphere2",
        logicWorld,
        false,
        0);


    //=================
    //       Air 
    //=================

    // -- Tube --
    G4Tubs* solidAirTube =
        new G4Tubs("AirSolidTube",
            0.0 * m,
            6.996 * m,
            35.0 * m,
            0.0 * deg,
            360.0 * deg);

    logicAirTube =
        new G4LogicalVolume(solidAirTube, air, "AirLogicalTube");

    G4ThreeVector positionAirTube = G4ThreeVector(0, 0, 0);

    G4RotationMatrix* rotAirTube = new G4RotationMatrix();
    rotAirTube->rotateX(180. * deg);

    G4VPhysicalVolume* physAirTube =
        new G4PVPlacement(rotAirTube,
            positionAirTube,
            logicAirTube,
            "airTube",
            logicWorld,
            false,
            0);
  
    // -- Sphere --
    // Curvature right side
    G4Sphere* solidAirSphere1 =
        new G4Sphere("AirSolidSphere1",     // name
            0 * m,                          // inner radius
            6.996 * m,                      // outer radius
            90.0 * deg,                     // starting Phi angle of the segment
            270.0 * deg,                    // delta Phi angle of the segment
            90.0 * deg,                     // starting Theta angle of the segment
            270.0 * deg);                   // delta Theta of the segment


    logicAirSphere1 =
        new G4LogicalVolume(solidAirSphere1, air, "AirLogicalSphere1");

    G4ThreeVector positionAirSphere1 = G4ThreeVector(0, 0, 0);

    new G4PVPlacement(0,
        positionAirSphere1,
        logicAirSphere1,
        "AirSphere1",
        logicAirplaneSphere1,
        false,
        0);

    // -- Sphere --
    // Curvature left side
    G4Sphere* solidAirSphere2 =
        new G4Sphere("AirSolidSphere2", // name
            0 * m,                      // inner radius
            6.996 * m,                  // outer radius
            90.0 * deg,                 // starting Phi angle of the segment
            270.0 * deg,                // delta Phi angle of the segment
            90.0 * deg,                 // starting Theta angle of the segment
            270.0 * deg);               // delta Theta of the segment


    logicAirSphere2 =
        new G4LogicalVolume(solidAirSphere2, air, "AirLogicalSphere2");

    G4ThreeVector positionAirSphere2 = G4ThreeVector(0, 0, 0);

    new G4PVPlacement(0,
        positionAirSphere2,
        logicAirSphere2,
        "AirSphere2",
        logicAirplaneSphere2,
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
                    logicAirTube,      
                    false,           
                    0);             

  
    //=====================================================================
    //                            Visualiation
    //=====================================================================

    G4VisAttributes* boxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    G4VisAttributes* airplaneTubeVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    G4VisAttributes* airplaneSphere1VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    G4VisAttributes* airplaneSphere2VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    G4VisAttributes* airTubeVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
    G4VisAttributes* airSphere1VisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
    G4VisAttributes* airSphere2VisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
    G4VisAttributes* phantomVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));

    airplaneTubeVisAtt->SetForceWireframe(true);
    airplaneSphere1VisAtt->SetForceWireframe(true);
    airplaneSphere2VisAtt->SetForceWireframe(true);
    airTubeVisAtt->SetForceWireframe(true);
    airSphere1VisAtt->SetForceWireframe(true);
    airSphere2VisAtt->SetForceWireframe(true);
    phantomVisAtt->SetForceWireframe(true);

    logicWorld->SetVisAttributes(boxVisAtt);
    logicAirplaneTube->SetVisAttributes(airplaneTubeVisAtt);
    logicAirplaneSphere1->SetVisAttributes(airplaneSphere1VisAtt);
    logicAirplaneSphere2->SetVisAttributes(airplaneSphere2VisAtt);
    logicAirTube->SetVisAttributes(airTubeVisAtt);
    logicAirSphere1->SetVisAttributes(airSphere1VisAtt);
    logicAirSphere2->SetVisAttributes(airSphere2VisAtt);
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
  pSDman->AddNewDetector(mFDet1);                        // Register SD to SDManager.
  logicAirplaneTube->SetSensitiveDetector(mFDet1);       // Assign SD to the logical volume.
  logicAirTube->SetSensitiveDetector(mFDet1);            // Assign SD to the logical volume.
  logicAirplaneSphere1->SetSensitiveDetector(mFDet1);    // Assign SD to the logical volume.
  logicAirplaneSphere2->SetSensitiveDetector(mFDet1);    // Assign SD to the logical volume.
  logicAirSphere1->SetSensitiveDetector(mFDet1);         // Assign SD to the logical volume.
  logicAirSphere2->SetSensitiveDetector(mFDet1);         // Assign SD to the logical volume.

  
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
  //   Total, by protons, by neutrons.

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

