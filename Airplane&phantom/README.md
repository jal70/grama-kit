# README

Cada carpeta tiene los archivos necesarios para la definición de un ejemplo. Se incluye un modelo cilíndrico y uno esférico.

En la carpeta src el archivo DetectorConstruction.cc es clave para parametrizar la simulación. Allí se definen el espacio, la aeronave, los materiales y el phantom.

* El espacio es un paralelepípedo de 30 m de lado

```c++
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
```

* La aeronave modelo en el ejemplo "avion_esferico_glare" es una esfera de 7 m de radio y está hecha a partir de esferas concéntricas de diferentes materiales. Por ejemplo, la esfera exterior es de aluminio 2024 y tiene un grosor de 0.4 mm.

```c++
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
```

## Uso
* El contenido de alguno de los ejemplos se copia en la carpeta Ejemplo/aeronaveG4/
* Las dimensiones del espacio y la geometría afectan los parámetros en las siguientes funciones:
    * Volumen define el espacio de referencia donde se generará la muestra (clase Volumen).
    ```python
    Volumen(Nx = 100, Ny = 100, Nz = 200, Lx = 5, Ly = 5, Lz = 5, Vx = 0, Vy = 0, Vz = 0,unidad='m',conversion2metros=1.0)
    ```
    * Rayos: prepara una muestra de direcciones de arribo de partículas de acuerdo a una distribución angular establecida (cos) (clase Rayos).
    ```python
    Rayos(Ro = np.array([0,0,0]), thetaMin = 0, thetaMax = math.pi/3, phiMin = 0, phiMax = 2*math.pi, M = 15, zFin = -1,verb=True, VolDetector=True)
    ```
    * ProcesosG4: genera el script para automatizar la simulación de GEANT4 (clase ProcesosG4).
    ```python
    ProcesosG4(archivo,r=...)
    ``` 