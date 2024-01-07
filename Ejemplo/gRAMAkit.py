import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import random as rnd
from mpl_toolkits.mplot3d import Axes3D
import os
#import pickle
import json
import codecs
import datetime
import matplotlib.colors as colors
#import seaborn as sb
from sklearn.cluster import KMeans
from scipy.stats import cosine
from statistics import mean
import time

##########################################
############## CLASS Volumen
##########################################

class Volumen:
    'Volumen físico'
    # Esta clase define el volumen donde vive el objeto de estudio (universo)
    # El universo es un paralelepipedo voxelizado, paralelo a los planos cartesianos en el primer octante
    # Nx,Ny,Nz son el numero de divisiones en cada eje del parelepipedo
    # Lx,Ly,Lz son las longitudes de arista de cada voxel
    # Vx,Vy,Vz son las coordenadas del origen del volumen medidas en longitud absoluta
    # Estos datos N,L y V estan resumidos en 'volumenInfo'
    ###
    # Propiedades
    ###
    # density es la cantidad de materia media por unidad de volumen en cada voxel. Inicialmente se fija en cero
    # unidad es el sistema de unidades para medir longitudes. Por defecto es metros
    # conversion2metros es el factor de conversion necesario para operar con las unidades establecidas por el usuario. Si se decide no usar metros, es necesario redefinir la escala manualmente
    def __init__(self,Nx = 100, Ny = 100, Nz = 200, Lx = 5, Ly = 5, Lz = 5, Vx = 0, Vy = 0, Vz = 0,unidad='m',conversion2metros=1.0):
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.Lx = Lx*conversion2metros
        self.Ly = Ly*conversion2metros
        self.Lz = Lz*conversion2metros
        self.Vx = Vx
        self.Vy = Vy
        self.Vz = Vz
        self.unidad = unidad
        self.conversion2metros=conversion2metros
        self.N = np.array([self.Nx,self.Ny,self.Nz])
        self.L = np.array([self.Lx,self.Ly,self.Lz])
        self.V = np.array([self.Vx,self.Vy,self.Vz])
        self.density=np.zeros([self.Nx,self.Ny,self.Nz])
        self.volumenInfo=[]
        self.volumenInfo.append('Volumen físico datos (N,L,V):')
        self.volumenInfo.append(self.N.tolist())
        self.volumenInfo.append(self.L.tolist())
        self.volumenInfo.append(self.V.tolist())
        self.densityInfo=[]
        self.densityRo=[]
        self.densityRf=[]
        self.densityR=0.0
        print("Volumen creado con celdas (Nx,Ny,Nz)=(", self.Nx, self.Ny, self.Nz,"), dimension de voxel (Lx,Ly,Lz)=(", self.Lx/self.conversion2metros, 
              self.Ly/self.conversion2metros, self.Lz/self.conversion2metros,
              "). Unidad de medida=",self.unidad,". Factor de conversion a metros (m)=",self.conversion2metros,". Todas distancias y coordenadas se suponen y operan en metros (m)")
    def celda(self,x: np.array):
        ###
        # La funcion celda indica el voxel al que pertenece un punto
        # Los voxels se indican con coordenadas enteras positivas (0,0,0) es el origen. (Nx-1,Ny-1,Nz-1) es el punto opuesto sobre la diagonal principal.
        ###
        # Primero se traslada el universo al origen de coordenadas (x -> x-V)
        desplazamiento=x-self.V
        # Luego se mide la parte entera de x /(longitud del voxel)
        desplazamiento=desplazamiento/self.L
        return [math.floor(xi) for xi in desplazamiento]
    def intersectBoundary(self,origin,director):
        ####
        # Esta funcion podria estar fuera de uso
        ####
        #Planos x=0 y x=NxLx
        listOfPoints = []
        if director[0]!=0:
            l=(self.V[0]-origin[0])/director[0]
            p=origin+l*director
            if 0<(p[1]-self.V[1])/(self.Ny*self.Ly)<1 and 0<(p[2]-self.V[2])/(self.Nz*self.Lz)<1:
                listOfPoints.append(p)
            l=(self.V[0]+self.Nx*self.Lx-origin[0])/director[0]
            p=origin+l*director
            if 0<(p[1]-self.V[1])/(self.Ny*self.Ly)<1 and 0<(p[2]-self.V[2])/(self.Nz*self.Lz)<1:
                listOfPoints.append(p)
        #Planos y=0 y y=NyLy
        if director[1]!=0:
            l=(self.V[1]-origin[1])/director[1]
            p=origin+l*director
            if 0<(p[0]-self.V[0])/(self.Nx*self.Lx)<1 and 0<(p[2]-self.V[2])/(self.Nz*self.Lz)<1:
                listOfPoints.append(p)
            l=(self.V[1]+self.Ny*self.Ly-origin[1])/director[1]
            p=origin+l*director
            if 0<(p[0]-self.V[0])/(self.Nx*self.Lx)<1 and 0<(p[2]-self.V[2])/(self.Nz*self.Lz)<1:
                listOfPoints.append(p)
        #Planos z=0 y z=NzLz
        if director[2]!=0:
            l=(self.V[2]-origin[2])/director[2]
            p=origin+l*director
            if 0<(p[0]-self.V[0])/(self.Nx*self.Lx)<1 and 0<(p[1]-self.V[1])/(self.Ny*self.Ly)<1:
                listOfPoints.append(p)
            l=(self.V[2]+self.Nz*self.Lz-origin[2])/director[2]
            p=origin+l*director
            if 0<(p[0]-self.V[0])/(self.Nx*self.Lx)<1 and 0<(p[1]-self.V[1])/(self.Ny*self.Ly)<1:
                listOfPoints.append(p)
        return listOfPoints
    def intersectVol(self,inOut=[],samples=10):
        # Esta funcion es necesaria en la clase 'muestra'
        # Su utilidad es dar una medida de cuantos voxels de materia de interes atraviesa una trayectoria recta definida por dos puntos incial y final
        # La variable inOut es un arreglo de 2 vectores que representan los puntos inicial y final de la trayectoria recta.
        # La variable samples es el número de muestras promedio que se toman dentro de un voxel, si el rayo lo atraviesa de forma perpendicular a una cara
        # materiaIntersectada es una medida de la cantidad de materia que atraviesa la linea recta entre los puntos inicial y final propuestos
        materiaIntersectada=0
        # Verificamos que haya puntos inicial y final y no más que eso
        if len(inOut)==2:
            # La variable N determina cuantas muestras se deben tomar sobre la trayectoria, de forma que haya un número de muestras = samples por cada voxel
            N=int(samples*np.linalg.norm(inOut[1]-inOut[0])/min(self.L))
            # Nos movemos desde inOut[0] hasta inOut[1] paso a paso N pasos
            for i in range(N):
                # Desplazamiento en sistema propio en metros
                desplazamiento=((N-i)*inOut[0]+i*inOut[1])/N-self.V
                # Desplazamiento medido en voxels
                desplazamiento=desplazamiento/self.L
                # celda indica las coordenadas de la celda específica donde se encuentra el punto desplazado
                celda=[math.floor(xi) for xi in desplazamiento]
                # Esta condición verifica que el voxel (la celda) esté dentro del universo
                if 0<=celda[0]<self.Nx and 0<=celda[1]<self.Ny and 0<=celda[2]<self.Nz:
                    # Se suma la densidad de materia en el voxel. si hay varias intersecciones dentro del voxel, se sumará varias veces.
                    materiaIntersectada=materiaIntersectada+self.density[celda[0]][celda[1]][celda[2]]
        # se regresa la cantidad de materia intersectada normalizada con el numero de muestras.
        # Si se atraviesa solo un voxel cúbico en dirección perpendicular a alguna cara, el resultado debe ser 1
        # Desviaciones de este valor indican que la intersección atravesó más o menos materia.
        # El uso más importante de este resultado es determinar si se atravesó materia o no
        # El peso no nulo que toma el valor retornado puede usarse para dar más importancia a unas trayectorias.
        return materiaIntersectada/samples
    #def cmVol(self):
        ####
        # Esta funcion podria estar fuera de uso
        ####
    #    cm=np.zeros([3])
    #    for i in range(self.Nx):
    #        for j in range(self.Ny):
    #            for k in range(self.Nz):
    #                cm=cm+np.array([(i+0.5)*self.Lx,(j+0.5)*self.Ly,(k+0.5)*self.Lz])
    #    return self.V+cm/(self.Nx*self.Ny*self.Nz)
    def set_density(self,LX=8,LY=8,PX0=250,PY0=250,PH=800,cxo=250,cyo=250,czo=200,cr=100,file=False,reset=False):
        # No usar. Ver set_voxel_mater()
        if reset:
            self.density=np.zeros([self.Nx,self.Ny,self.Nz])
        if not file:
            print('yes')
            for i in range(self.Nx):
                for j in range(self.Ny):
                    for k in range(self.Nz):
                        self.density[i][j][k]=materia_Par_Es(LX,LY,PX0,PY0,PH,cxo,cyo,czo,cr,(i+0.5)*self.Lx,(j+0.5)*self.Ly,(k+0.5)*self.Lz)
        else:
            for i in range(self.Nx):
                for j in range(self.Ny):
                    for k in range(self.Nz):
                        self.density[i][j][k]=1.0
    def set_voxel_mater(self,x: np.array,state=1,reset=False):
        # Define el valor de densidad de materia en el voxel correspondiente al punto x
        # se usa internamente en la función cigar()
        # x es el punto en coordenadas con respecto al origen exterior al universo
        # state es el valor de densidad que se definirá en el voxel x. Por defecto es 1 (hay materia)
        # reset es un valor lógico que determina si se debe eliminar toda la materia del universo
        if reset:
            # Si reset es True, se borran los valores previos de densidad
            self.density=np.zeros([self.Nx,self.Ny,self.Nz])
        # r es la coordenada del voxel que ocupa x relativa al universo (sistema propio). Ver función celda()
        r=self.celda(x)
        # se asigna el nuevo valor de densidad de materia al voxel
        self.density[r[0]][r[1]][r[2]]=state
    def cigar(self,ro: np.array, rf: np.array,r):
        # cigar ocupa con materia el volumen interior de un cilindro con dos tapas semiesféricas
        # ro y rf son los centros de las semiesferas en los extremos con respecto al observador externo
        # r es el radio en metros de las semiesferas y por lo tanto es el radio transversal del cilindro
        # Se hace un reset de la información de materia que podría haber
        self.densityInfo=[]
        self.densityRo=ro
        self.densityRf=rf
        self.densityR=r
        # Elevamos el radio al cudadrado y lo guardamos en una variable
        r2=r**2
        # Se define el desplazamiento longitudinal del cilindro entre sus extremos
        tAxis=rf-ro
        # Se define una varible para almacenar la longitud del cilindro
        lAxis=np.linalg.norm(tAxis)
        # Recorremos todos los voxelx del universo
        # se fija el plano x = ctte
        for i in range(self.Nx):
            # Se imprime el progreso del plano x
            #print("x=",i)
            # se intersecta con el cada plano y = ctte
            for j in range(self.Ny):
                # se intersecta con cada plano z =  ctte
                for k in range(self.Nz):
                    # (i,j,k) representan un voxel en el sistema de referencia propio del universo
                    # cell center of mass cmC
                    # Se calcula la posición del centro de masa del voxel cmC
                    cmC=self.V+np.array([(i+.5)*self.Lx,(j+.5)*self.Ly,(k+.5)*self.Lz])
                    # Si la longitud del cilindro es cero, entonces el volumes una esfera
                    # Si no es cero 
                    if lAxis!=0:
                        # Se calcula la proyección del vector cmC-ro sobre el desplazamiento longitudinal del cilindro
                        d=np.dot((cmC-ro),tAxis)/lAxis
                        # a es True si d es positivo y menor que la longitud del cilindro y si la distancia al desplazamiento del cilindro es menor que el radio
                        # o si la magnitud cmC-ro es menor que el radio 
                        # o si la magnitud de cmC*rf es menor que el radio
                        a=(d>0 and d<lAxis) and ((np.linalg.norm(cmC-ro))**2 - d**2 < r2) or np.linalg.norm(cmC-ro) < r or np.linalg.norm(cmC-rf) < r
                    # Sí la longitud del cilindr es cero
                    else:
                        # a es True si cmC está dentro de la esfera
                        a=np.linalg.norm(cmC-ro) < r                             
                    #print(a)
                    # Usamos la función set_voxel_mater para cambiar el estado de materia del voxel (i,j,k)
                    self.set_voxel_mater(cmC,(a))
        self.densityInfo.append('Density from function cigar on ' + str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        self.densityInfo.append('ro = ' + str(ro))
        self.densityInfo.append('rf = ' + str(rf))
        self.densityInfo.append('r = ' + str(r))
        self.densityInfo.append('Se supone que todas las coordenadas y las distancias están en metros')
        self.densityInfo.append(self.volumenInfo)
    def plotDensity(self):
        # Grafica la densidad de materia
        xx=[]
        yy=[]
        zz=[]
        for i in range(self.Nx):
                for j in range(self.Ny):
                    for k in range(self.Nz):
                        # Se selecciona un voxel. Si la densidad es 1, se agregan las coordenadas del centro a los arreglos xx, yy y zz en forma correspondiente
                        if self.density[i][j][k]:
                            xx.append((i+.5)*self.Lx+self.Vx)
                            yy.append((j+.5)*self.Ly+self.Vy)
                            zz.append((k+.5)*self.Lz+self.Vz)
        # Isometrica
        ##
        #Set colours and render
        mcd=np.gcd.reduce(list(map(int,np.multiply(self.N,self.L)))) #https://stackoverflow.com/questions/29194588/python-gcd-for-list
        fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        ax = fig.add_subplot(projection='3d')
        ax.set_box_aspect(aspect = list(map(int,np.multiply(self.N,self.L)))/mcd)  #(2,2,1))
        #ax.plot_surface(x, y, z,  rstride=1, cstride=1, color='c', alpha=0.1, linewidth=0)
        ax.scatter(xx,yy,zz,color="b", alpha=0.1,s=20)
        ax.set_xlim([self.Vx,self.Vx+self.Nx*self.Lx])
        ax.set_ylim([self.Vy,self.Vy+self.Ny*self.Ly])
        ax.set_zlim([self.Vz,self.Vz+self.Nz*self.Lz])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        #ax.set_aspect("equal")
        plt.tight_layout()
        ax.set_title("Densidad de materia en proyección isométrica", va='bottom')
        plt.show()
    def saveDensity(self,name,force=False):
        if not(force) and os.path.isfile(name):
                print("Warning: El nombre de archivo elegido ya existe")
                print("Warning: Sobreescribir el archivo puede eliminar datos costosos")
                return print("Use la opción force=True si desea sobreescribirlo")
        b=self.density.tolist()
        densityAndMetadata=[]
        densityAndMetadata.append(self.densityInfo)
        densityAndMetadata.append(self.volumenInfo)
        densityAndMetadata.append(b)
        json.dump(densityAndMetadata, codecs.open(name, 'w', encoding='utf-8'), 
          separators=(',', ':'), 
          sort_keys=True, 
          indent=4)
    def loadDensity(self,name):
        if os.path.isfile(name):
            obj_text = codecs.open(name, 'r', encoding='utf-8').read()
            b_new = json.loads(obj_text)
            print(b_new[0])
            self.densityInfo=b_new[0]
            linea=np.array(b_new[1][1])
            print(linea)
            print(self.N)
            if linea[0]!=self.N[0] or linea[1]!=self.N[1] or linea[2]!=self.N[2]:
                return print('Error: Las dimensiones del arreglo a restaurar no coinciden con el volumen existente')
            self.density=np.array(b_new[2])

#############
#### Ejemplo volumen
#############
# Generar un volumen: 40m x 40m x 14m  
# El grid tiene las siguientes divisione en x,y,z: 64 x 64 x 32  
# Las longitudes de cada voxel en x,y,z son L: .625m x .625m x .4375m
#
# v=Volumen(64,64,32,625,625,437.5,0,0,-14,unidad='mm',conversion2metros=0.001)
            

##########################################
############## CLASS Rayos
##########################################
            

class Rayos:
    'Conjunto de rayos que convergen en un punto del espacio'
    def __init__(self,Ro = np.array([0,0,0]), thetaMin = 0, thetaMax = math.pi/3, phiMin = 0, 
                 phiMax = 2*math.pi, M = 15, zFin = -1,verb=True,VolDetector=True):
        ###################
        ### Intro
        ###################
        # Modelo de trayectorias de radiación en un punto del espacio
        # En cada punto del espacio hay un campo de radiación caracterizado por un espectro de tasa de flujo diferencial
        # La dirección vertical es z
        # Se usan coordenadas esféricas:
        # x=r cos phi sin theta
        # y=r sin phi sin theta
        # z=r cos theta
        # Ángulo polar minimo (debería ser 0) 
        self.thetaMin=thetaMin
        # Ángulo polar máximo (opción por defecto 60°)
        self.thetaMax=thetaMax
        # Rango de ángulo axial (0, 2pi)
        self.phiMin=phiMin
        self.phiMax=phiMax
        # M determina el número de puntos del grid en la superficie esférica.
        # Si thetaMax=pi/4 => el hemisferio de la esfera se corta en M paralelos con
        # separación de arco de angulo polar constante
        self.M=M
        # Ro es la posición del punto de referencia. En ese punto converge la radiación y alternativamente
        # podemos considerarlo una fuente caracterizada por su espectro angular
        # La radiación se propaga a partir de Ro y viaja hasta alcanzar una coordenada de referencia en el
        # eje z
        self.Ro = Ro
        # zFin es la coordenada z de referencia, que delimita el volumen de interés
        self.zFin=zFin
        #self.sinThetaMax=math.sin(thetaMax)
        # Delta diferencial de ángulo polar
        self.DeltaTheta=math.pi/(2*M)
        # Número de paralelos entre thetaMin y thetaMax separados en DeltaTheta
        self.m=int((self.thetaMax-self.thetaMin)/self.DeltaTheta)
        # Arreglo con delta diferencial de ángulo axial en cada paralelo (son m paralelos)
        # El valor por defecto es 2pi
        self.DeltaPhi=[2*math.pi]*self.m
        # Arreglo con el número de diferenciales de ángulo axial en cada paralelo.
        # en el polo norte sera n=1 por elección. Ese es el valor por defecto.
        self.nTheta=[1]*self.m
        # Arreglo de diferenciales de ángulo solido ($d\Omega=\sin{\theta} d\theta d\varphi$) en cada paralelo
        # el valor por defecto considera $d\varphi=2\pi$
        self.deltaSolidAngle=[self.DeltaTheta*2*math.pi*math.sin(self.DeltaTheta)]*self.m
        # Este ciclo calcula el número de divisiones del rango phi en cada paralelo.
        for i in range(1,self.m):
            # El número de divisiones de phi para cada theta se calcula de forma que compense el factor
            # sin(theta) en el ángulo sólido
            # La opción VolDetector=False agrega un factor cos(theta) a la medida que da 
            # cuenta de la proyección del flujo sobre una superficie paralela al plano xy
            if VolDetector:
                self.nTheta[i]=int(round(math.sin((i+1)*self.DeltaTheta)/math.sin(self.DeltaTheta),0))
                # El DeltaPhi característico en cada paralelo: se divide 2pi entre el número nTheta
                self.DeltaPhi[i]=self.DeltaPhi[i]/self.nTheta[i]
                # Diferencial de ángulo sólido en cada franja
                # Se puede usar para verificar la elección del grid o calcular integrales
                self.deltaSolidAngle[i]=math.sin((i+1)*self.DeltaTheta)*self.DeltaPhi[i]*self.DeltaTheta
            else:
                ############################################
                # Ojo, Ojo: esta opción no está actualizada
                # esta opción no se usa en este momento
                # debe diseñarse de nuevo
                ############################################
                self.nTheta[i]=int(round((math.sin((i+1)*self.DeltaTheta)*math.cos((i+1)*self.DeltaTheta))/(
                    math.sin(self.DeltaTheta)*math.cos(self.DeltaTheta)),0))
                # El DeltaPhi característico en cada paralelo: se divide 2pi entre el número nTheta
                self.DeltaPhi[i]=self.DeltaPhi[i]/self.nTheta[i]
                # Diferencial de ángulo sólido en cada franja
                # Se puede usar para verificar la elección del grid o calcular integrales
                self.deltaSolidAngle[i]=math.sin((i+1)*self.DeltaTheta)*self.DeltaPhi[i]*self.DeltaTheta
                #self.deltaSolidAngle[i]=math.sin((i+1)*self.DeltaTheta)*self.DeltaPhi[i]*(self.DeltaTheta*math.cos((i+1)*self.DeltaTheta))
        ########
        ## Para graficar distribucion de vectores directores de rayos
        ########
        # r es coordenada radial en plano xy
        self.r=[] 
        # phi es angulo polar en plano xy (es phi en coord esféricas)
        self.phi=[] 
        # xx es una lista de coordenadas x del vector director del rayo
        self.xx=[]
        # yy es una lista de coordenadas y del vector director del rayo
        self.yy=[]
        # zz es una lista de coordenadas z del vector director del rayo
        self.zz=[]
        ###############
        ## Para operar analíticamente
        ###############
        ## uRay es una lista de vectores directores. Tiene la misma información que xx, yy, zz,
        ## pero es una estructura más natural para el álgebra
        self.uRay=[]
        # uRayTheta es el theta que corresponde al rayo unitario, para calcular las normalizaciones
        self.uRayTheta=[]
        #########
        ## Genera distribución de rayos
        #########
        for i in range(len(self.nTheta)): # rango en theta
            ## Se selecciona el origen angular axial al azar para cada theta (desconectado)
            # phiInicial=rnd.random()*8*math.pi
            for j in range(self.nTheta[i]): # rango de angulos axiales
                # print(math.sin(i*self.DeltaTheta),j*self.DeltaPhi[i])
                # Distribución angular axial para cada coordenada polar
                # el origen i*10*math.pi/self.M se fija para mejorar la isotropia en phi
                self.phi.append(j*self.DeltaPhi[i]+i*10*math.pi/self.M) # 5 vueltas
                ## origen aleatorio (desconectado, ver phiInicial)
                #self.phi.append(j*self.DeltaPhi[i]+phiInicial)
                ## Dirección axial aleatoria (desconectado)
                #self.phi.append(j*self.DeltaPhi[i]+rnd.random()*2*math.pi)
                # Coordenada radial en plano xy
                self.r.append(math.sin(i*self.DeltaTheta))
                # Arreglos con componentes de vectores directores de cada rayo
                self.xx.append(math.cos(self.phi[-1])*math.sin(i*self.DeltaTheta))
                self.yy.append(math.sin(self.phi[-1])*math.sin(i*self.DeltaTheta))
                self.zz.append(-math.cos(i*self.DeltaTheta))
                ## El append del vector director en la lista no funciona
                self.uRay.append(np.array([self.xx[-1],self.yy[-1],self.zz[-1]]))
                self.uRayTheta.append(i)
        #########
        ## Guarda componenetes de la muestra de rayos en un dataframe
        #########
        self.rayosDf=pd.DataFrame({'ux': self.xx,'uy': self.yy,'uz': self.zz})
        ############
        ## Espectro diferencial fenomenológico (cos^2)
        ############
        # factorDensidad es un factor corrector de densidad de rayos para cada dirección cenital (theta)
        # Se fija un valor arbitrario 100 para los rayos en dirección thetaMax
        self.factorDensidad=[100]*self.m
        # Este ciclo calcula los valores de número de rayos en los theta<thetaMax
        # el angulo acimutal está integrado!!!! no hay discriminación de dirección angular
        ##########
        # Se toma en cuenta un factor de corrección cos(theta)^2 fenomenológico
        ##########
        # se toma en cuenta un factor de corrección por la variabilidad del ángulo sólido en la red usada
        ##########
        for i in range(self.m):
            self.factorDensidad[i]=int(self.factorDensidad[self.m-1]*
                                       (self.deltaSolidAngle[i]/self.deltaSolidAngle[-1])*
                                       ((math.cos(i*self.DeltaTheta)/(math.cos((self.m-1)*self.DeltaTheta)))**2))
        # normalizacionDensidadRayos es un factor de normalización para el espectro diferencial 
        # en una dirección angular
        self.normalizacionDensidadRayos=sum(np.multiply(self.factorDensidad,self.nTheta))
        ######
        ## Si hay verbosity: muestra datos iniciales y gráficas
        ######
        self.rayosInfo=[]
        self.rayosInfo.append('Haz cónico de rayos (M,thetaMin,thetaMax):')
        self.rayosInfo.append(self.M)
        self.rayosInfo.append(self.thetaMin)
        self.rayosInfo.append(self.thetaMax)
        if verb:
            print("Creado un haz cónico de rayos centrado en (Xo,Yo,Zo)=",self.Ro,self.m)
            print(self.nTheta,sum(self.nTheta))
            print(self.DeltaPhi)
            print(self.deltaSolidAngle)
            # Proyecccion plana
            # https://matplotlib.org/stable/gallery/pie_and_polar_charts/polar_demo.html
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
            ax.scatter(self.phi, self.r,s=1)
            ax.set_rmax(1)
            ax.set_rticks([0.25, 0.5, 0.75, 1])  # Less radial ticks
            ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
            ax.grid(True)
            ax.set_title("Vectores directores de rayos en plano xy", va='bottom')
            plt.show()
            ##
            # Isometrica
            # https://stackoverflow.com/questions/31768031/plotting-points-on-the-surface-of-a-sphere-in-pythons-matplotlib
            # Create a sphere
            r = 1
            pi = np.pi
            cos = np.cos
            sin = np.sin
            phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
            x = r*sin(phi)*cos(theta)
            y = r*sin(phi)*sin(theta)
            z = r*cos(phi)
            #Set colours and render
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(x, y, z,  rstride=1, cstride=1, color='c', alpha=0.1, linewidth=0)
            ax.scatter(self.xx,self.yy,self.zz,color="k",s=20)
            ax.set_xlim([-1,1])
            ax.set_ylim([-1,1])
            ax.set_zlim([-1,1])
            #ax.set_aspect("equal")
            plt.tight_layout()
            ax.set_title("Vectores directores de rayos en proyección isométrica", va='bottom')
            plt.show()
    def dispersionAnguloSolido(self):
        # Determina la dispersión relativa en los ángulos sólidos
        return np.std(self.deltaSolidAngle)/np.mean(self.deltaSolidAngle)
    def rayo(self,ds,rayLabel=int):
        #################
        ## Parametriza un rayo
        #################
        # Rayo generado por el vector director etiquetado por rayLabel
        # ds es el delta de longitud de arco en unidades de los ejes cartesianos
        # ray es una lista con las coordenadas del rayo, separadas en unidades de ds
        ray=[]
        # r es un punto sobre el rayo
        r=self.Ro
        # La condición del while garantiza que se generen puntos del rayo entre las coordenadas z del
        # volumen físico de interés (self.Ro[2] y zFin)
        while (r[2]==self.Ro[2]) or (r[2]<self.Ro[2])!=(r[2]<=self.zFin):
            # unimos el último punto generado a la lista
            ray.append(r)
            # genera un punto del rayo usando la ecuación paramétrica de la recta
            r=r+ds*self.uRay[rayLabel]
        return ray
    def centrarCono(self,r: np.array):
        ##############
        ## Cambia el origen del cono
        ##############
        # Esta función traslada el haz de partículas al punto r 
        self.Ro=r

#############
#### Ejemplo
#############
# Divide la esfera en 2M paralelos a separación angular constante  
# zFin es la profundidad máxima de propagación de los rayos. Está en metros
#
# r=Rayos(M=40,zFin=-14)
        
##########################################
############## CLASS Particles
##########################################
        
class particles:
    "Identificación de partículas"
    __pIDs = pd.read_table('./particleIDs.csv', delimiter=',',header=8)
    def __init__(self,A=0):
        self.A = A
    def summary(self):
        return self.__pIDs
    def isname(self,x):
        if sum(self.__pIDs['name']==x)!=0:
            return True
        else:
            return False
    def isG4GPS(self,x):
        if sum(self.__pIDs['G4 GPS']==x)!=0:
            return True
        else:
            return False
    def isCORSIKA(self,x):
        x=int(x)
        if sum(self.__pIDs['CORSIKA']==x)!=0:
            return True
        else:
            return False
    def isPDG(self,x):
        x=int(x)
        if sum(self.__pIDs['PYTHIA HERWIG (PDG)']==x)!=0:
            return True
        else:
            return False
    def CORSIKA2G4GPS(self,x):
        x=int(x)
        if sum(self.__pIDs['CORSIKA']==x)!=0:
            return self.__pIDs['G4 GPS'][self.__pIDs['CORSIKA']==x].to_frame().reset_index(drop=True)['G4 GPS'][0]
        else:
            return False
    def CORSIKA2NAME(self,x):
        x=int(x)
        if sum(self.__pIDs['CORSIKA']==x)!=0:
            return self.__pIDs['name'][self.__pIDs['CORSIKA']==x].to_frame().reset_index(drop=True)['name'][0]
        else:
            return False
particles=particles()

##########################################
############## CLASS histograma
##########################################

class histograma:
    "Histograma afluencia total de partículas"
    __PartID = particles
    __pin2names = {'p_in_bin(GeV)':'E_MeV','N_phot':'gamma','N_e+':'positron','N_e-':'electron','N_mu+':'muon_p',
           'N_mu-':'muon_m','N_pi0':'pion_0','N_pi+':'pion_p','N_pi-':'pion_m','N_n':'neutron',
           'N_p':'proton','N_pbar':'anti_proton'}
    def __init__(self):
        #Lista vacia para todas las partículas de la clase particles
        self.histogram=[ [None] for i in range(len(self.__PartID.summary()['PYTHIA HERWIG (PDG)'])) ]
        for i in range(len(self.histogram)):
            #Inicia lista vacía para la partícula i
            self.histogram[i]=[]
            #Codigo PDG para la partícula i
            self.histogram[i].append(self.__PartID.summary()['PYTHIA HERWIG (PDG)'][i])
            #Afluencia promedio para la partícula i (0.0 por defecto)
            self.histogram[i].append(0.0)
            #Inicia lista de bines de energía
            self.histogram[i].append([])
            #Inicia lista de afluencia por bin
            self.histogram[i].append([])
    def loadHST_Pinilla(self,name,_skiprows=3,_sep=' ',rebaja=0.01,nBin=20):
        print('Carga un histograma generado por Sergio Pinilla en su trabajo de grado')
        print('Opciones:')
        print('_skiprows=3: es el número de lineas comentadas que omite al principio del archivo')
        print('_sep=' ': es el separador de columnas')
        print('rebaja=0.01: es la fracción de partículas que va a tomar')
        print('Los histogramas de Pinilla son partículas por unidad de área en m^2 para todo el vuelo')
        print('Hay que tomar en cuenta la duración de cada vuelo para ajustar el parámetro "rebaja" de forma adecuada')
        self.histo=pd.read_csv(name,skiprows=_skiprows,sep=_sep)
        vec=list(self.histo.columns)
        del(vec[0])
        del(vec[0])
        del(self.histo[vec[-1]])
        del(self.histo[vec[-2]])
        self.histo.columns=vec
        self.histo.columns = [self.__pin2names.get(x, x) for x in self.histo.columns]
        self.histo['E_MeV']=self.histo['E_MeV']*1000
        self.histo.iloc[:,1:]=self.histo.iloc[:,1:]*rebaja
        self.nBin=nBin
    def saveHST(self,name='./histogramas/histograma_' + str(datetime.datetime.now().strftime('%Y%m%d_%H%M%S')) + '.hst'):
        self.histo.to_csv(name,sep=' ', index=False)
    def CORSIKALAGO(self,name,headerLines=8,deli=' ',nBin=20):
        muestra = pd.read_table(name, delimiter=deli,header=headerLines)
        del(muestra['prm_theta'])
        del(muestra['prm_phi'])
        muestra = muestra.rename(columns={'#':'CorsikaId', '#.1': 'px', 'CorsikaId':'py', 'px':'pz', 'py':'x', 'pz':'y', 'x':'z', 'y':'shower_id', 
                       'z': 'prm_id', 'shower_id':'prm_energy', 'prm_id':'prm_theta', 'prm_energy':'prm_phi'})
        Energias = muestra['px']**2+muestra['py']**2+muestra['pz']**2
        LEmin=math.floor(math.log10(min(Energias))/2)
        LEmax=math.floor(math.log10(max(Energias))/2)+1
        Energias=np.floor((np.log10(Energias)/2-LEmin)*nBin).astype(int)
        self.histo = pd.DataFrame(
            np.zeros([(LEmax-LEmin)*nBin,2+len(muestra['CorsikaId'].unique())]),columns=['E'] + 
            muestra['CorsikaId'].unique().tolist() + ['Total']).astype(int)
        for i in range((LEmax-LEmin)*nBin):
            self.histo['E'][i]=10**(LEmin+i/nBin+3)
        for i in range(len(Energias)):
            self.histo[muestra['CorsikaId'][i]][Energias[i]]+=1
        self.histo=self.histo.rename(lambda x: particles.CORSIKA2NAME(x) if (x!='E' and x!='Total') else x, axis=1)
        for i in range(len(self.histo)):
            self.histo.iloc[i,-1]=sum(self.histo.iloc[i,1:-1])
        self.specieWeight=[]
        for i in range(1,len(muestra['CorsikaId'].unique())+1):
            self.specieWeight.append(sum(self.histo.iloc[:,i]))
        self.nBin=nBin

##########################################
############## CLASS muestra
##########################################
        
class muestra:
    "Muestra de rayos que inciden sobre la frontera de un volumen"
    def __init__(self,Volumen,Rayos):
        self.Volumen=Volumen
        self.Rayos=Rayos
        self.maparayos=[]
        self.maparayosInfo=[]
        self.maparayosInfo.append(self.Volumen.densityInfo)
        self.maparayosInfo.append(self.Rayos.rayosInfo)
        self.planoRefCDF = np.zeros((self.Volumen.Nx,self.Volumen.Ny))
    def test(self):
        print(self.Volumen.Nx)
    def mapaRayos(self,res=2,rapido=True):
        # mapaRayos genera 
        # Se inicializa como una lista de listas bidimiensional vacia
        # Las dimensiones se arreglan con los pixels del plano de inyección
        mapaRayosb = [ [None]*self.Volumen.Nx for i in range(self.Volumen.Ny) ]
        # Recorremos el plano de inyección xy
        # x
        for i in range(0,self.Volumen.Nx,1):
            # Imprimimos la fila procesada para monitorear el avance
            print(i)
            # y
            for j in range(0,self.Volumen.Ny,1):
                #planoRef[i][j]=0
                # Iniciamos una lista en el pixel
                mapaRayosb[i][j]=[]
                # El primer elemento de la lista es la probabilidad condicional al pixel de rayos que intersectan materia
                # inicializamos en cero
                mapaRayosb[i][j].append(0.0)
                # Agregamos dos listas vacías al sitio i,j
                # la primera lista tiene las etiquetas de los rayos seleccionados
                mapaRayosb[i][j].append([])
                # la segunda lista tiene las probabilidades condicionales al pixel de los rayos seleccionados
                mapaRayosb[i][j].append([])
                # Rayos es un objeto de la clase rayos. Es un conjunto de rayos que pasan por el origen del
                # plano cartesiano
                # La operación centrarCono traslada la convergencia del conjunto de rayos al punto seleccionado
                # Trasladamos el cono al centro del pixel
                self.Rayos.centrarCono(np.array([(i+.5)*self.Volumen.Lx+self.Volumen.Vx,(j+.5)*self.Volumen.Ly+self.Volumen.Vy,0]))
                #print(i,j,r.Ro)
                # Recorremos todos los rayos del cono
                for k in range(len(self.Rayos.uRay)):
                    # rayo es una lista de puntos sobre el rayo
                    # esta parametrización atraviesa todo el volumen físico
                    # Los puntos del rayo están a una distancia = (alto del volumen)/res
                    # el valor por defecto de res es 2 para una estimación rápida
                    rayo=self.Rayos.rayo(self.Volumen.Nz*self.Volumen.Lz/res,k)
                    #Se verifica si el rayo intersecta materia antes de registrarlo
                    if self.Volumen.intersectVol([rayo[0],rayo[-1]],2)!=0.0:
                        # La primera componente de la lista en en pixel es la densidad de probabilidad de los rayos del pixel
                        mapaRayosb[i][j][0]+=self.Rayos.factorDensidad[self.Rayos.uRayTheta[k]]/self.Rayos.normalizacionDensidadRayos
                        # La segunda componente es la lista de etiquetas de los rayos. Se suma el rayo k
                        mapaRayosb[i][j][1].append(k)
                        # La tercera componente es la probabilidad del rayo 
                        mapaRayosb[i][j][2].append(self.Rayos.factorDensidad[self.Rayos.uRayTheta[k]]/self.Rayos.normalizacionDensidadRayos)
        # Renombramos la lista de listas
        self.maparayos=mapaRayosb
        # El arreglo planoRef tiene las probabilidades condicionales de los rayos efectivos en cada pixel 
        # del plano de inyección
        # Inicializamos en probabilidad cero
        self.planoRef = np.zeros((self.Volumen.Nx,self.Volumen.Ny))
        # Recorremos el plano de inyección
        # x
        for i in range(self.Volumen.Nx):
            # y
            for j in range(self.Volumen.Ny):
                # Se igualan los arreglos para cada pixel
                self.planoRef[i][j]=self.maparayos[i][j][0]
    def mapaRayos2(self,res=2,rapido=True):
        # mapaRayos2 produce un arreglo de vectores lógicos que indican cuando un rayo intersecta materia en el volumen
        # la estructura que obtendremos se llama planoCono
        # Las dimensiones se arreglan con los pixels del plano de inyección
        mapaRayosb = [ [None]*self.Volumen.Nx for i in range(self.Volumen.Ny) ]
        #mapaRayosLin = np.array([None]*self.Volumen.Nx*self.Volumen.Ny)
        # Recorremos el plano de inyección xy
        # x
        for i in range(0,self.Volumen.Nx,1):
            # Imprimimos la fila procesada para monitorear el avance
            print(i)
            # y
            for j in range(0,self.Volumen.Ny,1):
                # Iniciamos una lista en el pixel
                mapaRayosb[i][j]=[]
                # Rayos es un objeto de la clase rayos. Es un conjunto de rayos que pasan por el origen del
                # plano cartesiano
                # La operación centrarCono traslada la convergencia del conjunto de rayos al punto seleccionado
                # Trasladamos el cono al centro del pixel
                self.Rayos.centrarCono(np.array([(i+.5)*self.Volumen.Lx+self.Volumen.Vx,(j+.5)*self.Volumen.Ly+self.Volumen.Vy,0]))
                #print(i,j,r.Ro)
                # Recorremos todos los rayos del cono
                for k in range(len(self.Rayos.uRay)):
                    # rayo es una lista de puntos sobre el rayo
                    # esta parametrización atraviesa todo el volumen físico
                    # Los puntos del rayo están a una distancia = (alto del volumen)/res
                    # el valor por defecto de res es 2 para una estimación rápida
                    rayo=self.Rayos.rayo(self.Volumen.Nz*self.Volumen.Lz/res,k)
                    #Se verifica si el rayo intersecta materia antes de registrarlo
                    mapaRayosb[i][j].append(int(self.Volumen.intersectVol([rayo[0],rayo[-1]],2)!=0.0))
                #mapaRayosLin[i*self.Volumen.Ny+j] = np.array(mapaRayos[i][j])
        #self.lineaCono=np.stack(mapaRayosLin)
        self.planoCono=mapaRayosb
    def mapaRayosEnEsfera(self,error=0.0):
        if(not((self.Volumen.densityRo==self.Volumen.densityRf).all())):
            print("El volumen objetivo no es una esfera")
            return "err"
        radio=self.Volumen.densityR
        centro=self.Volumen.densityRo
        # mapaRayosEnEsfera produce un arreglo de vectores lógicos que indican cuando un rayo intersecta materia en un volumen esférico
        # la estructura que obtendremos se llama planoCono
        # Las dimensiones se arreglan con los pixels del plano de inyección
        mapaRayosb = [ [None]*self.Volumen.Nx for i in range(self.Volumen.Ny) ]
        #mapaRayosLin = np.array([None]*self.Volumen.Nx*self.Volumen.Ny)
        # Recorremos el plano de inyección xy
        # x
        for i in range(0,self.Volumen.Nx,1):
            # Imprimimos la fila procesada para monitorear el avance
            print(i)
            # y
            for j in range(0,self.Volumen.Ny,1):
                # Iniciamos una lista en el pixel
                mapaRayosb[i][j]=[]
                # Rayos es un objeto de la clase rayos. Es un conjunto de rayos que pasan por el origen del
                # plano cartesiano
                # La operación centrarCono traslada la convergencia del conjunto de rayos al punto seleccionado
                # Trasladamos el cono al centro del pixel
                self.Rayos.centrarCono(np.array([(i+.5)*self.Volumen.Lx+self.Volumen.Vx,(j+.5)*self.Volumen.Ly+self.Volumen.Vy,0]))
                #print(i,j,r.Ro)
                # Recorremos todos los rayos del cono
                for k in range(len(self.Rayos.uRay)):                  
                    dist=np.linalg.norm(np.cross(self.Rayos.uRay[k],(self.Rayos.Ro-centro)))
                    #print(dist)
                    mapaRayosb[i][j].append(int(dist<(radio+error)))
                    if(dist<(radio+error)):
                        print(k,dist)
                    # rayo es una lista de puntos sobre el rayo
                    # esta parametrización atraviesa todo el volumen físico
                    # Los puntos del rayo están a una distancia = (alto del volumen)/res
                    # el valor por defecto de res es 2 para una estimación rápida
                    ##rayo=self.Rayos.rayo(self.Volumen.Nz*self.Volumen.Lz/res,k)
                    #Se verifica si el rayo intersecta materia antes de registrarlo
                    ##mapaRayosb[i][j].append(int(self.Volumen.intersectVol([rayo[0],rayo[-1]],2)!=0.0))
                #mapaRayosLin[i*self.Volumen.Ny+j] = np.array(mapaRayos[i][j])
        #self.lineaCono=np.stack(mapaRayosLin)
        self.planoCono=mapaRayosb
    def plano2linea(self):
        mapaRayosLin = np.array([None]*self.Volumen.Nx*self.Volumen.Ny)
        for i in range(self.Volumen.Nx):
            for j in range(self.Volumen.Ny):
                mapaRayosLin[i*self.Volumen.Ny+j] = np.array(self.planoCono[i][j])
        self.lineaCono=np.stack(mapaRayosLin)
    def rotaP(self,p=[]):
        # Esta funcion perturba ligeramente la direccion de un vector en un eje aleatorio sobre el plano xy y genera una rotacion
        # grande en el eje z.
        thetaz, epsilonx = [(math.pi)*2*rnd.random()]*2, self.Rayos.DeltaTheta*(rnd.random()-.5)
        c, s = np.cos(thetaz), np.sin(thetaz)
        c2, s2 = np.cos(sum(thetaz)), np.sin(sum(thetaz))
        R = np.array(((c2,-s2,-epsilonx*s[1]),(s2,c2,epsilonx*c[1]),(-epsilonx*s[0],-epsilonx*c[0],1.0)))
        return np.dot(R,p)
    def rotaPinf(self,rayo,p=[]):
        # Esta funcion perturba ligeramente la direccion de un vector en un eje aleatorio sobre el plano xy y en el eje z.
        epsilonz, epsilonx = 2*math.pi/self.Rayos.nTheta[self.Rayos.uRayTheta[rayo]]*(rnd.random()-.5), self.Rayos.DeltaTheta*(rnd.random()-.5)
        R = np.array(((1.0,-epsilonz,0.0),(epsilonz,1.0,epsilonx),(0.0,-epsilonx,1.0)))
        return np.dot(R,p)
    def momentum2energy(self,p,m=0.0):
        return math.sqrt(m**2 + np.dot(p,p))
    def energy2absp(self,E,m=0):
        return math.sqrt(E**2-m**2)
    def normalizeP(self,p):
        norm = np.linalg.norm(p)
        if norm == 0: 
           return p
        return p / norm 
    def disparoMixP(self,particulas=[],isCORSIKA=True):
        numRayos=sum(self.Rayos.nTheta)
        listaSalida=[]
        CambioDeBase = np.array(((0.0,1.0,0.0),(1.0,0.0,0.0),(0.0,0.0,-1.0)))
        for particula in range(len(particulas)):
            #print(particula)
            #print(particulas[particula])
            #print(particulas[particula][4])
            # Extraemos el momentum. Lo llamaremos p
            p=np.array([particulas[particula][1],particulas[particula][2],particulas[particula][3]])
            p=np.dot(CambioDeBase,p)
            # Seleccionamos la direccion del momentum en nuestra plantilla de rayos
            # test es una lista que indica la clase de rayos a la que pertenece p
            # primero se inicializa la lista en ceros
            test=[0]*numRayos
            # luego se cambia a 1 el valor correspondiente al indice de la clase de rayos que maximiza p dot rayos
            claseRayoMax=self.Rayos.rayosDf.dot(p).idxmax()
            test[claseRayoMax]=1
            # identificamos los lugares donde la clase de rayos seleccionada instesecta el volumen de interes
            # para eso hacemos el producto escalar de test (que indica la clase de rayos) con lineaCono, que tiene listas de clases de rayos que intersectan el volumen de inter'es en cada sitio
            # lugares es un arreglo sobre los sitios. Si 1, hay un rayo que intersecta el volumen.
            lugares=np.dot(self.lineaCono,test)
            #print(len(lugares),set(lugares),sum(lugares))
            #print('lugares=',lugares)
            #print('lugares filtrado=',np.where(lugares!=0))
            # Generar las coordenadas de inyeccion
            # ver https://numpy.org/doc/stable/reference/arrays.nditer.html para entender el uso de no.nditer
            for lugar in np.nditer(np.where(lugares!=0)):
                newPart=[None]*8
                # Cambiamos el id de la particula
                if(isCORSIKA):
                    newPart[0]=particles.CORSIKA2G4GPS(particulas[particula][0])
                else:
                    newPart[0]=particulas[particula][0]
                # Perturbamos el momentum de la particula
                pRotada=self.rotaPinf(claseRayoMax,p)
                # Se calcula la energia
                newPart[1]=self.momentum2energy(pRotada)
                # se normaliza el momentum perturbado
                up=self.normalizeP(pRotada)
                # se guarda
                newPart[2]=up[0]
                newPart[3]=up[1]
                newPart[4]=up[2]
                # se sortean los lugares y se guardan
                newPart[5]=(lugar//self.Volumen.Ny+rnd.random())*self.Volumen.Lx
                newPart[6]=(lugar%self.Volumen.Ny+rnd.random())*self.Volumen.Ly
                newPart[7]=0
                listaSalida.append(newPart)
        listaSalida=pd.DataFrame(listaSalida)
        listaSalida.columns =['C_ID','E','px','py','pz','x','y','z']
        return listaSalida
    def disparoParalelo(self,particulas=[]):
        numRayos=sum(self.Rayos.nTheta)
        for particula in range(len(particulas)):
            #print(particula)
            #print(particulas[particula])
            #print(particulas[particula][4])
            # Extraemos el momentum. Lo llamaremos p
            p=np.array([particulas[particula][1],particulas[particula][2],particulas[particula][3]])
            # Seleccionamos la direccion del momentum en nuestra plantilla de rayos
            # test es una lista que indica la clase de rayos a la que pertenece p
            # primero se inicializa la lista en ceros
            test=[0]*numRayos
            # luego se cambia a 1 el valor correspondiente al indice de la clase de rayos que maximiza p dot rayos
            test[self.Rayos.rayosDf.dot(p).idxmax()]=1
            # identificamos los lugares donde la clase de rayos seleccionada instesecta el volumen de interes
            # para eso hacemos el producto escalar de test (que indica la clase de rayos) con lineaCono, que tiene listas de clases de rayos que intersectan el volumen de inter'es en cada sitio
            # lugares es un arreglo sobre los sitios. Si 1, hay un rayo que intersecta el volumen.
            lugares=np.dot(self.lineaCono,test)
            #print('lugares=',lugares)
            #print('lugares filtrado=',np.where(lugares!=0))
            # Generar las coordenadas de inyeccion
            # ver https://numpy.org/doc/stable/reference/arrays.nditer.html para entender el uso de no.nditer
            for lugar in np.nditer(np.where(lugares!=0)):
                #print(lugar)
                particulas[particula][4]=(lugares[lugar]//self.Volumen.Ny+rnd.random())*self.Volumen.Lx
                particulas[particula][5]=(particula%self.Volumen.Ny+rnd.random())*self.Volumen.Ly
                particulas[particula][6]=0
                print(particulas[particula])
    def saveRayos(self,name,force=False):
        if not(force) and os.path.isfile(name):
                print("Warning: El nombre de archivo elegido ya existe")
                print("Warning: Sobreescribir el archivo puede eliminar datos costosos")
                return print("Use la opción force=True si desea sobreescribirlo")
        maparayosAndMetadata=[]
        maparayosAndMetadata.append(self.maparayosInfo)
        maparayosAndMetadata.append(self.maparayos)
        json.dump(maparayosAndMetadata, codecs.open(name, 'w', encoding='utf-8'), 
          separators=(',', ':'), 
          sort_keys=True, 
          indent=4)
    def saveRayos2(self,name,force=False):
        if not(force) and os.path.isfile(name):
                print("Warning: El nombre de archivo elegido ya existe")
                print("Warning: Sobreescribir el archivo puede eliminar datos costosos")
                return print("Use la opción force=True si desea sobreescribirlo")
        maparayosAndMetadata=[]
        maparayosAndMetadata.append(self.maparayosInfo)
        maparayosAndMetadata.append(self.planoCono)
        json.dump(maparayosAndMetadata, codecs.open(name, 'w', encoding='utf-8'), 
          separators=(',', ':'), 
          sort_keys=True, 
          indent=4)
    def loadRayos(self,name):
        if os.path.isfile(name):
            obj_text = codecs.open(name, 'r', encoding='utf-8').read()
            b_new = json.loads(obj_text)
            print(b_new[0][0][5][1])
            self.maparayosInfo=b_new[0]
            linea=np.array(b_new[0][0][5][1])
            print(linea)
            print(self.Volumen.N)
            if linea[0]!=self.Volumen.N[0] or linea[1]!=self.Volumen.N[1] or linea[2]!=self.Volumen.N[2]:
                return print('Error: Las dimensiones del arreglo a restaurar no coinciden con el volumen existente')
            self.maparayos=b_new[1]
            self.planoRef = np.nan * np.empty((self.Volumen.Nx,self.Volumen.Ny))
            for i in range(self.Volumen.Nx):
                for j in range(self.Volumen.Ny):
                    self.planoRef[i][j]=self.maparayos[i][j][0]
    def loadRayos2(self,name):
        if os.path.isfile(name):
            obj_text = codecs.open(name, 'r', encoding='utf-8').read()
            b_new = json.loads(obj_text)
            print(b_new[0][0][5][1])
            self.maparayosInfo=b_new[0]
            linea=np.array(b_new[0][0][5][1])
            print(linea)
            print(self.Volumen.N)
            if linea[0]!=self.Volumen.N[0] or linea[1]!=self.Volumen.N[1] or linea[2]!=self.Volumen.N[2]:
                return print('Error: Las dimensiones del arreglo a restaurar no coinciden con el volumen existente')
            self.planoCono=b_new[1]
            #self.planoRef = np.nan * np.empty((self.Volumen.Nx,self.Volumen.Ny))
            #for i in range(self.Volumen.Nx):
            #    for j in range(self.Volumen.Ny):
            #        self.planoRef[i][j]=self.maparayos[i][j][0]
    def plot(self):
        plt.imshow(self.planoRef, cmap='viridis', interpolation='nearest')
        plt.colorbar()
        plt.show()
    def CDF(self,clusters=10):
        # Cumulative distribution function
        # Divide el plano de inyección en zonas de probabilidad parecida
        # Número de zonas = clusters
        # Integra la probabilidad sumando desde la zona máxima al minimo
        # Actualiza un arrglo self.planoRefCDF
        # cada elemento i,j del arreglo tiene la probabilidad acumulada
        # El valor por defecto de self.planoRefCDF es cero (ver __init__)
        ######
        # Se normaliza la probabilidad de cada pixel 
        suma=0.0
        # Recorremos el arreglo
        for i in range(self.Volumen.Nx):
            for j in range(self.Volumen.Ny):
                # Se suma el valor planoREF de cada elemento
                # planoRef es la probabilidad condicional al pixel de los rayos que intersecten materia 
                suma+=self.planoRef[i][j]
        #print(suma)
        # Se define self.planoRefPDF (Probability Density Function) para tener la probabilidad normalizada
        self.planoRefPDF=self.planoRef/suma
        # Definimos la lista planoRefPDF_dataFrame con los elementos de planoRefPDF 
        # en orden (0,0), (0,1) ... (1,0), (1,1), ...
        self.planoRefPDF_dataFrame=[]
        # Recorremos el plano de inyección xy
        # x
        for i in range(self.Volumen.Nx):
            # y
            for j in range(self.Volumen.Ny):
                # Se asignan los elementos de planoRefPDF en orden (0,0), (0,1) ... (1,0), (1,1), ...
                self.planoRefPDF_dataFrame.append(self.planoRefPDF[i][j])
        # Convertimos la lista de probabilidades en un dataframe
        self.planoRefPDF_dataFrame=pd.DataFrame(self.planoRefPDF_dataFrame)
        # Agrupamos las probabilidades usando kmeans con un número de grupos igual a "clusters"
        # el valor por defecto de clusters es 10
        kmeans = KMeans(n_clusters=clusters).fit(self.planoRefPDF_dataFrame)
        # Obtenemos los centroides de cada grupo
        centroids = kmeans.cluster_centers_
        # Preparamos arreglos para calcular CDF, ordenar, etc
        # intProbabilidadClusterK es un arreglo que contiene la probabilidad total de cada cluster de densidades
        intProbabilidadClusterK=[0.0]*len(centroids)
        # CumProbabilidadClusterK es la CDF hasta el i-ésimo cluster ordenados por probabilidad
        CumProbabilidadClusterK=[0.0]*len(centroids)
        # planoRefCluster es un arreglo que almacena el centroide correspondiente al cluster de probabilidad
        # del pixel
        self.planoRefCluster = np.nan * np.empty((self.Volumen.Nx,self.Volumen.Ny))
        # Generamos la lista de los indices del arreglo centroids ordenada desde la probabilidad mayor a la menor
        listaDescendente = sorted(range(len(centroids)), reverse=True, key=lambda k: centroids[k])
        
        # Recorremos el plano de inyección
        # x
        for i in range(self.Volumen.Nx):
            # y
            for j in range(self.Volumen.Ny):
                # se asigna a planoRefCluster el valor del centroide del conjunto de probabilidad 
                # correspondiente al pixel  
                self.planoRefCluster[i][j]=centroids[kmeans.predict(self.planoRefPDF_dataFrame)[i*self.Volumen.Nx 
                                                                                                + j]]
                # Vamos sumando la densidad de probabilidad en cada cluster
                intProbabilidadClusterK[kmeans.predict(self.planoRefPDF_dataFrame)[i*self.Volumen.Nx 
                                                                                   + j]]+=self.planoRefPDF[i][j]
        # La suma de las probabilidades de todos los clusters debe ser 1
        sum(intProbabilidadClusterK)
        # Hacemos la CDF sobre los clusters ordenados de mayor a menor probabilidad
        # cum es la variable que mantiene la CDF. Inicializamos en cero
        cum=0.0
        # Recorremos la lista de clusters de mayor a menor
        for i in range(len(listaDescendente)):
            # Imprime el valor del centroide iésimo
            print(centroids[listaDescendente[i]])
            # Suma la probabilidad de los pixels del cluster a la CDF
            cum+=intProbabilidadClusterK[listaDescendente[i]]
            # Registra la probabilidad acumulada hasta el i-ésimo cluster inclusive
            CumProbabilidadClusterK[listaDescendente[i]]=cum
            # Imprime la CDF hasta el i-ésimo cluster
            print(CumProbabilidadClusterK[listaDescendente[i]])
        # Actualizamos el arreglo planoRefCDF que registra la CDF correspondiente a cada pixel, de acuerdo 
        # al cluster.
        # Recorremos el plano de inyección
        # x
        for i in range(self.Volumen.Nx):
            # y
            for j in range(self.Volumen.Ny):
                # Se actualiza el valor de planoRefCDF
                self.planoRefCDF[i][j]=CumProbabilidadClusterK[kmeans.predict(
                    self.planoRefPDF_dataFrame)[i*self.Volumen.Nx + j]]
    def generaMuestra(self, n, tipo=1,corte=1,Energia='NoEnergia',Especie='NoEspecie',
                      name='./muestras/muestra_' + str(datetime.datetime.now().strftime('%Y%m%d_%H%M%S')),
                     flujoIntegral2vertical=3.0/(2*math.pi)):
        # n es el número promedio de partículas por pixel
        # corte es el porcentaje de partículas simuladas. Se eliminan las partículas de las zonas con menor probabilidad de impactar en el volumen
        # tipo es el tipo de partícula (ver corsika)
        # E es la energía cinética
        name=str(name) + '_' + str(Energia) + '_' + str(Especie) + '.csv'
        print(name)
        if n!=0:
            n=n*self.Volumen.Lx*self.Volumen.Ly*flujoIntegral2vertical
            ### print(n,self.Volumen.Lx,self.Volumen.Ly)
            with open(name, 'w') as f1:
                f1.write('#######################'+ os.linesep)
                f1.write('## Rayos'+ os.linesep)
                f1.write('##' + str(self.maparayosInfo)+ os.linesep)
                f1.write('## Coordenadas en metros, momentum espacial es unitario ' + os.linesep)
                f1.write('## Energía cinética en MeV ' + os.linesep)
                f1.write('## Columnas:'+ os.linesep)
                f1.write('## C_ID,E,px,py,pz,x,y,z ' + os.linesep)
                self.totalPorPixel=np.random.poisson(n*self.planoRef)
                for i in range(self.Volumen.Nx):
                    for j in range(self.Volumen.Ny):
                        #selecciona direcciones de momentum
                        if self.totalPorPixel[i][j]!=0 and (self.planoRefCDF<=corte)[i][j]:
                            direccionesP=rnd.choices(population=self.maparayos[i][j][1],weights=np.divide(self.maparayos[i][j][2],self.maparayos[i][j][0]),k=self.totalPorPixel[i][j])
                            for l in range(self.totalPorPixel[i][j]):
                                #selecciona punto de inyección
                                x=(i+np.random.random())*self.Volumen.Lx + self.Volumen.Vx
                                y=(j+np.random.random())*self.Volumen.Ly + self.Volumen.Vy
                                z=self.Volumen.Nz*self.Volumen.Lz + self.Volumen.Vz
                                #momentum
                                px=self.Rayos.uRay[direccionesP[l]][0]
                                py=self.Rayos.uRay[direccionesP[l]][1]
                                pz=self.Rayos.uRay[direccionesP[l]][2]
                                # print(Especie,Energia,x,y,z,px,py,pz)
                                #f1.write(f'{tipo:04}' + ',' + str(px) + ',' + str(py) + ',' + str(pz) + ',' + str(x) + ',' + str(y) + ',' + str(z) + os.linesep)
                                f1.write(str(Especie)  + ',' + str(Energia) + ',' + str(px) + 
                                         ',' + str(py) + ',' + str(pz) + ',' + str(x) + ',' + 
                                         str(y) + ',' + str(z) + os.linesep)
                f1.write(os.linesep)     
    def rebajaBordes(self,corte=.99):
        # Esta función muestra un mapa de calor con los puntos de inyección hasta
        # una CDF igual al parámetro "corte"
        # Se inicializa en cero el arreglo a graficar
        planoRestringido=np.zeros((self.Volumen.Nx,self.Volumen.Ny))
        # Recorremos el plano xy
        # x
        for i in range(self.Volumen.Nx):
            # y
            for j in range(self.Volumen.Ny):
                # Se verifica si el valor de CDF de cada punto es inferior al corte
                # El valor por defecto de self.planoRefCDF es cero (ver __init__) no descarta ningún punto
                if (self.planoRefCDF<=corte)[i][j]:
                    # Si el punto se preserva, se asigna el valor de planoRef a planoRestrigido
                    # planoRef es la probabilidad condicional al pixel de los rayos que intersecten materia 
                    planoRestringido[i][j]=self.planoRef[i][j]
        # Se grafica el heatmap
        plt.imshow(planoRestringido, cmap='viridis', interpolation='nearest')
        plt.colorbar()
        plt.show()
    def generaMuestraHST(self,h,corte=1,name='',flujoIntegral2vertical=3.0/(2*math.pi)):
        # corte es el porcentaje de partículas simuladas. Se eliminan las partículas de las zonas con menor probabilidad de impactar en el volumen
        # h es un objeto clase histograma
        # flujoIntegral2vertical es el factor de corrección que convierte el flujo integral en flujo vertical
        # el valor de flujoIntegral2vertical por defecto supone que la distribución cenital es cos2 y el 
        # ángulo está integrado en media esfera (0,\pi/2) para detector volumétrico (\Omega=sin\theta d\thetad\varphi)
        if name=='':
            name='./muestras/muestra_h_' + str(datetime.datetime.now().strftime('%Y%m%d_%H%M%S')) + '.csv'
        print(name)
        # Abrimos el archivo
        # Formato "muestra_h_YYYYMMDD_HHMMSS.csv"
        with open(name, 'w') as f1:
            f1.write('#######################'+ os.linesep)
            f1.write('## Rayos'+ os.linesep)
            f1.write('##' + str(self.maparayosInfo)+ os.linesep)
            f1.write('## Coordenadas en metros, momentum espacial es unitario ' + os.linesep)
            f1.write('## Energía cinética en MeV ' + os.linesep)
            f1.write('## Columnas:'+ os.linesep)
            f1.write('## C_ID,E,px,py,pz,x,y,z ' + os.linesep)
            # Recorremos las columnas del histograma (especies)
            for Especie in h.histo.columns[1:]:
                # Verificamos si el nombre de la columna está en formato G4 GPS y si no, saltamos
                if particles.isname(Especie):
                    #Especie=k
                    # Recorremos las filas
                    for m in range(len(h.histo)):
                        # Energía está en columna 0
                        Energia=h.histo.iloc[m,0]
                        # El número esperado de partículas está en la columna de la especie seleccionada
                        n=h.histo[Especie][m]
                        # Verificamos que haya partículas esperadas y si no, saltamos
                        if n!=0:
                            # Renormalizamos el promedio de partículas
                            # n viene en m^{-2}
                            # Se multiplica por el área del pixel donde se inyectará (Lx\timesLy)
                            # Se multiplica por el factor de corrección que convierte el flujo integral
                            # en flujo vertical
                            n=n*self.Volumen.Lx*self.Volumen.Ly*flujoIntegral2vertical
                            #print(n,self.Volumen.Lx,self.Volumen.Ly)
                            #print(Especie,Energia,n)
                            # Se define el número de partículas a inyectar en cada pixel
                            # Se usa una distribución de poisson con promedio igual a n
                            self.totalPorPixel=np.random.poisson(n*self.planoRef)
                            # Recorremos el plano xy en la superficie de inyección
                            # eje x
                            for i in range(self.Volumen.Nx):
                                # eje y
                                for j in range(self.Volumen.Ny):
                                    #selecciona direcciones de momentum
                                    # Condiciones para considerar el pixel
                                    # 1. Hay partículas inyectadas
                                    # 2. El pixel está en la zona incluída en planoRefCDF (Cumulative Density Function)
                                    if self.totalPorPixel[i][j]!=0 and (self.planoRefCDF<=corte)[i][j]:
                                        # Se seleccionan las direcciones
                                        # Hay una colección de rayos que intersectan el volumen de interés
                                        # Se escogen tantos rayos como 
                                        # el número de partículas inyectadas 
                                        # tomando pesos basados en la probabilidad condicional de
                                        # tener un rayo particular en el cuando se tiene alguno
                                        # self.maparayos[i][j][0] es la probabilidad de tener algún rayo del set
                                        # self.maparayos[i][j][2] es un arreglo con las probabilidades de los 
                                        # rayos posibles
                                        # self.maparayos[i][j][1] es un arreglo con las etiquetas de los rayos 
                                        direccionesP=rnd.choices(population=self.maparayos[i][j][1],
                                                                 weights=np.divide(self.maparayos[i][j][2],
                                                                                   self.maparayos[i][j][0]),
                                                                 k=self.totalPorPixel[i][j])
                                        # Recorremos los rayos generados
                                        for l in range(self.totalPorPixel[i][j]):
                                            #selecciona punto de inyección
                                            # en un lugar al azar dentro del pixel
                                            x=(i+np.random.random())*self.Volumen.Lx + self.Volumen.Vx
                                            y=(j+np.random.random())*self.Volumen.Ly + self.Volumen.Vy
                                            z=self.Volumen.Nz*self.Volumen.Lz + self.Volumen.Vz
                                            #momentum
                                            # Se anota la dirección del momentum usando la etiqueta de cada 
                                            # rayo escogido
                                            px=self.Rayos.uRay[direccionesP[l]][0]
                                            py=self.Rayos.uRay[direccionesP[l]][1]
                                            pz=self.Rayos.uRay[direccionesP[l]][2]
                                            #print(Especie,Energia,x,y,z,px,py,pz)
                                            # Se escribe al archivo
                                            # Agrgamos especie, energía, dirección de momentum y coordenadas
                                            # La energía se toma distribuída uniformemente e el logaritmo de la
                                            # energía del bin
                                            f1.write(str(Especie)  + ',' + str(Energia*10**(np.random.random()/h.nBin)) + ',' + str(px) + 
                                                     ',' + str(py) + ',' + str(pz) + ',' + str(x) + ',' + 
                                                     str(y) + ',' + str(z) + os.linesep)
    def sampleFromARTI(self,inFile,outFile,sampleSize=0.01,ifChunk=True):
        pd.DataFrame(columns=['C_ID','E','px','py','pz','x','y','z']).to_csv(outFile, header=True, index = False)
        submuestraTotal=sampleSize
        factorArea=self.Volumen.Lx*self.Volumen.Ly
        submuestreo=submuestraTotal*factorArea
        chunksize = int(10/(submuestreo))
        if ifChunk:
            for chunk in pd.read_csv(inFile, chunksize=chunksize, skiprows=8, sep=' '):
                sample=chunk.sample(frac=submuestreo)
                self.disparoMixP(sample.values.tolist()).to_csv(outFile, mode='a', header=False, index = False)
    def sampleFromHistogram(self,inFile,outFile='test.csv',sampleSize=0.01):
        #pd.DataFrame(columns=['C_ID','E','px','py','pz','x','y','z']).to_csv(outFile, header=True, index = False)
        submuestraTotal=sampleSize
        factorArea=self.Volumen.Lx*self.Volumen.Ly
        submuestreo=submuestraTotal*factorArea
        #chunksize = int(10/(submuestreo))
        histograma=pd.read_csv(inFile,sep=',',index_col=False)
        #print(submuestreo)
        #print(histograma)
        #print(histograma.columns[1:].to_list())
        #print(histograma['E_GeV'].to_list())
        name=outFile
        with open(name, 'w') as f1:
            f1.write("Centro Volumen Ro"+ os.linesep)
            f1.write(np.array2string(self.Volumen.densityRo, separator=', ')+ os.linesep)
            f1.write("Centro Volumen Rf"+ os.linesep)
            f1.write(np.array2string(self.Volumen.densityRf, separator=', ')+ os.linesep)
            f1.write("Radio Volumen r"+ os.linesep)
            f1.write(str(self.Volumen.densityR)+ os.linesep)
            f1.write("Submuestreo"+ os.linesep)
            f1.write(str(sampleSize)+ os.linesep)
            f1.write("Info"+ os.linesep)
            f1.write(str(self.maparayosInfo)+ os.linesep)
            f1.write('#'+ os.linesep)
            f1.write('C_ID,E,px,py,pz,x,y,z'+ os.linesep)
        for i in histograma.columns[1:].to_list():
            for j in range(len(histograma['E_GeV'].to_list())):
                n_particles=int(histograma[i][j]*submuestreo)
                #print(i,j,histograma[i][j],n_particles)
                if(n_particles!=0):
                    muestra=pd.DataFrame(columns=['theta','phi'])  #
                    muestra['theta']=abs(cosine.rvs(size=n_particles)/(2))
                    muestra['phi']=2*math.pi*(0.5-np.random.random(n_particles))
                    sample=pd.DataFrame(columns=['C_ID','px','py','pz','x','y','z'])
                    #print(muestra)
                    sample['C_ID']=[i]*n_particles
                    #sample['E']=histograma['E_GeV'][j]
                    p=self.energy2absp(histograma['E_GeV'][j])
                    sample['x']=[0.0]*n_particles
                    sample['y']=[0.0]*n_particles
                    sample['z']=[0.0]*n_particles
                    sample['px']=p*np.sin(muestra['theta'])*np.cos(muestra['phi'])
                    sample['py']=p*np.sin(muestra['theta'])*np.sin(muestra['phi'])
                    sample['pz']=p*np.cos(muestra['theta'])
                    #print(sample)
                    self.disparoMixP(sample.values.tolist(),isCORSIKA=False).to_csv(outFile, mode='a', header=False, index = False)
                
        
##########################################
############## CLASS ProcesosG4
##########################################
                
class ProcesosG4:
    'Archivos para facilitar corridas'
    def __init__(self,archivo,r=1.0,frac=1):
        ###################
        ### Intro
        ###################
        self.frac=frac
        with open(archivo) as input_file:
            self.head = [next(input_file) for _ in range(11)]
        self.df=pd.read_csv(archivo, sep=',', skiprows=11).sample(frac=frac)
        self.particulas=np.unique(self.df['C_ID']).tolist()
        print(self.particulas)
        self.d_salida=pd.DataFrame({'E':[],'Particula':[],'Angulo':[],'N_part':[],'x':[],'y':[],'z':[],'px':[],'py':[],'pz':[]})
        self.d_salida_vert=pd.DataFrame({'E':[],'Particula':[],'N_part':[],'x':[],'y':[],'z':[],'px':[],'py':[],'pz':[]})
        #Particulas = [], angulos = [], Energia = [], p = np.array(([0,0,0])), x = np.array(([0,0,0])) N_particulas=[]
        # Modelo de trayectorias de radiación en un punto del espacio
        # En cada punto del espacio hay un campo de radiación caracterizado por un espectro de tasa de flujo diferencial
        # La dirección vertical es z
        # Se usan coordenadas esféricas:
        # x=r cos phi sin theta
        # y=r sin phi sin theta
        # z=r cos theta
        # Ángulo polar minimo (debería ser 0) 
        #self.thetaMin=thetaMin
        # Ángulo polar máximo (opción por defecto 60°)
        #self.thetaMax=thetaMax
        # Rango de ángulo axial (0, 2pi)
        #self.phiMin=phiMin
        #self.phiMax=phiMax
        # M determina el número de puntos del grid en la superficie esférica.
        # Si thetaMax=pi/4 => el hemisferio de la esfera se corta en M paralelos con
        # separación de arco de angulo polar constante
        #self.M=M
        # Ro es la posición del punto de referencia. En ese punto converge la radiación y alternativamente
        # podemos considerarlo una fuente caracterizada por su espectro angular
        # La radiación se propaga a partir de Ro y viaja hasta alcanzar una coordenada de referencia en el
        # eje z
        #self.Ro = Ro
        # zFin es la coordenada z de referencia, que delimita el volumen de interés
        #self.zFin=zFin
        #self.sinThetaMax=math.sin(thetaMax)
        # Delta diferencial de ángulo polar
        #self.DeltaTheta=math.pi/(2*M)
        # Número de paralelos entre thetaMin y thetaMax separados en DeltaTheta
        #self.m=int((self.thetaMax-self.thetaMin)/self.DeltaTheta)
        # Arreglo con delta diferencial de ángulo axial en cada paralelo (son m paralelos)
        # El valor por defecto es 2pi
        #self.DeltaPhi=[2*math.pi]*self.m
        # Arreglo con el número de diferenciales de ángulo axial en cada paralelo.
        # en el polo norte sera n=1 por elección. Ese es el valor por defecto.
        #self.nTheta=[1]*self.m
        # Arreglo de diferenciales de ángulo solido ($d\Omega=\sin{\theta} d\theta d\varphi$) en cada paralelo
        # el valor por defecto considera $d\varphi=2\pi$
        #self.deltaSolidAngle=[self.DeltaTheta*2*math.pi*math.sin(self.DeltaTheta)]*self.m
        # Este ciclo calcula el número de divisiones del rango phi en cada paralelo.
        for particula in self.particulas:
            uVertical=np.array(([0,1,0]))
            print(particula)
            d_salida_parcial=pd.DataFrame({'E':[],'Particula':[],'Angulo':[],'N_part':[],'x':[],'y':[],'z':[],'px':[],'py':[],'pz':[]})
            d_salida_parcial_vert=pd.DataFrame({'E':[],'Particula':[],'N_part':[],'x':[],'y':[],'z':[],'px':[],'py':[],'pz':[]})
            df_part=self.df[self.df['C_ID']==particula].reset_index(drop=True)
            #print(df_part)
            # Versi'on rapida
            df_partRap=self.df[self.df['C_ID']==particula].reset_index(drop=True)
            countsEn_rapida, binsEn_rapida = np.histogram(np.log(df_part['E']*1000))
            labelsEn_rapida=binsEn_rapida[:-1] + np.diff(binsEn_rapida)/2
            d_salida_parcial_vert['E']=np.exp(labelsEn_rapida)
            d_salida_parcial_vert['Particula']=[particula]*len(countsEn_rapida)
            d_salida_parcial_vert['N_part']=countsEn_rapida
            d_salida_parcial_vert['x']=[r*uVertical[0]]*len(countsEn_rapida)
            d_salida_parcial_vert['y']=[r*uVertical[1]]*len(countsEn_rapida)
            d_salida_parcial_vert['z']=[r*uVertical[2]]*len(countsEn_rapida)
            d_salida_parcial_vert['px']=[float(-uVertical[0])]*len(countsEn_rapida)
            d_salida_parcial_vert['py']=[float(-uVertical[1])]*len(countsEn_rapida)
            d_salida_parcial_vert['pz']=[float(-uVertical[2])]*len(countsEn_rapida)
            self.d_salida_vert=pd.concat([self.d_salida_vert,d_salida_parcial_vert],axis=0, ignore_index = True)
            ##################
            countsCen, binsCen = np.histogram(np.arccos(-df_part['pz'])/math.pi)
            labelsCen=binsCen[:-1] + np.diff(binsCen)/2
            #print(len(labelsCen))
            df_part['Angulo'] = pd.cut(np.arccos(-df_part['pz'])/math.pi, bins=binsCen, labels=labelsCen)
            df_part.dropna(inplace=True)
            self.angulos=np.unique(df_part['Angulo']).tolist()
            for angulo in self.angulos:
                uVector=np.array(([0,math.cos(angulo*math.pi),math.sin(angulo*math.pi)]))
                df_ang=df_part[df_part['Angulo']==angulo].reset_index(drop=True)
                countsEn, binsEn = np.histogram(np.log(df_ang['E']*1000))
                labelsEn=binsEn[:-1] + np.diff(binsEn)/2
                #print(len(labelsEn))
                df_ang['E_binned'] = pd.cut(np.log(df_ang['E']*1000), bins=binsEn, labels=labelsEn)
                df_ang.dropna(inplace=True)
                d_salida_parcial['E']=np.exp(labelsEn)
                d_salida_parcial['Particula']=[particula]*len(countsEn)
                d_salida_parcial['Angulo']=[angulo]*len(countsEn)
                d_salida_parcial['N_part']=countsEn
                d_salida_parcial['x']=[r*uVector[0]]*len(countsEn)
                d_salida_parcial['y']=[r*uVector[1]]*len(countsEn)
                d_salida_parcial['z']=[r*uVector[2]]*len(countsEn)
                d_salida_parcial['px']=[float(-uVector[0])]*len(countsEn)
                d_salida_parcial['py']=[float(-uVector[1])]*len(countsEn)
                d_salida_parcial['pz']=[float(-uVector[2])]*len(countsEn)
                #print((d_salida_parcial))
                self.d_salida=pd.concat([self.d_salida,d_salida_parcial],axis=0, ignore_index = True)
        #self.d_salida.drop(self.d_salida[self.d_salida.N_part == 0].index, inplace=True)
    def macGun(self,file):
        file1 = open(file, "w")  # write mode
        file1.write("Class,Data \n")
        file1.close()
        file1 = open(file, "a")  # append mode
        for i in range(len(self.d_salida)):
            if self.d_salida['N_part'][i]!=0:
                linea='Class ' + str(i)+ ',' + '/gun/particle ' + str(self.d_salida['Particula'][i])
                file1.write(linea + "\n")
                linea='Class ' + str(i)+ ',' + '/gun/energy ' + str(self.d_salida['E'][i]) + ' MeV'
                file1.write(linea + "\n")
                linea='Class ' + str(i)+ ',' + '/gun/position ' + str(self.d_salida['x'][i]) + ' ' + str(self.d_salida['y'][i]) + ' '  + str(self.d_salida['z'][i]) +' m'
                file1.write(linea + "\n")
                linea='Class ' + str(i)+ ',' + '/gun/direction ' + str(self.d_salida['px'][i]) + ' ' + str(self.d_salida['py'][i]) + ' '  + str(self.d_salida['pz'][i])
                file1.write(linea + "\n")
                linea='Class ' + str(i)+ ',' + '/run/beamOn ' + str(int(self.d_salida['N_part'][i]))
                file1.write(linea + "\n")
        file1.close()
    def macGun_vert(self,file):
        file1 = open(file, "w")  # write mode
        for i in range(11):
            file1.write(str(self.head[i].strip())+ os.linesep)
        file1.write('Segundo submuestreo'+ os.linesep)
        file1.write(str(self.frac)+ os.linesep)
        file1.write('#'+ os.linesep)
        file1.write("Class,Data \n")
        file1.close()
        file1 = open(file, "a")  # append mode
        for i in range(len(self.d_salida_vert)):
            if self.d_salida_vert['N_part'][i]!=0:
                linea='Class ' + str(i)+ ',' + '/gun/particle ' + str(self.d_salida_vert['Particula'][i])
                file1.write(linea + "\n")
                linea='Class ' + str(i)+ ',' + '/gun/energy ' + str(self.d_salida_vert['E'][i]) + ' MeV'
                file1.write(linea + "\n")
                linea='Class ' + str(i)+ ',' + '/gun/position ' + str(self.d_salida_vert['x'][i]) + ' ' + str(self.d_salida_vert['y'][i]) + ' '  + str(self.d_salida_vert['z'][i]) +' m'
                file1.write(linea + "\n")
                linea='Class ' + str(i)+ ',' + '/gun/direction ' + str(self.d_salida_vert['px'][i]) + ' ' + str(self.d_salida_vert['py'][i]) + ' '  + str(self.d_salida_vert['pz'][i])
                file1.write(linea + "\n")
                linea='Class ' + str(i)+ ',' + '/run/beamOn ' + str(int(self.d_salida_vert['N_part'][i]))
                file1.write(linea + "\n")
        file1.close()

##########################################
############## CLASS visuales
##########################################
        
class visuales:
    'Visualizacion'
    def __init__(self,inputFile):
        self.Lx=float('nan')
        self.Ly=float('nan')
        self.Lz=float('nan')
        self.Y = pd.read_csv(inputFile,sep =' ', header=None,skiprows=15)
        self.Y['Etotal']=self.Y.iloc[:, 3:].sum(axis=1)
        self.Nx=int(max(self.Y[0])+1)
        self.Ny=int(max(self.Y[1])+1)
        self.Nz=int(max(self.Y[2])+1)
        ## en Z
        self.capasP=[]
        for i in range(0,self.Nx):
            capa=[]
            for j in range(self.Ny):
                capa.append(self.Y.loc[self.Y[0] == i].loc[self.Y[1] == j]['Etotal'].to_list())
            capa=pd.DataFrame(capa)
            self.capasP.append(capa)
        ## en Y
        self.capasL=[]
        for i in range(0,self.Ny):
            capa=[]
            for j in range(self.Nz):
                capa.append(self.Y.loc[self.Y[1] == i].loc[self.Y[2] == j]['Etotal'].to_list())
            capa=pd.DataFrame(capa)
            self.capasL.append(capa)
    def factors(self,n):
        _factors = []
        _factors.append(1)
        _factors.append(n)
        med = int(np.sqrt(n)) + 1
        for i in range(2,med):
            if n%i == 0:
                _factors.append(i)
                _factors.append(int(n/i))
        return(np.unique(_factors))
    def CortesEnZ(self):
        factores=self.factors(self.Nx)
        if len(factores[factores<=5])>1:
            col=factores[factores<=5][-1]
        else:
            col=factores[1]
        row=int(self.Nx/col)
        if self.Nx>30:
            col=5
            row=4
            muestreo=int(self.Nx/20)
        else:
            muestreo=1
        ## Mantiene la misma escala lineal de color para todos los cortes (excluye el cero)
        #norm = plt.Normalize(np.nanmin(Y['Etotal']), max(Y['Etotal']))
        ## Mantiene la misma escala log de color para todos los cortes
        norm=colors.LogNorm(np.nanmin(self.Y[self.Y['Etotal']>0]['Etotal']), max(self.Y['Etotal']))
        fig, axs = plt.subplots(row, col,figsize=(20,20))
        for i in range(row):
            for j in range(col):
                im = axs[i, j].contourf(self.capasP[(i*col+j)*muestreo],norm=norm)
                #im = axs[i, j].imshow(self.capasP[(i*col+j)*muestreo],norm=norm)
                axs[i, j].tick_params(
                    axis='both',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom=False,      # ticks along the bottom edge are off
                    top=False,  # ticks along the top edge are off
                    right=False,
                    left=False,
                    labelbottom=True,
                    labelleft=True)
                axs[i, j].set_xlabel('x')
                axs[i, j].set_ylabel('y')
                axs[i, j].title.set_text('z = ' + str((i*col+j)*muestreo))
                plt.colorbar(im, ax=axs[i, j],orientation='horizontal')
    def CortesEnY(self):
        factores=self.factors(self.Ny)
        if len(factores[factores<=5])>1:
            col=factores[factores<=5][-1]
        else:
            col=factores[1]
        row=int(self.Ny/col)
        if self.Ny>30:
            col=5
            row=4
            muestreo=int(self.Ny/20)
        else:
            muestreo=1
        ## Mantiene la misma escala lineal de color para todos los cortes (excluye el cero)
        #norm = plt.Normalize(np.nanmin(Y['Etotal']), max(Y['Etotal']))
        ## Mantiene la misma escala log de color para todos los cortes
        norm=colors.LogNorm(np.nanmin(self.Y[self.Y['Etotal']>0]['Etotal']), max(self.Y['Etotal']))
        fig, axs = plt.subplots(row, col,figsize=(20,20))
        for i in range(row):
            for j in range(col):
                im = axs[i, j].contourf(self.capasL[(i*col+j)*muestreo],norm=norm)
                #im = axs[i, j].imshow(self.capasL[(i*col+j)*muestreo],norm=norm)
                axs[i, j].tick_params(
                    axis='both',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom=False,      # ticks along the bottom edge are off
                    top=False,  # ticks along the top edge are off
                    right=False,
                    left=False,
                    labelbottom=True,
                    labelleft=True)
                axs[i, j].set_xlabel('x')
                axs[i, j].set_ylabel('z')
                axs[i, j].title.set_text('y = ' + str((i*col+j)*muestreo))
                plt.colorbar(im, ax=axs[i, j],orientation='horizontal')
        

def wrneutron(E):
    if(E<1):
        return 2.5+18.2*math.exp((-(math.log(float(E)))**2)/6)
    elif E<=50:
        return 5.0+17.0*math.exp((-(math.log(2.0*float(E)))**2)/6)
    else:
        return 2.5+3.25*math.exp((-(math.log(.04*float(E)))**2)/6)

def dosisMediaEqEsfera(inputfile,radio,sample=1.0,MeV=True,Neutron=False):
    # Dosis equivalente en mu Sievert
    filePD=pd.read_csv(inputfile,sep =' ',skiprows=14)
    delta_equ_proton=filePD.filter(regex=("proton.*")).sum().sum()*(1.0)
    #print(delta_equ_proton)
    energy_n=np.array(filePD.filter(regex=("neutron.*")).sum().to_list())
    equivalente=[]
    for i in range(len(filePD.filter(regex=("neutron.*")).columns)):
        Deltaeq=wrneutron(float(filePD.filter(regex=("neutron.*")).columns[i].split("_")[1]))-1
        #print(filePD.filter(regex=("neutron.*")).columns[i].split("_")[1],Deltaeq)
        equivalente.append(Deltaeq)
    equivalente=np.array(equivalente)
    delta_equ_neutron=np.dot(equivalente,energy_n)
    #print(delta_equ_neutron)
    Etot=filePD.filter(regex=("._.")).sum().sum()
    #print(Etot)
    H=(Etot+delta_equ_neutron+delta_equ_proton)*1.602e-13*1e6/((4*math.pi*radio**3/3)*sample)
    if(Neutron):
        H=delta_equ_neutron=np.dot((equivalente+1),energy_n)*1.602e-13*1e6/((4*math.pi*radio**3/3)*sample)
    #H=(Etot+delta_equ_neutron+delta_equ_proton)/(q.Nx*q.Nz*q.Ny*.1*.1*.1)*1.602e-13*1e5*1e6/1.0
    return H

def dosisMediaEqEsferaErr(inputfile,radio,sample=1.0,MeV=True,Neutron=False):
    # Desviación estandar para dosis equivalente en mu Sievert
    filePDAbs=pd.read_csv(inputfile,sep =' ',skiprows=14)
    filePD=filePDAbs[['X','Y','Z']]
    for i in filePDAbs.filter(regex=("._.")).columns:
        filePD = pd.concat([filePD, filePDAbs[i]/math.sqrt(float(i.split("_")[2]))], axis=1)
    delta_equ_proton=filePD.filter(regex=("proton.*")).sum().sum()*(1.0)
    #print(delta_equ_proton)
    energy_n=np.array(filePD.filter(regex=("neutron.*")).sum().to_list())
    equivalente=[]
    for i in range(len(filePD.filter(regex=("neutron.*")).columns)):
        Deltaeq=wrneutron(float(filePD.filter(regex=("neutron.*")).columns[i].split("_")[1]))-1
        #print(filePD.filter(regex=("neutron.*")).columns[i].split("_")[1],Deltaeq)
        equivalente.append(Deltaeq)
    equivalente=np.array(equivalente)
    delta_equ_neutron=np.dot(equivalente,energy_n)
    #print(delta_equ_neutron)
    Etot=filePD.filter(regex=("._.")).sum().sum()
    #print(Etot)
    H=(Etot+delta_equ_neutron+delta_equ_proton)*1.602e-13*1e6/((4*math.pi*radio**3/3)*sample)
    if(Neutron):
        H=delta_equ_neutron=np.dot((equivalente+1),energy_n)*1.602e-13*1e6/((4*math.pi*radio**3/3)*sample)
    #H=(Etot+delta_equ_neutron+delta_equ_proton)/(q.Nx*q.Nz*q.Ny*.1*.1*.1)*1.602e-13*1e5*1e6/1.0
    return H

def lee(file,line_numbers = []):
    with open(file, 'r') as fp:
        # lines to read
        # line_numbers = [7]
        # To store lines
        lines = []
        for i, line in enumerate(fp):
            # read line 4 and 7
            if i in line_numbers:
                lines.append(line.strip())
            elif i > max(line_numbers):
                # don't read after max to save time
                break
    return lines
