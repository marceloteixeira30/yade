# File generated by kdevelop's qmake manager. 
# ------------------------------------------- 
# Subdir relative project main directory: ./plugins/CollisionFunctor/CollisionFunctor4ClosestFeatures/Sphere2Mesh2D4ClosestFeatures
# Target is a library:  

LIBS += -lMesh2D \
        -lSphere \
        -lSerialization \
        -lClosestFeatures \
        -lDistances \
        -lMath \
        -lGeometry \
        -lInteraction \
        -lMultiMethods \
        -rdynamic 
INCLUDEPATH = $(YADEINCLUDEPATH) 
MOC_DIR = $(YADECOMPILATIONPATH) 
UI_DIR = $(YADECOMPILATIONPATH) 
OBJECTS_DIR = $(YADECOMPILATIONPATH) 
QMAKE_LIBDIR = ../../../../plugins/GeometricalModel/Mesh2D/$(YADEDYNLIBPATH) \
               ../../../../plugins/GeometricalModel/Sphere/$(YADEDYNLIBPATH) \
               ../../../../toolboxes/Libraries/Serialization/$(YADEDYNLIBPATH) \
               ../../../../plugins/InteractionGeometry/ClosestFeatures/$(YADEDYNLIBPATH) \
               ../../../../toolboxes/ComputationalGeometry/Distances/$(YADEDYNLIBPATH) \
               ../../../../toolboxes/Libraries/Math/$(YADEDYNLIBPATH) \
               ../../../../yade/Geometry/$(YADEDYNLIBPATH) \
               ../../../../yade/Interaction/$(YADEDYNLIBPATH) \
               ../../../../toolboxes/Libraries/MultiMethods/$(YADEDYNLIBPATH) \
               $(YADEDYNLIBPATH) 
QMAKE_CXXFLAGS_RELEASE += -lpthread \
                          -pthread 
QMAKE_CXXFLAGS_DEBUG += -lpthread \
                        -pthread 
DESTDIR = $(YADEDYNLIBPATH) 
CONFIG += debug \
          warn_on \
          dll 
TEMPLATE = lib 
HEADERS += Sphere2Mesh2D4ClosestFeatures.hpp 
SOURCES += Sphere2Mesh2D4ClosestFeatures.cpp 
