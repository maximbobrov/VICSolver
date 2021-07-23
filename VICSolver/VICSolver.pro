TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11
SOURCES += \
        globals.cpp \
        leastsquaressolver.cpp \
        main.cpp
#QMAKE_CXXFLAGS_RELEASE += -O3 -ffast-math  -msse -std=c++11
QMAKE_CXXFLAGS_RELEASE += -O3

#QMAKE_LFLAGS += -O3 -ffast-math  -msse -std=c++11
QMAKE_LFLAGS += -O3

unix{
LIBS+=  -lGL -lGLU -lglut -lm
}
win32{
LIBS += -lopenGL32 -lGLU32 -lm
LIBS += -L$$PWD/my_lib -lglut32
}

HEADERS += \
    globals.h \
    leastsquaressolver.h
