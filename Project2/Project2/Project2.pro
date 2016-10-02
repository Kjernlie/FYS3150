TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    jacobi.cpp \
    eigsym.cpp \
    two_electrons.cpp \
    schrodinger.cpp

HEADERS += \
    jacobi.h \
    eigsym.h \
    schrodinger.h \
    two_electrons.h

LIBS += -llapack -lblas -larmadillo
