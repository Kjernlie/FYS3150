TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    jacobi.cpp \
    eigsym.cpp \
    schrodinger.cpp

HEADERS += \
    jacobi.h \
    eigsym.h \
    schrodinger.h

LIBS += -llapack -lblas -larmadillo
