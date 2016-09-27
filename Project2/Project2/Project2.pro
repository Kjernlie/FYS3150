TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    jacobi.cpp

HEADERS += \
    jacobi.h

LIBS += -llapack -lblas -larmadillo
