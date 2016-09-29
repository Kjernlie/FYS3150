TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    jacobi.cpp \
    exb.cpp \
    eigsym.cpp

HEADERS += \
    jacobi.h \
    exb.h \
    eigsym.h

LIBS += -llapack -lblas -larmadillo
