TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    system.cpp

HEADERS += \
    system.h

LIBS += -llapack -lblas -larmadillo
