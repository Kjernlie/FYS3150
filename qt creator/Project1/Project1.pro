TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    exb.cpp \
    exc.cpp \
    exd.cpp \
    exe.cpp

LIBS += -llapack -lblas -larmadillo
