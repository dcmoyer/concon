
CPP_COMPILER=/ifshome/dmoyer/gcc-5.2.0/bin/g++
FLAGS=-fPIC -std=c++11

default: all

dirs:
	mkdir -p bin

spherical_harmonics: dirs sh/spherical_harmonics.cc
	$(CPP_COMPILER) $(FLAGS) -c sh/spherical_harmonics.cc -I./sh -o bin/spherical_harmonics.o
	$(CPP_COMPILER) $(FLAGS) -shared sh/spherical_harmonics.cc -I./sh -o bin/spherical_harmonics.so
	ar rcs bin/spherical_harmonics.a bin/spherical_harmonics.o

all: spherical_harmonics




