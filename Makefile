#enable all warnings
FLAGS = -Wall
#use optimization level 3
FLAGS += -O3
#enable openmp for paralell processing
FLAGS += -fopenmp
#enable fast math
#FLAGS += -ffast-math

#to use for valgrind --leak-check=yes ./ComBinE
#or ddd ./ComBinE
FLAGS += -g

#ComBinE_SN: ComBinE.cpp objectfiles/tables.o objectfiles/structures.o objectfiles/setup.o objectfiles/evolution.o objectfiles/radii.o objectfiles/tracks.o objectfiles/gravitationalwave.o objectfiles/statistics.o objectfiles/galactic.o objectfiles/user.o
#	c++ $(FLAGS) ComBinE.cpp objectfiles/*.o -o SNShell/ComBinE_SN -std=gnu++0x

ComBinE: ComBinE.cpp objectfiles/tables.o objectfiles/structures.o objectfiles/setup.o objectfiles/evolution.o objectfiles/radii.o objectfiles/tracks.o objectfiles/gravitationalwave.o objectfiles/statistics.o objectfiles/galactic.o objectfiles/user.o
	c++ $(FLAGS) ComBinE.cpp objectfiles/*.o -o ComBinE -std=gnu++0x

objectfiles/tables.o: ComBinElib/ComBinElib.h ComBinElib/tables.cpp
	c++ $(FLAGS) -c ComBinElib/tables.cpp -o objectfiles/tables.o

objectfiles/structures.o: ComBinElib/ComBinElib.h ComBinElib/structures.cpp
	c++ $(FLAGS) -c ComBinElib/structures.cpp -o objectfiles/structures.o

objectfiles/setup.o: ComBinElib/ComBinElib.h ComBinElib/setup.cpp
	c++ $(FLAGS) -c ComBinElib/setup.cpp -o objectfiles/setup.o -std=gnu++0x

objectfiles/evolution.o: ComBinElib/ComBinElib.h ComBinElib/evolution.cpp
	c++ $(FLAGS) -c ComBinElib/evolution.cpp -o objectfiles/evolution.o

objectfiles/radii.o: ComBinElib/ComBinElib.h ComBinElib/radii.cpp
	c++ $(FLAGS) -c ComBinElib/radii.cpp -o objectfiles/radii.o

objectfiles/tracks.o: ComBinElib/ComBinElib.h ComBinElib/tracks.cpp
	c++ $(FLAGS) -c ComBinElib/tracks.cpp -o objectfiles/tracks.o

objectfiles/gravitationalwave.o: ComBinElib/ComBinElib.h ComBinElib/gravitationalwave.cpp
	c++ $(FLAGS) -c ComBinElib/gravitationalwave.cpp -o objectfiles/gravitationalwave.o

objectfiles/statistics.o: ComBinElib/ComBinElib.h ComBinElib/statistics.cpp
	c++ $(FLAGS) -c ComBinElib/statistics.cpp -o objectfiles/statistics.o

objectfiles/galactic.o: ComBinElib/ComBinElib.h ComBinElib/galactic.cpp
	c++ $(FLAGS) -c ComBinElib/galactic.cpp -o objectfiles/galactic.o

objectfiles/user.o: ComBinElib/ComBinElib.h ComBinElib/user.cpp
	c++ $(FLAGS) -c ComBinElib/user.cpp -o objectfiles/user.o

plot: plot.cpp
	c++ $(FLAGS) plot.cpp -o plot -lcpgplot -lpgplot -lX11 -lm

profile:
	rm -f gmon.out
	c++ $(FLAGS) -pg -c ComBinElib/tables.cpp -o objectfiles/tables.o
	c++ $(FLAGS) -pg -c ComBinElib/structures.cpp -o objectfiles/structures.o
	c++ $(FLAGS) -pg -c ComBinElib/setup.cpp -o objectfiles/setup.o -std=gnu++0x
	c++ $(FLAGS) -pg -c ComBinElib/evolution.cpp -o objectfiles/evolution.o
	c++ $(FLAGS) -pg -c ComBinElib/radii.cpp -o objectfiles/radii.o
	c++ $(FLAGS) -pg -c ComBinElib/tracks.cpp -o objectfiles/tracks.o
	c++ $(FLAGS) -pg -c ComBinElib/gravitationalwave.cpp -o objectfiles/gravitationalwave.o
	c++ $(FLAGS) -pg -c ComBinElib/statistics.cpp -o objectfiles/statistics.o
	c++ $(FLAGS) -pg -c ComBinElib/galactic.cpp -o objectfiles/galactic.o
	c++ $(FLAGS) -pg -c ComBinElib/user.cpp -o objectfiles/user.o
	c++ $(FLAGS) -pg ComBinE.cpp objectfiles/*.o -o ComBinE -std=gnu++0x
	ComBinE 12345678901 > output_profile.txt
	gprof ComBinE gmon.out > analysis.txt

clean:
	rm -f objectfiles/*.o
