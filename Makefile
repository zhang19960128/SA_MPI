CXX=mpiicc
CXXFLAGS=
DEPS= interface.h image.h atom.h readion.h readpara.h sa.h penalty.h simann.h
LIBPATH =
CXXFLAGS +=
CXXFLAGS +=-I./include
vpath %.cpp src
vpath %.h include
vpath %.o obj
sa.x: atom.o ewald.o lj12.o bv.o bvv.o image.o stress.o readion.o readpara.o main.o sa.o simann.o penalty.o output.o
	mkdir -p obj bin
	$(CXX) -o sa.x $(LIBPATH) atom.o main.o ewald.o lj12.o bv.o bvv.o image.o stress.o readion.o readpara.o sa.o simann.o penalty.o output.o
	mv *.o bin/
%.o: %.c $(DEPS)
	$(CXX)  $(CXXFLAGS) -c -o $@ $^
clean:
	rm -rf *.o obj bin
