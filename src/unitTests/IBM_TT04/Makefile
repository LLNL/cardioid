HEART_FILES = IBM_TT04.cc IBM_TT04.hh IBM_TT04_LUT.cc IBM_TT04_LUT.hh

BGPHOME = /bgsys/drivers/ppcfloor

INC = -I. 
#-I$(BGPHOME)/comm/include

LIB = 
#-L$(BGPHOME)/comm/lib -lmpich.cnk -ldcmfcoll.cnk -ldcmf.cnk
#L2 = -L$(BGPHOME)/runtime/SPI -lSPI.cna -lpthread -lrt



# Compiler to be used - blrts_xlC for BGD/BGW
#COMP=blrts_xlC
COMP=g++

# Compile flags
CFLAGS = -g

SRC=SingleCell.cc
OBJS= ./obj/SingleCell.o ./obj/IBM_TT04.o ./obj/IBM_TT04_LUT.o
#OBJS= ./obj/SingleCell.o ./obj/IBM_TT04.o
TARGET= SingleCell



##########################################
#               main program             #
##########################################

./$(TARGET) : $(HEART_FILES) $(OBJS)
	$(COMP) $(CFLAGS) $(INC) $(OBJS) $(LIB) -o $(TARGET)

$(HEART_FILES):
	../tools/mkLinks_heart.sh $@

links: $(HEART_FILES)
	@echo "links made"

clean:
	rm -f $(OBJS) $(TARGET)

cleanall:
	rm -f -r $(OBJS) $(TARGET)

###########################################
#               subroutines               #
###########################################
./obj: 
	mkdir -p ./obj

./obj/SingleCell.o : $(SRC) ./obj
	$(COMP) $(CFLAGS) $(INC) $(LIB) -c $(SRC) -o ./obj/SingleCell.o

./obj/IBM_TT04.o: ./IBM_TT04.cc ./IBM_TT04.hh  ./obj
	$(COMP) $(CFLAGS) $(DEFINES) $(INC) $(LIB) -c ./IBM_TT04.cc -o ./obj/IBM_TT04.o


./obj/IBM_TT04_LUT.o: ./IBM_TT04_LUT.cc ./IBM_TT04_LUT.hh  ./obj
	$(COMP) $(CFLAGS) $(DEFINES) $(INC) $(LIB) -c ./IBM_TT04_LUT.cc -o ./obj/IBM_TT04_LUT.o


distclean:
	rm -rf obj $(HEART_FILES) $(TARGET)
