CXXARGS=-std=c++17 -O3
OUTPUT=a.out

all: $(OUTPUT)

$(OUTPUT): logodds.cpp
	$(CXX) -o $(OUTPUT) logodds.cpp $(CXXARGS)

clean:
	@rm -f $(OUTPUT)

