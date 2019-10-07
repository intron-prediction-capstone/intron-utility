CXXARGS=-std=c++17 -O3
OUTPUT=a.out
INPUT=$(wildcard *.cpp)

all: $(OUTPUT)

$(OUTPUT): $(INPUT)
	$(CXX) -o $(OUTPUT) $(INPUT) $(CXXARGS)

clean:
	@rm -f $(OUTPUT)

