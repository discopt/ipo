#include <string>
#include <cstdio>
using namespace std;

class Foo {
	public:
		Foo();
		std::string printFoo();

		~Foo(){
			printf("Destructing\n");
		}
};
