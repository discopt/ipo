#include <string>
#include <cstdio>
using namespace std;

class Foo2 {
	public:
		Foo2();
		std::string printFoo2();

		~Foo2(){
			printf("Destructing\n");
		}
};
