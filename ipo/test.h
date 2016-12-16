#include <string>
#include <cstdio>
using namespace std;

namespace ipo{
class Foo {
	public:
		Foo();
		std::string printFoo();

		~Foo(){
			printf("Destructing\n");
		}
};
}
