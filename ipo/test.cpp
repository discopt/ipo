#include "test.h"

namespace ipo {
Foo::Foo(){
    printf("Create Foo\n");
}

std::string Foo::printFoo(){
	return "I am a Test!";
}
}
