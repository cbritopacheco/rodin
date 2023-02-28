#include <Corrade/TestSuite/Tester.h>

namespace Rodin::Tests
{
  struct Sum : TestSuite::Tester {
    explicit MyTest();

    void commutativity();
    void associativity();
    void pi();
    void sin();
    void bigEndian();

    void prepend1kItemsVector();
    void prepend1kItemsList();
};
