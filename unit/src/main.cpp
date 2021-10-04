
#include <cppunit/CompilerOutputter.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextTestProgressListener.h>
#include <cppunit/TextTestRunner.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <fstream>
#include <iostream>

using namespace CppUnit;

class CustomTestProgressListener : public TextTestProgressListener {
 public:
  virtual void startTest(Test *test)
  {
    std::cerr << "Running " << test->getName().c_str() << "() .";
    std::cerr.flush();
  }
  virtual void endTest(Test *test)
  {
    std::cerr << std::endl;
    std::cerr.flush();
  }
};

int
main(int argc, char *argv[])
{
  // Switch the standard TITANIA output to a temporary file
  std::fstream output;
  output.open("TITANIA.out", std::ios::binary | std::ios::in);

  std::cin.tie(&output);

  TextTestRunner runner;
  TestFactoryRegistry &registry = TestFactoryRegistry::getRegistry();

  if (argc > 2)
  {
    // the simple lookup below doesn't work for more than one argument
    std::cerr << "ERROR: too many arguments given to " << argv[0] << std::endl;
    return 1;
  }
  // run all tests if none specified on command line
  Test *test_to_run = registry.makeTest();
  if (argc > 1)
    test_to_run = test_to_run->findTest(argv[1]);

  runner.addTest(test_to_run);

  TestResult controller;
  TestResultCollector result;
  controller.addListener(&result);

  CustomTestProgressListener progress;
  controller.addListener(&progress);

  try
  {
    runner.run(controller);
  }
  catch (std::invalid_argument &e)
  {
    // Test path not resolved
    std::cerr << std::endl << "ERROR: " << e.what() << std::endl;
    return 0;
  }
  CompilerOutputter outputter(&result, std::cerr);
  outputter.write();

  // TODO(vsz) add xml output to be parsed by GitLab
  // std::ofstream outStream( "test_results.xml" );
  // XmlOutputter xmloutputter( &result, outStream );
  // xmloutputter.write();

  return result.wasSuccessful() ? 0 : 1;
}
