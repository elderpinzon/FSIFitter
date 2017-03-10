#include "FSIFitterConfig.hxx"

// Declare the global configuration data table:
const ConfigParams configParams("app_configs.ini");

// The constructor:
ConfigParams::ConfigParams(const std::string & configFileName)
  : param1(getConfigFromFile(configFileName, "param1"))
  , param2(getConfigFromFile(configFileName, "param2"))
    // and so on ...
{
}

// And this function would open the file and lookup the requested parameter key, returning its value:
int getConfigFromFile(const std::string & configFileName, const std::string & key)
{
  return 1;
}
